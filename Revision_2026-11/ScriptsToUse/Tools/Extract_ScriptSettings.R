################################################################################
# Extract_ScriptSettings.R -- record which settings each script copy was run with
################################################################################
# Author:   Jonathon Sewell (with Claude)
# Version:  3
# Created:  2026-10-22
# v1 -> v2: v1 only recognised settings at the left margin, so it missed scripts
#   whose "Variables to Change" block is indented (10.1, 12.1, 12.1.1). Now also
#   records indented  name <- value  lines when the value is a plain constant
#   (text, number, TRUE/FALSE, NULL, c(...) of those), and stops at the first
#   "#### Read in ..." section header.
# v2 -> v3: that stop rule cut off settings in scripts with an early "Load" section
#   (e.g. 16_Barplot), so it was removed. Each setting now records the section
#   header it sits under (column "section"), so late, unimportant lines can be
#   filtered out instead.
#
# PURPOSE
#   Your figure folders hold many copies of the same script (e.g. 13 copies of
#   one LineGraph script with 10 different contents). Each copy has the
#   settings that were typed in for one particular figure panel, so the copies
#   are the only record of how each panel was made. This tool reads every .R
#   script under a folder and writes down the settings near the top of each one
#   (lines such as  Set_Threshold <- 0.95  or  TreatmentToAnalyze <- c(...) ).
#   It only READS your files; nothing is changed or uploaded.
#
#   Identical copies are listed once ("copies" says where they live).
#
# WHAT IS WRITTEN (into OUTPUT_FOLDER)
#   script_settings_<time>.csv   one row per setting per distinct script:
#       script_name, md5, n_copies, copies, line_no, variable, value
#   Paths are shortened: "/Users/<name>/..." becomes "~/...".
#   No data are read, only the script text.
#
# HOW TO USE
#   1. Fill in USER SETTINGS (or leave ROOT_FOLDER blank to be asked).
#   2. Click "Source".
#   3. Send the CSV to Claude.
#
# HOW IT DECIDES WHAT IS A SETTING
#   Starting at the top of the script, a line counts as a setting if
#     (a) it starts at the left margin with "name <- value" or "name = value", or
#     (b) it is indented, uses "<-", and the value is a plain constant
#         (text, number, TRUE/FALSE, NULL, c(...) of those).
#   Scanning stops at the first function definition or after MAX_LINES lines.
#   Multi-line values are joined; trailing comments are dropped.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
ROOT_FOLDER   <- ""   # folder to scan, including sub-folders ("" = ask)
OUTPUT_FOLDER <- ""   # where to save the CSV ("" = inside ROOT_FOLDER)

SCRIPT_EXTENSIONS <- c("R", "r")    # which files to read
MAX_LINES <- 800                    # stop looking for settings after this line
SKIP_FOLDER_PATTERN <- "^(\\.git|\\.Rproj\\.user|renv)$"
#### END OF USER SETTINGS ######################################################


#### Set up ####################################################################
if (!requireNamespace("svDialogs", quietly = TRUE)) {
  stop("Please install the 'svDialogs' package: install.packages(\"svDialogs\")",
       call. = FALSE)
}
if (!nzchar(ROOT_FOLDER)) {
  ROOT_FOLDER <- svDialogs::dlg_dir(default = getwd(),
                                    title = "Select the folder with the scripts")$res
  if (length(ROOT_FOLDER) == 0 || !nzchar(ROOT_FOLDER)) stop("No folder selected.", call. = FALSE)
}
if (!dir.exists(ROOT_FOLDER)) stop("ROOT_FOLDER not found: ", ROOT_FOLDER, call. = FALSE)
ROOT_FOLDER <- normalizePath(ROOT_FOLDER, mustWork = TRUE)
if (!nzchar(OUTPUT_FOLDER)) OUTPUT_FOLDER <- ROOT_FOLDER
if (!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER, recursive = TRUE)
tag <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")


#### Helpers ###################################################################
# For one line of R code: bracket depth change, and the code before any comment.
# Characters inside quotes are ignored.
scan_line <- function(x) {
  chars <- strsplit(x, "")[[1]]
  in_q <- ""
  depth <- 0
  end <- length(chars)
  i <- 1
  while (i <= length(chars)) {
    ch <- chars[i]
    if (nzchar(in_q)) {
      if (ch == "\\") i <- i + 1
      else if (ch == in_q) in_q <- ""
    } else {
      if (ch %in% c("\"", "'", "`")) in_q <- ch
      else if (ch == "#") { end <- i - 1; break }
      else if (ch %in% c("(", "[", "{")) depth <- depth + 1
      else if (ch %in% c(")", "]", "}")) depth <- depth - 1
    }
    i <- i + 1
  }
  list(depth = depth, code = trimws(paste(chars[seq_len(max(end, 0))], collapse = "")))
}

# TRUE when a value is a plain constant (no function calls other than c())
is_constant <- function(rhs) {
  x <- rhs
  x <- gsub("\"[^\"]*\"|'[^']*'", "S", x)
  x <- gsub("\\b(TRUE|FALSE|NULL|NA|Inf|T|F)\\b", "", x)
  x <- gsub("[A-Za-z_.][A-Za-z0-9_.-]*\\s*=\\s*", "", x)
  x <- gsub("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?L?", "", x)
  x <- gsub("\\bc\\(", "", x)
  x <- gsub("[(),:[:space:]S]", "", x)
  nchar(x) == 0
}

shorten_path <- function(x) gsub("/Users/[^/\"' ]+", "~", x)

extract_settings <- function(path) {
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  lines <- utils::head(lines, MAX_LINES)
  out <- list()
  section <- ""
  i <- 1
  while (i <= length(lines)) {
    ln <- lines[i]
    # stop at the first function definition at the left margin
    if (grepl("^[A-Za-z_.][A-Za-z0-9_.]*\\s*(<-|=)\\s*function\\s*\\(", ln, perl = TRUE)) break
    if (grepl("^#{3,}\\s*\\S", ln)) {
      section <- trimws(gsub("^#+|#+\\s*$", "", ln))
    }
    at_margin <- grepl("^[A-Za-z_.][A-Za-z0-9_.$]*\\s*(<-|=)(?!=)", ln, perl = TRUE)
    indented  <- grepl("^\\s+[A-Za-z_.][A-Za-z0-9_.$]*\\s*<-", ln, perl = TRUE)
    if ((at_margin || indented) && !grepl("^\\s*(library|require)\\(", ln)) {
      first <- i
      sc <- scan_line(ln)
      code <- sc$code
      depth <- sc$depth
      while (depth > 0 && i < length(lines)) {
        i <- i + 1
        sc2 <- scan_line(lines[i])
        code <- paste(code, sc2$code)
        depth <- depth + sc2$depth
      }
      parts <- regmatches(code, regexec("^([^<=]+?)\\s*(<-|=)\\s*(.*)$", code, perl = TRUE))[[1]]
      keep_it <- length(parts) == 4 && !grepl("^function\\s*\\(", parts[4]) &&
        (at_margin || is_constant(parts[4]))
      if (isTRUE(keep_it)) {
        out[[length(out) + 1]] <- data.frame(
          line_no  = first,
          indented = !at_margin,
          section  = substr(section, 1, 80),
          variable = trimws(parts[2]),
          value    = shorten_path(substr(gsub("\\s+", " ", parts[4]), 1, 500)),
          stringsAsFactors = FALSE)
      }
    }
    i <- i + 1
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}


#### Find scripts ##############################################################
all_paths <- list.files(ROOT_FOLDER, recursive = TRUE, all.files = TRUE, include.dirs = FALSE)
keep <- vapply(all_paths, function(p) {
  parts <- strsplit(p, "/", fixed = TRUE)[[1]]
  if (!(tools::file_ext(p) %in% SCRIPT_EXTENSIONS)) return(FALSE)
  if (length(parts) > 1 && any(grepl(SKIP_FOLDER_PATTERN, parts[-length(parts)]))) return(FALSE)
  TRUE
}, logical(1))
all_paths <- all_paths[keep]
if (length(all_paths) == 0) stop("No R scripts found in ", ROOT_FOLDER, call. = FALSE)

md5 <- unname(tools::md5sum(file.path(ROOT_FOLDER, all_paths)))
cat("Found ", length(all_paths), " script files (", length(unique(md5)),
    " distinct contents). Reading...\n", sep = "")


#### Read each distinct script once ############################################
rows <- list()
for (m in unique(md5)) {
  idx <- which(md5 == m)
  copies <- all_paths[idx]
  settings <- extract_settings(file.path(ROOT_FOLDER, copies[1]))
  if (is.null(settings)) next
  rows[[length(rows) + 1]] <- data.frame(
    script_name = basename(copies[1]),
    md5 = m,
    n_copies = length(copies),
    copies = paste(copies, collapse = " ; "),
    settings,
    stringsAsFactors = FALSE
  )
}
result <- do.call(rbind, rows)
result <- result[order(result$script_name, result$md5, result$line_no), ]
row.names(result) <- NULL

out_file <- file.path(OUTPUT_FOLDER, paste0("script_settings_", tag, ".csv"))
write.csv(result, out_file, row.names = FALSE)

cat("\nDistinct scripts with settings: ", length(unique(result$md5)), "\n", sep = "")
cat("Setting lines recorded: ", nrow(result), "\n", sep = "")
cat("Written: ", out_file, "\n", sep = "")
