################################################################################
# Make_FileManifest.R -- list and fingerprint every file in your figure folders
################################################################################
# Author:   Jonathon Sewell (with Claude)
# Version:  1
# Created:  2026-10-22
#
# PURPOSE
#   Before reorganizing the analysis scripts, we need to know exactly what the
#   published analysis folders contain. This tool walks through a folder (e.g.
#   your "Figures_v12,13" folder, including all sub-folders) and writes a small
#   table describing every file. It does NOT change, move or upload anything.
#
#   The table answers questions such as:
#     * Are the cluster_R1_assigns.csv files in Figures 1, 2 and 3 identical?
#       (same MD5 = byte-for-byte the same file)
#     * Does Figure 4 hold a *different* cluster_R1_assigns.csv than Figure 1?
#     * Which script versions are sitting in which folder?
#     * How many rows/columns does each CSV have, and what are its columns?
#   Later, the same table is the "answer key" for the golden-file comparison:
#   after the scripts are reorganized, the new outputs are checked against it.
#
# WHAT IS WRITTEN (into OUTPUT_FOLDER)
#   manifest_<time>.csv           one row per file (see columns below)
#   manifest_duplicates_<time>.csv  files with identical content in >1 place
#   manifest_summary_<time>.txt   plain-language summary (read this first)
#   manifest_record_<time>.json   settings, R version, time (run record)
#   Only RELATIVE paths are written, so the files can be shared without
#   revealing your folder names above the chosen folder.
#
# MANIFEST COLUMNS
#   path, top_folder, file_name, ext, size_bytes, modified, md5,
#   script_version   (first "Version" line found in the header of .R scripts)
#   For .csv files also:
#   n_rows, n_cols, columns (first 40 names, separated by |),
#   n_numeric_cols, numeric_sum (a fingerprint of all numbers, 10 digits),
#   n_missing        (only for files up to DEEP_READ_MAX_MB; larger files only
#                     get the row/column counts)
#
# HOW TO USE
#   1. Fill in USER SETTINGS (or leave ROOT_FOLDER blank to be asked).
#   2. Click "Source". Big folders can take a few minutes (MD5 of large CSVs).
#   3. Send the manifest_*.csv and manifest_summary_*.txt to Claude. They are
#      small and contain no data values (only totals).
################################################################################


#### USER SETTINGS -- edit this section only ###################################
ROOT_FOLDER   <- ""   # folder to scan, including sub-folders ("" = ask)
OUTPUT_FOLDER <- ""   # where to save the manifest ("" = inside ROOT_FOLDER)

# CSV files up to this size (in MB) are read in full to count missing values and
# fingerprint the numbers. Larger CSVs get row/column counts only.
DEEP_READ_MAX_MB <- 150

# Files and folders to ignore
SKIP_FILE_PATTERN   <- "^(\\.DS_Store|\\.Rhistory|\\.Rapp\\.history|~\\$.*)$"
SKIP_FOLDER_PATTERN <- "^(\\.git|\\.Rproj\\.user|renv)$"

# Which key CSV names to single out in the summary (cross-folder comparison)
KEY_FILES <- c("cluster_R1_assigns.csv", "cluster_R1_groups.csv",
               "cluster_R2_assigns.csv", "cluster_R2_groups.csv",
               "cluster_R2_assigns_merged.csv",
               "cluster_R1_merge.csv", "cluster_R2_merge.csv",
               "expression_matrix_analysis.csv", "expression_matrix_other.csv",
               "filenums.csv", "metadata.csv", "panel.csv",
               "umap_xy_R1.csv", "umap_xy_R2.csv", "umap_xyz_3D_R1.csv")
#### END OF USER SETTINGS ######################################################


#### Set up ####################################################################
if (!requireNamespace("svDialogs", quietly = TRUE)) {
  stop("Please install the 'svDialogs' package: install.packages(\"svDialogs\")",
       call. = FALSE)
}
have_fread <- requireNamespace("data.table", quietly = TRUE)
have_json  <- requireNamespace("jsonlite", quietly = TRUE)

if (!nzchar(ROOT_FOLDER)) {
  ROOT_FOLDER <- svDialogs::dlg_dir(default = getwd(),
                                    title = "Select the folder to scan")$res
  if (length(ROOT_FOLDER) == 0 || !nzchar(ROOT_FOLDER)) stop("No folder selected.", call. = FALSE)
}
if (!dir.exists(ROOT_FOLDER)) stop("ROOT_FOLDER not found: ", ROOT_FOLDER, call. = FALSE)
ROOT_FOLDER <- normalizePath(ROOT_FOLDER, mustWork = TRUE)
if (!nzchar(OUTPUT_FOLDER)) OUTPUT_FOLDER <- ROOT_FOLDER
if (!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER, recursive = TRUE)

run_time <- Sys.time()
tag <- format(run_time, "%Y-%m-%d_%H-%M-%S")

report <- character(0)
say <- function(...) {
  line <- paste0(...)
  report <<- c(report, line)
  cat(line, "\n", sep = "")
}


#### Find files ################################################################
all_paths <- list.files(ROOT_FOLDER, recursive = TRUE, all.files = TRUE,
                        full.names = FALSE, include.dirs = FALSE)
keep <- vapply(all_paths, function(p) {
  parts <- strsplit(p, "/", fixed = TRUE)[[1]]
  if (grepl(SKIP_FILE_PATTERN, parts[length(parts)])) return(FALSE)
  if (length(parts) > 1 &&
      any(grepl(SKIP_FOLDER_PATTERN, parts[-length(parts)]))) return(FALSE)
  TRUE
}, logical(1))
all_paths <- all_paths[keep]
if (length(all_paths) == 0) stop("No files found in ", ROOT_FOLDER, call. = FALSE)
cat("Found ", length(all_paths), " files. Fingerprinting...\n", sep = "")


#### Per-file helpers ##########################################################
# Count lines without loading the whole file
count_lines <- function(path, chunk = 100000L) {
  con <- file(path, "r")
  on.exit(close(con))
  n <- 0
  repeat {
    x <- readLines(con, n = chunk, warn = FALSE)
    if (length(x) == 0) break
    n <- n + length(x)
  }
  n
}

# First "Version" line in the top of a script
script_version <- function(path) {
  x <- tryCatch(readLines(path, n = 40, warn = FALSE), error = function(e) character(0))
  hit <- grep("^\\s*#+\\s*(Script\\s+)?Version\\s*[:=]", x, ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  trimws(sub("^\\s*#+\\s*", "", hit[1]))
}

csv_info <- function(path, size_bytes) {
  out <- list(n_rows = NA_real_, n_cols = NA_integer_, columns = NA_character_,
              n_numeric_cols = NA_integer_, numeric_sum = NA_real_, n_missing = NA_real_)
  header <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character(0))
  if (length(header) == 0) return(out)
  cols <- tryCatch(names(read.csv(text = header, check.names = FALSE, nrows = 0)),
                   error = function(e) character(0))
  out$n_cols  <- length(cols)
  out$columns <- paste(utils::head(cols, 40), collapse = "|")
  out$n_rows  <- max(count_lines(path) - 1, 0)

  if (size_bytes <= DEEP_READ_MAX_MB * 1024^2) {
    df <- tryCatch({
      if (have_fread) data.table::fread(path, data.table = FALSE, check.names = FALSE)
      else read.csv(path, check.names = FALSE, comment.char = "")
    }, error = function(e) NULL)
    if (!is.null(df)) {
      is_num <- vapply(df, is.numeric, logical(1))
      out$n_numeric_cols <- sum(is_num)
      out$n_missing      <- sum(is.na(df))
      if (any(is_num)) {
        out$numeric_sum <- signif(sum(vapply(df[is_num], function(v) sum(as.numeric(v), na.rm = TRUE),
                                             numeric(1))), 10)
      }
    }
  }
  out
}


#### Build the manifest ########################################################
rows <- vector("list", length(all_paths))
for (i in seq_along(all_paths)) {
  p    <- all_paths[i]
  full <- file.path(ROOT_FOLDER, p)
  info <- file.info(full)
  ext  <- tolower(tools::file_ext(p))
  parts <- strsplit(p, "/", fixed = TRUE)[[1]]
  row <- data.frame(
    path        = p,
    top_folder  = if (length(parts) > 1) parts[1] else "(root)",
    file_name   = parts[length(parts)],
    ext         = ext,
    size_bytes  = info$size,
    modified    = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    md5         = unname(tools::md5sum(full)),
    script_version = if (ext %in% c("r")) script_version(full) else NA_character_,
    n_rows = NA_real_, n_cols = NA_integer_, columns = NA_character_,
    n_numeric_cols = NA_integer_, numeric_sum = NA_real_, n_missing = NA_real_,
    stringsAsFactors = FALSE
  )
  if (ext == "csv") {
    ci <- csv_info(full, info$size)
    for (nm in names(ci)) row[[nm]] <- ci[[nm]]
  }
  rows[[i]] <- row
  if (i %% 25 == 0) cat("  ", i, " / ", length(all_paths), "\n", sep = "")
}
manifest <- do.call(rbind, rows)
manifest <- manifest[order(manifest$path), ]
row.names(manifest) <- NULL

manifest_file <- file.path(OUTPUT_FOLDER, paste0("manifest_", tag, ".csv"))
write.csv(manifest, manifest_file, row.names = FALSE)


#### Duplicates (same content in more than one place) ##########################
dup_md5 <- names(which(table(manifest$md5) > 1))
dup_tbl <- do.call(rbind, lapply(dup_md5, function(m) {
  sub <- manifest[manifest$md5 == m, ]
  data.frame(md5 = m, n_copies = nrow(sub),
             file_names = paste(unique(sub$file_name), collapse = " | "),
             size_bytes = sub$size_bytes[1],
             locations = paste(sub$path, collapse = " ; "),
             stringsAsFactors = FALSE)
}))
if (is.null(dup_tbl)) dup_tbl <- data.frame(md5 = character(), n_copies = integer(),
                                            file_names = character(), size_bytes = numeric(),
                                            locations = character())
dup_file <- file.path(OUTPUT_FOLDER, paste0("manifest_duplicates_", tag, ".csv"))
write.csv(dup_tbl, dup_file, row.names = FALSE)


#### Plain-language summary ####################################################
say("FILE MANIFEST SUMMARY  ", format(run_time, "%Y-%m-%d %H:%M:%S"))
say("Scanned folder: ", basename(ROOT_FOLDER), "   (", nrow(manifest), " files)")
say(strrep("-", 72))

say("")
say("FILES PER TOP-LEVEL FOLDER")
per_top <- table(manifest$top_folder)
for (nm in names(per_top)) say("  ", format(per_top[[nm]], width = 4), "  ", nm)

# Key data files: same name, same content?
say("")
say("KEY DATA FILES -- are the copies identical across folders?")
say("  (copies sharing a letter have identical content; different letters = the")
say("   file was changed or regenerated between folders)")
for (kf in KEY_FILES) {
  sub <- manifest[manifest$file_name == kf, ]
  if (nrow(sub) == 0) next
  groups <- unique(sub$md5)
  say("")
  say("  ", kf, "  (", nrow(sub), " copies, ", length(groups), " distinct version",
      ifelse(length(groups) == 1, "", "s"), ")")
  for (g in seq_along(groups)) {
    s2 <- sub[sub$md5 == groups[g], ]
    dims <- if (!is.na(s2$n_rows[1])) paste0(format(s2$n_rows[1], big.mark = ","), " rows x ",
                                             s2$n_cols[1], " cols") else ""
    say("    ", LETTERS[g], ": ", dims, "   in: ",
        paste(unique(s2$top_folder), collapse = "; "))
  }
}

# Script versions
scripts <- manifest[manifest$ext == "r", ]
if (nrow(scripts) > 0) {
  say("")
  say("SCRIPTS (distinct file names; where each copy lives)")
  for (nm in sort(unique(scripts$file_name))) {
    s2 <- scripts[scripts$file_name == nm, ]
    n_ver <- length(unique(s2$md5))
    say("  ", nm, "   [", nrow(s2), " cop", ifelse(nrow(s2) == 1, "y", "ies"), ", ",
        n_ver, " distinct content", ifelse(n_ver == 1, "", "s"), "]",
        if (!is.na(s2$script_version[1])) paste0("   ", s2$script_version[1]) else "")
  }
}

# Output-like CSVs (everything that is not a key data file)
other_csv <- manifest[manifest$ext == "csv" & !(manifest$file_name %in% KEY_FILES), ]
say("")
say("OTHER CSV FILES (likely script outputs / candidate golden files): ", nrow(other_csv))
if (nrow(other_csv) > 0) {
  for (tf in unique(other_csv$top_folder)) {
    s2 <- other_csv[other_csv$top_folder == tf, ]
    say("  ", tf, ": ", nrow(s2), " files")
  }
}

say("")
say("Identical-content groups (copies in >1 place): ", nrow(dup_tbl),
    "   -> ", basename(dup_file))
say("Total size of files that are exact duplicates of another: ",
    format(round(sum(manifest$size_bytes[duplicated(manifest$md5)]) / 1024^2, 1), nsmall = 1),
    " MB")
say("")
say("Files written to: ", OUTPUT_FOLDER)

summary_file <- file.path(OUTPUT_FOLDER, paste0("manifest_summary_", tag, ".txt"))
writeLines(report, summary_file)


#### Run record ################################################################
if (have_json) {
  rec <- list(
    Module = "Make_FileManifest", Script_Version = 1,
    Run_At = format(run_time, "%Y-%m-%d %H:%M:%S"),
    Scanned_Folder_Name = basename(ROOT_FOLDER),
    Files = nrow(manifest),
    Settings = list(DEEP_READ_MAX_MB = DEEP_READ_MAX_MB,
                    SKIP_FILE_PATTERN = SKIP_FILE_PATTERN,
                    SKIP_FOLDER_PATTERN = SKIP_FOLDER_PATTERN),
    Environment = list(R = R.version.string, platform = R.version$platform,
                       data.table = if (have_fread) as.character(utils::packageVersion("data.table")) else NA,
                       jsonlite = as.character(utils::packageVersion("jsonlite")))
  )
  jsonlite::write_json(rec, file.path(OUTPUT_FOLDER, paste0("manifest_record_", tag, ".json")),
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
}
