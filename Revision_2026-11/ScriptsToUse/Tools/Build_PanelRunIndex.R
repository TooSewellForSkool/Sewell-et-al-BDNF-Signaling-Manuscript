################################################################################
# Build_PanelRunIndex.R -- fill the panel-to-run index from your run records
################################################################################
# Layer:    Tools (reads run_parameters.json files and one CSV; writes new CSVs; changes nothing else)
# Author:   Jonathon Sewell (with Claude)
# Version:  1
# Created:  2026-10-29
#
# PURPOSE
#   Panel_Run_Index.csv says, for every figure panel, which script chain and which run folder made it.
#   After you regenerate the figures, type the run folder name into the column "run_folder(s)" (one
#   or more folder names, separated by ";"). This script then looks each folder up, reads its
#   run_parameters.json, and writes the settings and script versions into the index. Nothing is
#   guessed: a folder that cannot be found is reported and left empty.
#
# HOW TO USE
#   1. Set ROOT_FOLDER (the project folder that holds your run folders) or leave "" to pick it.
#   2. Set INDEX_FILE (leave "" to pick Panel_Run_Index.csv).
#   3. Click Source. Two files are written next to the index:
#        Run_Records_Scanned.csv       one row per run folder and module (all settings, as text)
#        Panel_Run_Index_FILLED.csv    the index with run_record_file(s), script_versions, output_md5 filled
#
# NOTE: written without being able to run R here. If it stops with an error, send me the message.
################################################################################

#### USER SETTINGS -- edit this section only ###################################
ROOT_FOLDER <- ""     # folder that contains the run folders (searched in all sub-folders); "" = choose
INDEX_FILE  <- ""     # Panel_Run_Index.csv; "" = choose

#### Setup #####################################################################
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Please install the package jsonlite first.", call. = FALSE)
# Pop-ups that work on Mac and Windows without XQuartz (RStudio's own picker, else the svDialogs package)
pick <- function(what, caption) {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(if (what == "dir") rstudioapi::selectDirectory(caption = caption) else rstudioapi::selectFile(caption = caption))
  }
  if (!requireNamespace("svDialogs", quietly = TRUE)) stop("Set the path in USER SETTINGS (no pop-up available).", call. = FALSE)
  res <- if (what == "dir") svDialogs::dlg_dir(title = caption)$res else svDialogs::dlg_open(title = caption)$res
  if (length(res) == 0 || !nzchar(res)) NULL else res
}
if (ROOT_FOLDER == "") ROOT_FOLDER <- pick("dir",  "Select the project folder that holds your run folders")
if (INDEX_FILE  == "") INDEX_FILE  <- pick("file", "Select Panel_Run_Index.csv")
if (is.null(ROOT_FOLDER) || is.null(INDEX_FILE)) stop("No folder or file selected. Run the script again.", call. = FALSE)
shorten <- function(x) gsub("/Users/[^/\"' ]+|[A-Za-z]:[/\\\\]Users[/\\\\][^/\\\\\"' ]+", "~", x)

#### Scan run records ##########################################################
files <- list.files(ROOT_FOLDER, pattern = "^run_parameters\\.json$", recursive = TRUE, full.names = TRUE)
cat("[INFO] ", length(files), " run_parameters.json files found\n", sep = "")
rows <- list()
for (f in files) {
  rec <- tryCatch(jsonlite::read_json(f, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(rec)) { cat("[WARN] could not read ", f, "\n", sep = ""); next }
  folder <- basename(dirname(f))
  info <- rec[["Folder_Info"]]
  for (m in setdiff(names(rec), "Folder_Info")) {
    sec <- rec[[m]]
    ver <- if (is.list(sec)) { v <- sec[["Script_Version"]]; if (is.null(v)) sec[["SCRIPT_VERSION"]] else v } else NULL
    md5 <- if (is.list(sec)) sec[["Output_MD5"]] else NULL
    rows[[length(rows) + 1]] <- data.frame(
      run_folder = folder, module = m,
      last_updated = if (is.null(info[["Last_Updated"]])) "" else info[["Last_Updated"]],
      script_version = if (is.null(ver)) "" else paste(unlist(ver), collapse = " "),
      output_md5 = if (is.null(md5)) "" else paste(names(md5), unlist(md5), sep = "=", collapse = "; "),
      settings = shorten(as.character(jsonlite::toJSON(sec, auto_unbox = TRUE))),
      stringsAsFactors = FALSE)
  }
}
scan <- if (length(rows)) do.call(rbind, rows) else data.frame()
out_dir <- dirname(INDEX_FILE)
utils::write.csv(scan, file.path(out_dir, "Run_Records_Scanned.csv"), row.names = FALSE)

#### Fill the index ############################################################
idx <- utils::read.csv(INDEX_FILE, check.names = FALSE, stringsAsFactors = FALSE)
col_run <- "run_folder(s)"
if (!(col_run %in% names(idx))) stop("The index has no column '", col_run, "'.", call. = FALSE)
for (i in seq_len(nrow(idx))) {
  fl <- trimws(strsplit(idx[[col_run]][i], ";")[[1]]); fl <- fl[nzchar(fl)]
  if (!length(fl)) next
  rec_files <- character(); vers <- character(); md5s <- character()
  for (fo in fl) {
    hit <- scan[scan$run_folder == fo, , drop = FALSE]
    if (nrow(hit) == 0) { cat("[WARN] ", idx$paper_panel[i], ": run folder not found: ", fo, "\n", sep = ""); next }
    rec_files <- c(rec_files, file.path(fo, "run_parameters.json"))
    vers <- c(vers, paste0(hit$module, " ", hit$script_version, collapse = ", "))
    md5s <- c(md5s, hit$output_md5[nzchar(hit$output_md5)])
  }
  idx[["run_record_file(s)"]][i] <- paste(rec_files, collapse = "; ")
  idx[["script_versions"]][i]    <- paste(vers, collapse = " | ")
  idx[["output_md5"]][i]         <- paste(md5s, collapse = " | ")
}
f_out <- file.path(out_dir, "Panel_Run_Index_FILLED.csv")
utils::write.csv(idx, f_out, row.names = FALSE)
cat("[DONE] Written to ", out_dir, ":\n  Run_Records_Scanned.csv\n  Panel_Run_Index_FILLED.csv\n", sep = "")
