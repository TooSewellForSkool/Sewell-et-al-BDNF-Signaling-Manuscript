################################################################################
# 02_RunBatch.R
################################################################################
# Layer:    DataCleanup / BatchCorrection -- OPTIONAL step (only if you will do
#           batch correction). Run after 01.3_Update_ChannelsForBatchCorr.
# Author:   Jonathon Sewell
# Created:  2026-07-04
# Version:  6   (input version: 02_RunBatch_v5)
# Edited:   2026-10-21
#
# Purpose:
# * Correct for batch effects across barcode sets stained with the same antibody
#   cocktail, using the anchor (universal) sample in every batch.
# * Asks for the settings (method, keywords, ...), checks that your file names fit
#   them, runs 03_BatchAdjust, and saves the corrected FCS files, a matching
#   metadata.csv / panel.csv and a run record in a NEW folder next to your FCS
#   files:   <date>_<time>_BatchCorrection_<method>
#
# Origin:
# * RunBatch > 02_RunBatch.R
#     - CytofBatchAdjust: Schuyler RP, et al. "Minimizing Batch Effects in Mass Cytometry Data." Front Immunol 2019;10:2367
#     - The copy this file came from credited the Zunder Lab, University of Virginia (kept as found; not verified)
#     - Script to input parameters for batch adjustment by 03_BatchAdjust.R
#     - Schuyler, et al. 2019 (https://10.3389/fimmu.2019.02367)
#     - https://github.com/CUHIMSR/CytofBatchAdjust
#
# ==============================================================================
# CHANGE LOG: v5 -> v6
# ==============================================================================
# Convention:  [RESULT]   = can change the corrected files or metadata
#              [BEHAVIOR] = changes how the script behaves, not the data
#              [SAME]     = output unchanged
#
# A. Needs your attention
#  A1 [BEHAVIOR] Output folder. v5 wrote <date>_BatchCorrection_<method> inside a
#     temporary link (/tmp/...) to your FCS folder, with a date only (a second run
#     on the same day with the same method wrote into the same folder). v6 makes
#     <date>_<time>_BatchCorrection_<method> as a sibling of your other stage
#     folders, directly inside the FCS folder. Its run_parameters.json starts from
#     a copy of the FCS folder's record (full history).
#  A2 [BEHAVIOR] /tmp symlink removed (it does not exist on Windows). 03_BatchAdjust
#     v2 handles folder names with spaces itself.
#  A3 [BEHAVIOR] New CHECKS before anything is corrected (v5 failed part-way, or
#     silently skipped files, if names did not fit):
#       - every FCS file in metadata.csv exists in the folder;
#       - the batch numbers found from the anchor files, including batch 1;
#       - exactly ONE anchor file per batch;
#       - every FCS file belongs to exactly one batch. 03_BatchAdjust only writes
#         files it can assign to a batch, so an unmatched file would be missing
#         from the corrected folder and the file numbers would no longer line up.
#       - every channel in ChannelsToAdjust.txt exists in panel.csv;
#       - after the run, every expected corrected file exists.
#     You see a summary table and confirm before the correction starts.
#  A4 [BEHAVIOR] "Has 01.2 been run?" is no longer a question. It is read from
#     run_parameters.json. "Did you modify ChannelsToAdjust.txt?" is asked only if
#     the file differs from what 01.3 saved (or 01.3 was not run).
#  A5 [RESULT, bug fix] Corrected-file names in the new metadata.csv: v5 used
#     sub(".fcs", ...), where "." matches any character, so a name containing
#     "fcs" earlier than the extension would be altered. v6 only changes the end of
#     the name. Names are identical for your data.
#  A6 [BEHAVIOR] Method "Quantile" is turned into "quantile". 03_BatchAdjust
#     accepts "Quantile" at the start but then applies the scaling branch, which
#     does not match the quantile maps (an error).
#  A7 [BEHAVIOR] Default batch keyword: v5 removed the digits 1-9 from each
#     BarcodeSet and offered the list of results as the default. v6 uses the text
#     all BarcodeSet values have in common (Sample001...Sample006 -> "Sample00",
#     same as v5). Only a single value can be used by 03_BatchAdjust.
#  A8 [BEHAVIOR] 03_BatchAdjust is found by name (03_BatchAdjust.R, else the highest
#     version number among 03_BatchAdjust*.R, ignoring README files and other file types). v5 took the
#     last file alphabetically containing "BatchAdjust", which would pick a README.
#  A9 [BEHAVIOR] The pop-up defaults are the settings from the last batch
#     correction in this folder's run record, if there was one.
#
# B. Efficiency / readability / reproducibility
#  B1 [SAME] The settings you answer in pop-ups are unchanged: Method, Batch
#     Keyword, Anchor Keyword, Transformation, Plot Diagnostics, Channels File.
#  B2 Shared helper file (00_Helpers.R); run record written with record_params().
#  B3 v5 keys in the run record keep their names (Method, BatchKeyword,
#     AnchorKeyword, ChannelsFile, Transformation, PlotDiagnostics, BaseDir,
#     OutputDir, AddExtension); BaseDir is now your FCS folder. New: Corrected_
#     FCS_Folder, Channels, Batches_Found, Anchor_Files, N_Files_Corrected, the
#     03_BatchAdjust version used, checksums, session information.
#  B4 Metadata and panel for the corrected files are written with the same
#     columns as the originals; ChannelsToAdjust.txt is copied in for the record.
#  B5 Removed: dead commented-out code, duplicated yes/no blocks, the local copy
#     of record_params (now in the helpers), setwd() into the output folder.
#
# C. Left as in v5 on purpose (flagged to you, not applied)
#  * The anchor keyword default is a guess from file names containing "niv"
#    (Universal). Check it in the pop-up.
#  * The reference batch is batch 1 (set by 03_BatchAdjust); to change it, rename.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
# (none: everything is chosen in the pop-ups)

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "02_RunBatch.R"
SCRIPT_VERSION           <- 6
MODULE_NAME              <- "02_RunBatch"
PREREQUISITE_MODULES     <- c("01_Get_Metadata", "01.2_Update_Metadata")
REQUIRED_HELPERS_VERSION <- 2
CHANNELS_FILENAME        <- "ChannelsToAdjust.txt"
BATCH_SCRIPT_PATTERN     <- "^03_BatchAdjust.*\\.R$"



#### Load shared helper functions (do not edit this block) #####################
# Looks for ScriptsToUse/Helpers/00_Helpers.R starting at this script's folder
# and moving up to 4 folders.
local({
  start_dir <- tryCatch({
    p <- rstudioapi::getSourceEditorContext()$path
    if (is.character(p) && nzchar(p)) dirname(p) else getwd()
  }, error = function(e) getwd())
  search_dir   <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  helpers_file <- NA_character_
  for (i in 1:5) {
    candidates <- c(file.path(search_dir, "Helpers", "00_Helpers.R"),
                    file.path(search_dir, "00_Helpers.R"))
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates) > 0) { helpers_file <- candidates[1]; break }
    search_dir <- dirname(search_dir)
  }
  if (is.na(helpers_file)) {
    stop("Could not find 00_Helpers.R.\n",
         "It should be in a folder called 'Helpers' inside ScriptsToUse.\n",
         "Looked upward from: ", start_dir,
         "\n(Is this script saved inside your ScriptsToUse folder, and is the ",
         "'rstudioapi' package installed? RStudio needs it to find this script.)",
         call. = FALSE)
  }
  assign("HELPERS_DIR", dirname(helpers_file), envir = .GlobalEnv)
  source(helpers_file, local = FALSE)
})
check_helpers_version(REQUIRED_HELPERS_VERSION)



#### Set up environment ########################################################
pkg_versions <- dc_load_packages(
  cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite"),
  bioc = c("flowCore"))

cat("\n========================================================\n")
cat("[START] Perform batch correction on dataset:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

# Folder of THIS script (03_BatchAdjust must be in the same folder)
script_dir <- tryCatch({
  p <- rstudioapi::getSourceEditorContext()$path
  if (is.character(p) && nzchar(p)) dirname(p) else getwd()
}, error = function(e) getwd())

Original_WD <- get_working_dir()

# Sequence check: 01 and 01.2 must have been run in this folder
run_record <- require_modules(Original_WD, PREREQUISITE_MODULES,
                              how_to_fix = "Run 01_Get_Metadata and 01.2_Update_Metadata first.")
if (!isTRUE(get_record_value(run_record, "01.2_Update_Metadata", "Reviewed", default = FALSE))) {
  stop("01.2_Update_Metadata is recorded but not marked as reviewed. ",
       "Run 01.2_Update_Metadata and click Done.", call. = FALSE)
}



#### Load 03_BatchAdjust (the Zunder lab code) #################################
batch_scripts <- list.files(script_dir, pattern = BATCH_SCRIPT_PATTERN)
if (length(batch_scripts) == 0) {
  stop("Could not find a file named 03_BatchAdjust*.R in:\n  ", script_dir,
       "\n02_RunBatch and 03_BatchAdjust must be in the same folder.", call. = FALSE)
}
ver_num <- suppressWarnings(as.numeric(sub("^.*_v([0-9]+)\\.R$", "\\1", batch_scripts)))
batch_script <- if ("03_BatchAdjust.R" %in% batch_scripts) {
  "03_BatchAdjust.R"                       # the normal file name
} else {
  batch_scripts[if (all(is.na(ver_num))) length(batch_scripts) else which.max(ver_num)]
}
batch_script_path <- file.path(script_dir, batch_script)
cat("[INFO] Using", batch_script, "\n")
source(batch_script_path, local = globalenv())



#### Read metadata and panel from disk #########################################
metadata <- read_metadata(Original_WD)
panel    <- read_panel(Original_WD)

missing_fcs <- metadata$File_Name[!file.exists(file.path(Original_WD, metadata$File_Name))]
if (length(missing_fcs) > 0) {
  stop(length(missing_fcs), " FCS file(s) listed in metadata.csv are not in the folder:\n  ",
       paste(head(missing_fcs, 5), collapse = "\n  "),
       if (length(missing_fcs) > 5) "\n  ..." else "", call. = FALSE)
}



#### Set up user-defined variables #############################################
# Defaults: settings from the last batch correction in this folder, if any
prev <- run_record[[MODULE_NAME]]
default_or <- function(key, fallback) {
  v <- prev[[key]]
  if (is.null(v) || length(v) != 1) fallback else as.character(v)
}

# Text that all BarcodeSet values start with (Sample001..Sample006 -> "Sample00")
common_prefix <- function(x) {
  x <- unique(as.character(x))
  if (length(x) == 0) return("")
  chars <- strsplit(x, "")
  n <- min(lengths(chars))
  keep <- 0
  for (i in seq_len(n)) {
    if (length(unique(vapply(chars, `[`, character(1), i))) == 1) keep <- i else break
  }
  substr(x[1], 1, keep)
}

# Guess the anchor keyword from file names containing "niv" (Universal), as in v5:
# the "_"-separated pieces containing "orm" or "niv", joined by "_"
guess_anchor_keyword <- function(file_names) {
  hits <- file_names[grepl("niv", file_names)]
  if (length(hits) == 0) return("")
  pieces <- unique(unlist(strsplit(hits, split = "_")))
  paste0(pieces[grepl("orm|niv", pieces)], collapse = "_")
}

# Pop-up asking for text; stops if the window is closed
ask_text <- function(message, default = "") {
  res <- svDialogs::dlg_input(message = message, default = default)$res
  if (length(res) != 1 || is.na(res)) {
    stop("The pop-up window was closed without an answer. Please run the script again.",
         call. = FALSE)
  }
  trimws(res)
}

example_anchor <- metadata$File_Name[grep("niv", metadata$File_Name)][1]
method_input        <- ask_text("Method\n\ni.e. quantile, SD, percentile [50p, 90p, 95p, etc.]",
                                default_or("Method", "95p"))
batchKeyword_input  <- ask_text("Batch Keyword\n\ni.e. how to ID barcode sets",
                                default_or("BatchKeyword", common_prefix(metadata$BarcodeSet)))
anchorKeyword_input <- ask_text(paste0("Anchor Keyword\n\ni.e. how to ID the universal sample across barcode sets (batches)\n\nEx: ",
                                       ifelse(is.na(example_anchor), "(no file name contains 'niv')", example_anchor)),
                                default_or("AnchorKeyword", guess_anchor_keyword(metadata$File_Name)))
transformation_input <- as.logical(ask_text("Transformation\n\ni.e. apply an asinh transformation to the data (TRUE/FALSE)",
                                            default_or("Transformation", "FALSE")))
plotDiagnostics_input <- as.logical(ask_text("Plot Diagnostics\n\ni.e. Generate distribution plots for each adjusted channel (TRUE/FALSE)",
                                             default_or("PlotDiagnostics", "TRUE")))
channelsFile_input  <- ask_text("Channels File\n\ni.e. file name containing channels on which to perform batch correction",
                                default_or("ChannelsFile", CHANNELS_FILENAME))

# Check the answers
if (tolower(method_input) %in% c("quantile", "qn")) method_input <- tolower(method_input)
if (tolower(method_input) == "qn") method_input <- "QN"
if (!grepl("^(quantile|QN|SD|sd|hybrid|[0-9]+p)$", method_input)) {
  stop("Method '", method_input, "' is not recognised. Use quantile, SD, or a ",
       "percentile such as 50p or 95p.", call. = FALSE)
}
if (!nzchar(batchKeyword_input) || grepl(" ", batchKeyword_input, fixed = TRUE)) {
  stop("The batch keyword must not be empty or contain spaces.", call. = FALSE)
}
if (!nzchar(anchorKeyword_input)) {
  stop("The anchor keyword is empty. Enter the text that identifies the universal ",
       "sample in the file names.", call. = FALSE)
}
if (is.na(transformation_input) || is.na(plotDiagnostics_input)) {
  stop("Transformation and Plot Diagnostics must be TRUE or FALSE.", call. = FALSE)
}

# Channels file (in the FCS folder)
channels_path <- file.path(Original_WD, channelsFile_input)
if (!file.exists(channels_path)) {
  stop(channelsFile_input, " was not found in the FCS folder. Run ",
       "01.3_Update_ChannelsForBatchCorr, or copy the file there.", call. = FALSE)
}
channels_to_adjust <- unique(trimws(readLines(channels_path, warn = FALSE)))
channels_to_adjust <- channels_to_adjust[nzchar(channels_to_adjust)]
not_in_panel <- setdiff(channels_to_adjust, panel$Name)
if (length(channels_to_adjust) == 0) stop(channelsFile_input, " is empty.", call. = FALSE)
if (length(not_in_panel) > 0) {
  stop("These channels in ", channelsFile_input, " are not in panel.csv (check the ",
       "spelling, e.g. Y89Di):\n  ", paste(not_in_panel, collapse = ", "), call. = FALSE)
}

# Was the channel list confirmed in 01.3?
saved_md5  <- get_record_value(run_record, "01.3_Update_ChannelsForBatchCorr", "Output_MD5")
current_md5 <- unname(tools::md5sum(channels_path))
if (!identical(channelsFile_input, CHANNELS_FILENAME) ||
    is.null(saved_md5) || !identical(as.character(unlist(saved_md5))[1], as.character(current_md5))) {
  if (!ask_yes_no(paste0(channelsFile_input, " was not saved by 01.3_Update_ChannelsForBatchCorr ",
                         "(or was changed afterwards).\n\nDid you modify it for this correction ",
                         "and is it correct (", length(channels_to_adjust), " channels)?"))) {
    stop("Please update ", channelsFile_input, " and run this script again.", call. = FALSE)
  }
}



#### Check that file names fit the settings (nothing is corrected yet) #########
batches_found <- tryCatch(
  suppressWarnings(listBatchesPresent(Original_WD, batchKeyword = batchKeyword_input,
                                      anchorKeyword = anchorKeyword_input)),
  error = function(e) numeric(0))
if (length(batches_found) == 0 || anyNA(batches_found)) {
  stop("No batches could be read from the anchor files.\n",
       "Anchor keyword: '", anchorKeyword_input, "'   Batch keyword: '", batchKeyword_input, "'\n",
       "Each batch needs ONE file whose name has the form  <batch keyword><number>_<anchor keyword>...\n",
       "Example: with batch keyword 'Sample00' and anchor keyword 'normalized_Universal',\n",
       "  Sample001_normalized_Universal_1_...fcs  is the anchor of batch 1.",
       call. = FALSE)
}
if (!(1 %in% batches_found)) {
  stop("Batch 1 (the reference batch) has no anchor file. Batches found: ",
       paste(batches_found, collapse = ", "), ". Check the file naming.", call. = FALSE)
}
if (anyDuplicated(batches_found)) {
  dup <- unique(batches_found[duplicated(batches_found)])
  stop("More than one anchor file was found for batch(es): ", paste(dup, collapse = ", "),
       ". Each batch needs exactly one. Tighten the anchor keyword.", call. = FALSE)
}

batch_tag   <- vapply(batches_found, function(b) sprintf("%s%i_", batchKeyword_input, b), character(1))
file_batch  <- vapply(metadata$File_Name, function(fn) {
  hit <- batches_found[vapply(batch_tag, function(tg) grepl(tg, fn, fixed = TRUE), logical(1))]
  if (length(hit) == 1) hit else NA_real_
}, numeric(1))
unassigned <- metadata$File_Name[is.na(file_batch)]
if (length(unassigned) > 0) {
  stop(length(unassigned), " file(s) do not belong to exactly one batch, so they would be ",
       "left out of the corrected files:\n  ",
       paste(head(unassigned, 5), collapse = "\n  "),
       if (length(unassigned) > 5) "\n  ..." else "",
       "\nFile names must contain '<batch keyword><number>_' (here e.g. '",
       batch_tag[1], "'). Zero-padded numbers need the zeros in the batch keyword ",
       "(Sample001 -> 'Sample00').", call. = FALSE)
}
anchor_files <- vapply(batches_found, function(b) {
  f <- list_fcs_files(Original_WD, name_contains = anchorKeyword_input)
  f <- basename(f)[vapply(basename(f), getBatchNumFromFilename, numeric(1),
                          batchKeyword = batchKeyword_input,
                          anchorKeyword = anchorKeyword_input) == b]
  f[1]
}, character(1))

summary_tbl <- data.frame(Batch = batches_found,
                          Files = as.integer(table(factor(file_batch, levels = batches_found))),
                          Anchor_File = anchor_files, stringsAsFactors = FALSE)
cat("\n[CHECK] Batches found:\n")
print(summary_tbl, row.names = FALSE)

if (!ask_yes_no(paste0(
  "Batch correction will use:\n\n",
  "Method: ", method_input, "\nBatch keyword: ", batchKeyword_input,
  "\nAnchor keyword: ", anchorKeyword_input, "\nTransformation: ", transformation_input,
  "\nPlot diagnostics: ", plotDiagnostics_input,
  "\nChannels: ", length(channels_to_adjust), " (", channelsFile_input, ")",
  "\n\nBatches found: ", paste(batches_found, collapse = ", "),
  "\nFiles to correct: ", nrow(metadata),
  "\n\nThe anchor files are listed in the R console. Continue?"))) {
  stop("Batch correction was cancelled. Nothing was written.", call. = FALSE)
}



#### Run 03_BatchAdjust using the user-defined variables #######################
Corrected_WD  <- make_stage_dir(Original_WD, paste0("BatchCorrection_", method_input))
add_extension <- paste0("_", method_input)
start_time    <- Sys.time()

cat("\n[INFO] Corrected files will be written to:\n       ", Corrected_WD, "\n", sep = "")

BatchAdjust(basedir         = Original_WD,
            outdir          = Corrected_WD,
            channelsFile    = channels_path,
            batchKeyword    = batchKeyword_input,
            anchorKeyword   = anchorKeyword_input,
            method          = method_input,
            transformation  = transformation_input,
            addExt          = add_extension,
            plotDiagnostics = plotDiagnostics_input)



#### Copy over necessary data files for downstream analysis ####################
# metadata: file names get the batch-correction ending, e.g. ..._95p.fcs
metadata_corrected <- metadata
metadata_corrected$File_Name <- sub("\\.fcs$", paste0(add_extension, ".fcs"), metadata$File_Name)

missing_out <- metadata_corrected$File_Name[!file.exists(file.path(Corrected_WD, metadata_corrected$File_Name))]
if (length(missing_out) > 0) {
  stop(length(missing_out), " corrected file(s) were not created, for example:\n  ",
       paste(head(missing_out, 5), collapse = "\n  "),
       "\nSee the LOG_BatchAdjust file in ", Corrected_WD, call. = FALSE)
}

data.table::fwrite(metadata_corrected, file = file.path(Corrected_WD, METADATA_NAME))
file.copy(file.path(Original_WD, PANEL_NAME),       file.path(Corrected_WD, PANEL_NAME),       overwrite = TRUE)
file.copy(channels_path, file.path(Corrected_WD, basename(channels_path)), overwrite = TRUE)



#### Keep objects in memory for the next scripts ###############################
# Original_WD is unchanged (the FCS folder). The corrected files are in:
Batch_Corrected_Dir <- Corrected_WD
batchCorrection_method <- method_input



#### Record parameters used in this run ########################################
Inputs <- list(
  # --- keys that existed in v5 (names unchanged) ---
  Method          = method_input,
  BatchKeyword    = batchKeyword_input,
  AnchorKeyword   = anchorKeyword_input,
  ChannelsFile    = channelsFile_input,
  Transformation  = transformation_input,
  PlotDiagnostics = plotDiagnostics_input,
  BaseDir         = Original_WD,
  OutputDir       = Corrected_WD,
  AddExtension    = add_extension,
  # --- new in v6 ---
  Script_Version        = SCRIPT_VERSION,
  Batch_Adjust_Script   = batch_script,
  Batch_Adjust_MD5      = unname(tools::md5sum(batch_script_path)),
  Completed_At          = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  Duration_Minutes      = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2),
  Corrected_FCS_Folder  = basename(Corrected_WD),
  Source_FCS_Folder     = Original_WD,
  N_Files_Corrected     = nrow(metadata_corrected),
  Batches_Found         = as.list(batches_found),
  Reference_Batch       = 1,
  Anchor_Files          = as.list(anchor_files),
  N_Channels_Adjusted   = length(channels_to_adjust),
  Channels              = as.list(channels_to_adjust),
  Channels_File_MD5     = as.character(current_md5),
  Random_Seed           = "not used in this script",
  Session               = session_info_list(names(pkg_versions))
)
record_params(MODULE_NAME, Inputs, write_dir = Corrected_WD, parent_dir = Original_WD)



#### Script completed message ##################################################
print_banner(" Batch correction complete ",
             c(paste0("Corrected FCS files, metadata.csv and panel.csv are in:  ", basename(Corrected_WD)),
               "Next: 02_FCS_to_CSV -- extract data from the FCS files and export to CSV files for downstream analysis."))
