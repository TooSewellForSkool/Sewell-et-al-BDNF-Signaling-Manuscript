################################################################################
# 02_FCS_to_CSV.R
################################################################################
# Layer:    DataCleanup -- LAST required script (before Leiden clustering in Python)
# Author:   Jonathon Sewell
# Created:  2026-06-08
# Version:  7   (input version: 02_FCS_to_CSV_v4; v5 was an interim test version)
# Edited:   2026-10-21
#
# Purpose:
# * Read the FCS files (original, or the batch-corrected copies), optionally
#   subsample cells, apply the transforms set in panel.csv, and write the CSV
#   files the Python clustering scripts use:
#     expression_matrix_analysis.csv   channels with Analysis = 1
#     expression_matrix_other.csv      all other channels
#     filenums.csv                     which file (File_Num) each cell came from
#     (subsample_idxs.csv              which cells were kept, if subsampled)
# * Everything goes in a NEW folder next to your FCS files:
#     <date>_<time>_CompiledData_<FullDataSet | Subsampled-N>[...]
#   together with metadata.csv, panel.csv, the clustering scripts and a
#   run_parameters.json that holds the history of every earlier step.
#
# Origin:
# * Subsample and make csv file containing all points asinh transformed
#     > 02_FCS_to_CSV.R - Corey Williams, University of Virginia, May 17, 2019
#     > Updates from 02_FCS_to_CSV.R: formatting and directory selection
#
# ==============================================================================
# CHANGE LOG: v4 -> v6
# ==============================================================================
# Convention:  [RESULT]   = can change the content of the output files
#              [BEHAVIOR] = changes how the script behaves, not the data
#              [SAME]     = output unchanged
#
# v5 -> v6 (after your test run)
#  G1 [BEHAVIOR] No more XQuartz. v5 used svDialogs::dlg_list() for the FCS-folder
#     and local/HPC choices, which needs XQuartz on macOS. All choices are now a
#     normal typed window listing numbered options ("Type the number") -- Helpers v4.
#  G2 [BEHAVIOR] Clustering files by environment, based on your Python/SLURM files:
#       locally : 03_Leiden.py, 04_UMAP.py, 04_UMAP_3D.py
#       HPC     : 03_Leiden_UpdateUMAPImport.py, 04_UMAP_UpdateUMAPImport.py,
#                 04_UMAP_3D_UpdateUMAPImport.py and the *.slurm files that run them
#     v4 copied the local files for both and the SLURM files for HPC, but the SLURM
#     files call the *_UpdateUMAPImport.py files, which v4 never copied. The files
#     are not edited. If a file is not found you are told which one (in your test
#     run none were found: no LeidenClustering folder was inside ScriptsToUse).
#
# v6 -> v7 (LeidenClustering update)
#  G3 [BEHAVIOR] Clustering files now use the single set of Python scripts (local AND
#     HPC) plus Pipeline2_03_Leiden.slurm, Pipeline2_04_UMAP.slurm and
#     Pipeline2_04_UMAP_3D.slurm. The copy step now calls copy_clustering_scripts()
#     in Helpers v6 (shared with 02.1) instead of its own lists. If the new slurm
#     files are not found, the older *_UpdateUMAPImport files are used.
#
# A. Needs your attention
#  A1 [SAME for your data] Cell order / file numbers. v4 read the files in the
#     order of metadata.csv and labelled them 1, 2, 3... by POSITION. v5 reads
#     them in File_Num order and writes the File_Num itself to filenums.csv. These
#     are identical while metadata.csv is in File_Num order starting at 1 (as 01
#     writes it). They differ only if rows were reordered or File_Num values
#     removed (the case 02.1 will handle).
#  A2 [RESULT, safety] Channel names are CHECKED. v4 labelled columns by position
#     (colnames <- panel$Fixed_Param) without checking. v5 stops unless the
#     channel names in every FCS file are exactly panel.csv's Name column, in the
#     same order. (This is the check that protects against a panel whose rows
#     were reordered.)
#  A3 [RESULT, safety] Factor. v4 asked for a transformation factor when any Factor
#     was empty, but then never used the answer (and tested the wrong variable).
#     v5 stops with a message to fix panel.csv (01.2). A Transform other than
#     Asinh or Linear also stops the script (v4 treated it as Linear).
#  A4 [BEHAVIOR] Event counts per file are compared with metadata.csv; the script
#     stops if they differ.
#  A5 [BEHAVIOR] Input choice. If batch-corrected folders exist next to your FCS
#     files you pick which to use (original or a corrected folder) in a list.
#     v4 guessed from folder names containing "Batch". The metadata.csv and
#     panel.csv used are the ones in the chosen folder.
#  A6 [BEHAVIOR] Output folder. v4 created <date>_CompiledData_... inside a /tmp
#     link (the /tmp link does not exist on Windows) with a date only. v5 creates
#     <date>_<time>_CompiledData_..._ as a sibling of your other stage folders,
#     directly in the FCS folder. Name parts: FullDataSet / Subsampled-N,
#     _FilesWithFewerCells, _BatchCorrected-<method>.
#  A7 [BEHAVIOR] Random seed. Subsampling is now reproducible: you may enter a seed
#     (blank = chosen automatically); it is saved in the run record. Cells are
#     drawn the same way as v4 (sample() per file, kept in random order).
#  A8 [BEHAVIOR] "Original_WD" in memory now stays the FCS folder. v4 changed it to
#     the new output folder. The output folder is in 'Compiled_WD'. Tell me if a
#     later script expects the v4 behaviour.
#  A9 [BEHAVIOR] 04_UMAP_3D.py: v4 copied 04_UMAP.py under the name 04_UMAP_3D.py
#     (a typo). v5 copies 04_UMAP_3D.py, and tells you if it is not found.
#     The clustering scripts are located by file name inside ScriptsToUse (any
#     folder with "LeidenClustering" in its name is searched first).
#  A10 [BEHAVIOR] Run record keys from v4 keep their names, with these changes:
#     Original_WD is the FCS folder; BatchCorrected is TRUE/FALSE (v4 stored the
#     position of "Batch" in the folder name); BatchMethod is "none" when not
#     corrected; INPUT_FOLDER / OUTPUT_FOLDER are the folders used.
#
# B. Efficiency / readability / reproducibility
#  B1 [SAME] Files are read one at a time with read.FCS(); v4 loaded all of them
#     into a flowSet first. Matrices are joined once with rbind() (v4 grew the
#     matrix file by file, which gets slow with many cells).
#  B2 [SAME] Division by Factor and asinh vectorised: sweep() instead of a loop
#     over every cell (v4: t(apply(...))). Same values.
#  B3 [SAME] Written with fwrite() as before; drop = FALSE so a single Analysis
#     channel still writes a one-column file (v4 would have failed).
#  B4 Shared helper file; sequence checks (01 and 01.2 recorded as run);
#     progress shown while files are read; a summary of cells per file.
#  B5 Run record gains: seed, input folder, per-file cells kept, channels in each
#     output, scripts copied with checksums, checksums of the CSV files, session.
#  B6 Removed: command-line arguments (INPUT_FOLDER/OUTPUT_FOLDER by args), the
#     unused 'analysis_markers', commented-out code, the local record_params().
#
# C. Suggested, NOT applied
#  * Nothing at present (the HPC file question was settled in G2).
################################################################################


#### USER SETTINGS -- edit this section only ###################################
# (none: everything is chosen in the pop-ups)

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "02_FCS_to_CSV.R"
SCRIPT_VERSION           <- 7
MODULE_NAME              <- "02_FCS_to_CSV"
PREREQUISITE_MODULES     <- c("01_Get_Metadata", "01.2_Update_Metadata")
REQUIRED_HELPERS_VERSION <- 6
RAW_CHOICE               <- "Original FCS files (no batch correction)"



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
cat("[START] Convert preprocessed FCS files for downstream analysis:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

Original_WD <- get_working_dir()   # the folder with your original FCS files



#### Choose the input: original or batch-corrected FCS files ###################
# Batch-corrected folders: stage folders in the FCS folder whose run record has
# a 02_RunBatch section but no 02_FCS_to_CSV section, and that contain FCS files
find_batch_folders <- function(parent) {
  dirs <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[!grepl("^[_.]", basename(dirs))]
  ok <- vapply(dirs, function(d) {
    rec <- read_run_record(d)
    !is.null(rec) && !is.null(rec[["02_RunBatch"]]) && is.null(rec[["02_FCS_to_CSV"]]) &&
      file.exists(file.path(d, METADATA_NAME)) &&
      length(list.files(d, pattern = "\\.fcs$")) > 0
  }, logical(1))
  sort(dirs[ok], decreasing = TRUE)      # newest first (names start with the date)
}

batch_dirs   <- find_batch_folders(Original_WD)
batch_wanted <- identical(get_record_value(read_run_record(Original_WD),
                                           "01_Get_Metadata", "batchCorrection_res"), "yes")
input_dir <- Original_WD
if (length(batch_dirs) > 0) {
  mem_dir <- get0("Batch_Corrected_Dir", envir = .GlobalEnv, inherits = FALSE)
  default_choice <- if (!is.null(mem_dir) && basename(mem_dir) %in% basename(batch_dirs)) {
    basename(mem_dir)
  } else if (batch_wanted) basename(batch_dirs[1]) else RAW_CHOICE
  choice <- ask_choice("Which FCS files should be converted?",
                       c(basename(batch_dirs), RAW_CHOICE), preselect = default_choice)
  if (choice != RAW_CHOICE) {
    input_dir <- batch_dirs[basename(batch_dirs) == choice]
  } else if (batch_wanted &&
             !ask_yes_no("You said in 01 that you would do batch correction, but chose the ORIGINAL files.\n\nContinue with the original files?")) {
    stop("Please run again and choose a batch-corrected folder.", call. = FALSE)
  }
} else if (batch_wanted) {
  if (!ask_yes_no("You said in 01 that you would do batch correction, but no corrected folder was found.\n\nContinue with the original FCS files?")) {
    stop("Run BatchCorrection/02_RunBatch first.", call. = FALSE)
  }
}
input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
is_corrected <- !identical(input_dir, Original_WD)
cat("[INFO] Reading FCS files from:\n       ", input_dir, "\n", sep = "")

# Sequence check and records are read from the chosen folder
run_record <- require_modules(input_dir, PREREQUISITE_MODULES,
                              how_to_fix = "Run 01_Get_Metadata and 01.2_Update_Metadata first.")
if (!isTRUE(get_record_value(run_record, "01.2_Update_Metadata", "Reviewed", default = FALSE))) {
  stop("01.2_Update_Metadata is recorded but not marked as reviewed. ",
       "Run 01.2_Update_Metadata and click Done.", call. = FALSE)
}
batch_method <- if (is_corrected) get_record_value(run_record, "02_RunBatch", "Method", default = "unknown") else "none"



#### Read metadata and panel, and check them ###################################
metadata <- read_metadata(input_dir)
panel    <- read_panel(input_dir)
metadata <- metadata[order(metadata$File_Num), , drop = FALSE]      # File_Num order
n_files  <- nrow(metadata)

if (anyDuplicated(metadata$File_Num)) stop("File_Num is not unique in metadata.csv.", call. = FALSE)
missing_fcs <- metadata$File_Name[!file.exists(file.path(input_dir, metadata$File_Name))]
if (length(missing_fcs) > 0) {
  stop(length(missing_fcs), " FCS file(s) listed in metadata.csv are not in:\n  ", input_dir, "\n  ",
       paste(head(missing_fcs, 5), collapse = "\n  "), call. = FALSE)
}
if (!"Event_Count" %in% colnames(metadata)) {
  stop("metadata.csv has no Event_Count column. Run 01_Get_Metadata again.", call. = FALSE)
}
if (anyNA(panel$Factor)) {
  stop("panel.csv has an empty or non-numeric Factor for: ",
       paste(panel$Name[is.na(panel$Factor)], collapse = ", "),
       "\nSet it in 01.2_Update_Metadata (5 for Asinh channels, 1 for Linear).", call. = FALSE)
}
if (any(panel$Factor == 0)) {
  stop("Factor is 0 for: ", paste(panel$Name[panel$Factor == 0], collapse = ", "), call. = FALSE)
}
if (!all(panel$Transform %in% c("Asinh", "Linear"))) {
  stop("Transform must be Asinh or Linear. Check: ",
       paste(panel$Name[!(panel$Transform %in% c("Asinh", "Linear"))], collapse = ", "), call. = FALSE)
}
if (anyDuplicated(panel$Fixed_Param)) stop("Fixed_Param is not unique in panel.csv.", call. = FALSE)
if (!any(panel$Analysis == 1)) stop("No channel has Analysis = 1 in panel.csv.", call. = FALSE)



#### Define input information (pop-ups) ########################################
SUBSAMPLE_BOOL <- !ask_yes_no("Do you want to use the full dataset? \nYes = use entire dataset \nNo = subsample the dataset")
NUM_SUBSAMPLES <- NA
update_res     <- "Empty"
seed_info      <- NULL

if (SUBSAMPLE_BOOL) {
  ask_number <- function(message, default) {
    res <- svDialogs::dlg_input(message = message, default = default)$res
    if (length(res) != 1) stop("The pop-up window was closed without an answer.", call. = FALSE)
    suppressWarnings(as.numeric(res))
  }
  NUM_SUBSAMPLES <- ask_number("Enter number of cells per file (subsampling):", "5000")
  if (is.na(NUM_SUBSAMPLES) || NUM_SUBSAMPLES < 1) stop("Invalid number entered for NUM_SUBSAMPLES")
  NUM_SUBSAMPLES <- floor(NUM_SUBSAMPLES)

  # Files with fewer cells than requested: keep going or change the number
  repeat {
    n_cells  <- metadata$Event_Count
    too_few  <- which(n_cells < NUM_SUBSAMPLES)
    if (length(too_few) == 0) break
    msg <- paste0(length(too_few), " file(s) contain fewer than NUM_SUBSAMPLES (", NUM_SUBSAMPLES,
                  ") cells:\n",
                  paste0(metadata$File_Name[too_few], " (", n_cells[too_few], " cells)", collapse = "\n"),
                  "\n\nWould you like to update NUM_SUBSAMPLES?")
    if (ask_yes_no(msg)) {
      new_n <- ask_number("Enter new value for NUM_SUBSAMPLES:", as.character(NUM_SUBSAMPLES))
      if (is.na(new_n) || new_n < 1) {
        svDialogs::dlg_message("Invalid entry - keeping previous NUM_SUBSAMPLES value.", type = "ok")
      } else {
        NUM_SUBSAMPLES <- floor(new_n)
      }
    } else {
      update_res <- "no"      # keep all cells of the small files
      break
    }
  }

  seed_text <- svDialogs::dlg_input(
    message = "Random seed for subsampling.\n\nLeave blank to choose one automatically (it is saved in the run record).\nEnter a number to repeat an earlier subsample exactly.",
    default = "")$res
  if (length(seed_text) != 1) stop("The pop-up window was closed without an answer.", call. = FALSE)
  seed_value <- if (nzchar(trimws(seed_text))) suppressWarnings(as.numeric(seed_text)) else NULL
  if (!is.null(seed_value) && is.na(seed_value)) stop("The seed must be a number.", call. = FALSE)
  seed_info <- setup_seed(seed_value)
}

cluster_res <- ask_choice("Will you be clustering locally or on an HPC?",
                          c("locally", "HPC"), preselect = "locally")



#### Read FCS files, subsample, join ###########################################
cat("\n[INFO] Reading ", n_files, " FCS files", if (SUBSAMPLE_BOOL) " and subsampling" else "",
    ": this may take a moment\n", sep = "")
mats        <- vector("list", n_files)
sample_idx  <- vector("list", n_files)
for (i in seq_len(n_files)) {
  fr <- flowCore::read.FCS(file.path(input_dir, metadata$File_Name[i]),
                           transformation = FALSE, truncate_max_range = FALSE)
  m  <- flowCore::exprs(fr)

  # Checks: channel names/order and number of events
  if (!identical(unname(colnames(m)), panel$Name)) {
    stop("The channels in ", metadata$File_Name[i], " do not match panel.csv (Name column), ",
         "or are in a different order.\n  FCS file : ", paste(head(colnames(m), 6), collapse = ", "), " ...",
         "\n  panel.csv: ", paste(head(panel$Name, 6), collapse = ", "), " ...",
         "\nColumns are named by position, so the order must be identical.", call. = FALSE)
  }
  if (nrow(m) != metadata$Event_Count[i]) {
    stop(metadata$File_Name[i], " has ", nrow(m), " events but metadata.csv says ",
         metadata$Event_Count[i], ". Run 01_Get_Metadata again.", call. = FALSE)
  }

  if (SUBSAMPLE_BOOL) {
    idx <- sample.int(nrow(m), min(NUM_SUBSAMPLES, nrow(m)))
    sample_idx[[i]] <- idx
    m <- m[idx, , drop = FALSE]
  }
  mats[[i]] <- m
  if (i %% 10 == 0 || i == n_files) cat("       ", i, "of", n_files, "files read\n")
}
cells_per_file <- vapply(mats, nrow, integer(1))
data_mat <- do.call(rbind, mats)
rm(mats, fr, m)
filenums <- rep(metadata$File_Num, times = cells_per_file)



#### Transform #################################################################
colnames(data_mat) <- panel$Fixed_Param
# Divide each channel by its Factor (1 for linear channels)
data_mat <- sweep(data_mat, 2, panel$Factor, "/")
# asinh transform for channels marked Asinh
asinh_channels <- which(panel$Transform == "Asinh")
data_mat[, asinh_channels] <- asinh(data_mat[, asinh_channels])



#### Save concatenated exprs arrays and filenums ###############################
Compiled_WD <- make_stage_dir(
  Original_WD,
  paste0("CompiledData",
         if (SUBSAMPLE_BOOL) paste0("_Subsampled-", NUM_SUBSAMPLES) else "_FullDataSet",
         if (update_res == "no") "_FilesWithFewerCells",
         if (is_corrected) paste0("_BatchCorrected-", batch_method)))
cat("\n[INFO] Writing files to:\n       ", Compiled_WD, "\n", sep = "")

is_analysis <- panel$Analysis == 1
data.table::fwrite(data_mat[, is_analysis,  drop = FALSE], file = file.path(Compiled_WD, "expression_matrix_analysis.csv"))
data.table::fwrite(data_mat[, !is_analysis, drop = FALSE], file = file.path(Compiled_WD, "expression_matrix_other.csv"))
data.table::fwrite(matrix(filenums, ncol = 1),             file = file.path(Compiled_WD, "filenums.csv"))

if (SUBSAMPLE_BOOL) {
  # Same layout as v4: one column per file (File_Num order), padded with blanks
  max_len <- max(lengths(sample_idx))
  subsample_idxs <- do.call(cbind, lapply(sample_idx, function(idx) { length(idx) <- max_len; idx }))
  colnames(subsample_idxs) <- metadata$File_Name
  data.table::fwrite(subsample_idxs, file = file.path(Compiled_WD, "subsample_idxs.csv"),
                     col.names = FALSE, sep = ",")
}

## Copy over metadata and panel (from the folder the FCS files came from)
file.copy(file.path(input_dir, METADATA_NAME), file.path(Compiled_WD, METADATA_NAME), overwrite = TRUE)
file.copy(file.path(input_dir, PANEL_NAME),    file.path(Compiled_WD, PANEL_NAME),    overwrite = TRUE)



#### Copy over scripts for downstream analysis #################################
# File lists and the HPC/legacy fallback live in Helpers (copy_clustering_scripts).
script_res <- copy_clustering_scripts(Compiled_WD, cluster_res)
py_res        <- script_res
not_found_all <- script_res$not_found
copied_files  <- script_res$copied



#### Keep objects in memory for the next scripts ###############################
# Original_WD is unchanged (the FCS folder). The new folder is:
metadata_compiled <- metadata
panel_compiled    <- panel
# Also in memory: Compiled_WD



#### Record parameters used in this run ########################################
csv_files <- file.path(Compiled_WD, c("expression_matrix_analysis.csv", "expression_matrix_other.csv",
                                       "filenums.csv", if (SUBSAMPLE_BOOL) "subsample_idxs.csv"))
Inputs <- list(
  # --- keys that existed in v4 (names unchanged; see change log A10) ---
  Original_WD             = Original_WD,
  SUBSAMPLE_BOOL          = SUBSAMPLE_BOOL,
  NUM_SUBSAMPLES          = if (SUBSAMPLE_BOOL) NUM_SUBSAMPLES else NA,
  UsedFilesWithFewerCells = update_res == "no",
  CopyClusteringScripts   = length(py_res$copied) > 0,
  ComputeEnvironment      = cluster_res,
  INPUT_FOLDER            = input_dir,
  OUTPUT_FOLDER           = Compiled_WD,
  BatchCorrected          = is_corrected,
  BatchMethod             = batch_method,
  # --- new in v5 ---
  Script_Version          = SCRIPT_VERSION,
  Completed_At            = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  Input_Folder_Name       = basename(input_dir),
  N_Files                 = n_files,
  N_Cells_Total           = nrow(data_mat),
  Cells_Per_File          = lapply(seq_len(n_files), function(i)
    list(File_Num = metadata$File_Num[i], File_Name = metadata$File_Name[i],
         Cells_In_File = metadata$Event_Count[i], Cells_Kept = cells_per_file[i])),
  Analysis_Channels       = as.list(panel$Fixed_Param[is_analysis]),
  Other_Channels          = as.list(panel$Fixed_Param[!is_analysis]),
  Asinh_Channels          = as.list(panel$Fixed_Param[asinh_channels]),
  Random_Seed             = if (SUBSAMPLE_BOOL) seed_info else "not used (full dataset)",
  Scripts_Copied          = as.list(copied_files),
  Scripts_Not_Found       = as.list(not_found_all),
  Output_MD5              = file_md5(csv_files),
  Session                 = session_info_list(names(pkg_versions))
)
record_params(MODULE_NAME, Inputs, write_dir = Compiled_WD, parent_dir = input_dir)



#### Script completed message ##################################################
print_banner(" Expression files output complete ",
             c(paste0("Files are in:  ", basename(Compiled_WD)),
               "Time to get multilingual on this analysis: switch over to Python for 03_Leiden.py and 04_UMAP.py."))
