################################################################################
# 01_PercAboveThreshold.R
################################################################################
# Layer:    Analysis -- data processing (no plots; plots are in the Visualization layer)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-23
# Edited:   2026-10-27
# Version:  4   (1 = NEW: replaces the table-making half of 10.1, 12.1 and 12.1.1;
#                2 = optional expression categories (12.1.1 'ThreshCat' mode), optional Inf -> NA;
#                3 = optional rows for sample/cluster combinations with no cells (12.1.1 Clusters mode),
#                4 = optional MEASURE = "mean_expression" (mean of a marker per sample and cluster, as Fig 5F))
#
# Purpose:
# * For every sample (FCS file) and cluster, calculate the PERCENT OF CELLS ABOVE A
#   THRESHOLD for every marker, then the FOLD CHANGE (and log2 fold change) against
#   a control condition. This is the table behind the signaling heatmaps/dot plots.
#     - threshold = the 95th percentile (THRESHOLD_PROB) of the baseline samples
#       (the samples at THRESHOLD_TIMEPOINT, e.g. NoTreat 0 hr): the percentile is
#       taken per baseline file, and the files are averaged
#     - ANALYSIS_BY = "PseudoBulk": one threshold from all cells (Cluster = 0)
#       ANALYSIS_BY = "Clusters":   one threshold per cluster
# * One script replaces the table-making code that was repeated in
#   10.1_..._PseudoBulk_..., 12.1_..._ByCluster_... and 12.1.1_..._PseudoBulk/Clusters.
#
# Origin (the calculations are those of these published scripts; see README):
# * 10.1_PlotSignalingHeatmap_PseudoBulk_PercAboveThreshold_v6 (f34b3f85)
# * 12.1_PlotSignalingHeatmap_ByCluster_BySample_PercAboveThreshold_NormToCtrl_v6 (c713ed46)
# * 12.1.1_PlotDotPlots_ByCluster_BySample_PercAboveThreshold_NormToCtrl_v9 (PseudoBulk mode;
#   v2 also its ThresholdCategory mode, see THRESHOLD_CATEGORY_FOLDER)
#
# ==============================================================================
# CHANGE LOG: published scripts -> v1 (v2 additions: A9-A11, B4; v3 addition: A12; v4 addition: A13)
# ==============================================================================
# Convention:  [RESULT]   = can change the content of the output tables
#              [BEHAVIOR] = changes how the script behaves, not the numbers
#
# A. Needs your attention (nothing here is applied silently; each item has a check)
#  A1 [RESULT: none expected] Same calculation as the published scripts. The golden
#     check (Tools/Compare_GoldenFiles.R, settings in README_Analysis.txt) must show
#     identical numbers for the published RAWvalues tables before this is used for
#     anything new. Numeric differences up to 1e-9 are rounding only.
#  A2 [RESULT: none for your data] Samples and results are paired by File_Num. 10.1
#     pasted the metadata next to the results BY ROW POSITION, which is only correct
#     when metadata.csv is sorted by File_Num (yours is).
#  A3 [BEHAVIOR] Samples excluded with FILES_TO_EXCLUDE still contribute to the
#     threshold if they are baseline samples (as in the published scripts), and are
#     removed from the tables afterwards. This is recorded in the run record.
#  A4 [RESULT: none] The marker-class rule that decides which markers go into the
#     fold-change tables (published: names starting with "p", plus H3K9ac, TrkB_s,
#     yH2AX, CC3, TrkB_i, minus DNA) is kept. Set MARKERS_FOR_DERIVED to choose them
#     yourself. A "Class" column in panel.csv would be more robust (suggestion only).
#  A5 [RESULT, declared] log2 fold change is -Inf when a sample has 0 % cells above
#     threshold and NaN when the control mean is 0 (as published). RepsAvg tables
#     therefore also show -Inf. Not changed; listed so it is not a surprise.
#  A6 [RESULT, declared] The "RepsAvg" tables no longer average the Replicate
#     number (published: a meaningless mean of 1,2,3); they have N_Replicates.
#  A7 [BEHAVIOR, declared] Column names are one set for every table:
#     File_Num, File_Name, Sample_Condition, Reordered, TimePoint, Replicate,
#     Cluster, N_Cells, then markers. (Published: Condition / Sample_Condition and
#     TimePoint / TimePoints, depending on the script; Cluster NA in 10.1, 0 in
#     12.1.1.) Extra columns: N_Cells (cells behind each percentage), File_Num in
#     every table. Pseudobulk Cluster = 0.
#  A8 [BEHAVIOR] New output Thresholds.csv: the threshold value used for every
#     marker (and cluster), with the number of baseline cells and files.
#
#  A9 [BEHAVIOR, new in v2] THRESHOLD_CATEGORY_FOLDER (default "" = off). When it points at a
#     folder made by 02_AssignThresholdCategories.R, every cell carries its expression
#     category and the whole calculation is done WITHIN each category (and cluster), as
#     12.1.1 does: the threshold comes from the baseline cells of that category only,
#     percentages are of that category's cells in each sample, and the fold-change control
#     is the baseline mean of the SAME category. Tables get a ThreshCat column. Samples with
#     no cells in a category have no row (as published); they are listed in the run record.
#     With the setting off, v2 gives the same tables as v1 (checked with the golden runs).
#     NOTE: 12.1 (heatmap script) in category mode took the threshold from category cells
#     but the percentages from ALL cells; that mixed mode is not reproduced (never used for
#     a published figure that I could find).
#  A10 [RESULT, optional] LOG2FC_INF_AS_NA (default FALSE = as 10.1/12.1 wrote their tables:
#     -Inf stays in the table, and the plot code turns it into NA later). 12.1.1 turned -Inf
#     into NA BEFORE averaging replicates, so its RepsAvg tables average only the replicates
#     that had a value. The two can differ for any condition where a replicate has 0 % cells
#     above threshold. TRUE reproduces 12.1.1. None of the tables checked so far contained -Inf.
#  A11 [BEHAVIOR] The run record now also lists sample/cluster/category combinations with no cells.
#
#  A12 [RESULT, optional, new in v3] FILL_EMPTY_CLUSTERS (default FALSE = as v2). In ANALYSIS_BY = "Clusters" mode a
#     sample that has NO cells in a cluster has no row. 12.1.1 added a row for every sample x cluster
#     (complete(Cluster, File_Num)) with NA for every marker; with 20 clusters (two very small ones)
#     that matters, the published 20-cluster tables have 45 x 20 = 900 rows. TRUE reproduces that: the
#     missing rows are added with N_Cells = 0 and NA values (their fold changes are NA, and a replicate-average
#     over only NA values is NaN, as published). Not available together with expression categories.
#     With the setting FALSE v3 gives the same tables as v2 (only this option was added).
#
#  A13 [RESULT, optional, new in v4] MEASURE (default "percent_above" = exactly as v3). With MEASURE = "mean_expression"
#     the table value is the MEAN of the marker over the cells of each sample (and cluster) instead of the percent of
#     cells above a threshold; fold change and log2 fold change are then taken against the baseline mean exactly as
#     before, and the tables have the same layout, so 03/04/05 and the plotting scripts work on them unchanged.
#     This is the quantity of the published Fig 5F (SelectPlots_v9: mean per sample and cluster, divided by the mean
#     of the 0 hr samples of that cluster). SelectPlots averaged the values AS STORED in expression_matrix_analysis.csv
#     (the arcsinh scale, no back-transformation); MEAN_SCALE "as_stored" (default) does the same. "back_transformed"
#     averages 5 * sinh(x) instead (a different quantity, not published). No threshold is calculated in this mode
#     (no Thresholds.csv); the baseline samples are still those at THRESHOLD_TIMEPOINT. The tables are named
#     Expt<set>_MeanExpr_<...>.csv. The golden check is in Inventory (Fig5F_points_extracted_from_PDFs.csv).
#
# B. Efficiency / readability / reproducibility
#  B1 Cells are read once; metadata (Expt/Expt_ID, label columns, optional Reordered)
#     is mapped automatically for both the published and the new metadata layout.
#  B2 Run record: all settings, files and cells kept at each filter, thresholds,
#     marker lists, warnings, checksums, session.
#  B3 No random numbers are used, so no seed is needed.
#  B4 (v2) Category files are read as text and checked against the cell count; nothing
#     is copied into the output folder.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
INPUT_FOLDER        <- ""         # folder with the expression matrices ("" = choose in a pop-up)
CLUSTER_ROUND       <- 1          # which cluster_R<N>_assigns.csv to use
CLUSTERS_MERGED     <- FALSE      # TRUE = use cluster_R<N>_assigns_merged.csv

ANALYSIS_BY         <- "Clusters"   # "PseudoBulk" (all cells together) or "Clusters"

# Which samples are analysed
SAMPLE_SET_COLUMN   <- "Expt"     # metadata column that defines the sample set ("Expt" or "Expt_ID")
SAMPLE_SET          <- c(1,2,3,4)
CONDITIONS_TO_EXCLUDE <- c()      # base treatments to leave out, e.g. c("CompMed") (c() = none)
CLUSTERS_TO_REMOVE    <- c()      # clusters to leave out of everything (c() = none)
FILES_TO_EXCLUDE      <- c()      # File_Name values removed from the tables (see change log A3)

# What is measured (see change log A13)
MEASURE             <- "percent_above"   # "percent_above" (default) or "mean_expression" (mean per sample and cluster, Fig 5F)
MEAN_SCALE          <- "as_stored"       # only for "mean_expression": "as_stored" (arcsinh scale, as Fig 5F) or "back_transformed"

# How "responsive" is defined
THRESHOLD_PROB      <- 0.95       # 0.95 = 95th percentile of the baseline samples
THRESHOLD_TIMEPOINT <- 0          # baseline = samples with this TimePoint (hours)

# What fold change is calculated against
CTRL_PER_TIMEPOINT  <- FALSE      # FALSE: every sample vs the baseline condition(s)
                                  # TRUE:  every sample vs the control at its OWN timepoint
NORM_CONDITIONS     <- NULL       # NULL = the baseline condition(s) (CTRL_PER_TIMEPOINT FALSE);
                                  # REQUIRED when CTRL_PER_TIMEPOINT is TRUE,
                                  # e.g. c("NoTreat_0hr", "ContDepr_5min", ..., "ContDepr_16hr")

# Expression categories (optional): analyse WITHIN categories made by 02_AssignThresholdCategories.R
THRESHOLD_CATEGORY_FOLDER <- ""   # "" = no categories; or the folder that 02_AssignThresholdCategories.R made
CATEGORIES_TO_ANALYSE     <- c()  # e.g. c("1", "2") or c("11", "41", "14", "44"); c() = every category

# Advanced (leave as they are unless you know you need to change them)
LOG2FC_INF_AS_NA    <- TRUE      # FALSE = as 10.1/12.1 tables; TRUE = as 12.1.1 tables (see change log A10)
FILL_EMPTY_CLUSTERS <- TRUE      # TRUE = add NA rows for sample x cluster pairs with no cells, as 12.1.1 (see change log A12)
LABEL_COLUMN        <- NULL       # metadata column holding labels like "BDNF_15min" (NULL = detected)
MARKERS_FOR_DERIVED <- NULL       # NULL = published marker-class rule (see A4)
OUTPUT_LABEL        <- "CL20"         # optional text added to the output folder name

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "01_PercAboveThreshold.R"
SCRIPT_VERSION           <- 4
MODULE_BASE_NAME         <- "01_PercAboveThreshold"
REQUIRED_HELPERS_VERSION <- 7
EXPR_ANALYSIS_NAME       <- "expression_matrix_analysis.csv"
EXPR_OTHER_NAME          <- "expression_matrix_other.csv"
FILENUMS_NAME            <- "filenums.csv"
CATEGORY_FILE_NAME       <- "ThresholdCategory.csv"
ASINH_FACTOR             <- 5      # the published scripts back-transform every channel with 5 * sinh(x)



#### Load shared helper functions (do not edit this block) #####################
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
pkg_versions <- dc_load_packages(cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite"))

cat("\n========================================================\n")
cat("[START] Percent above threshold and fold change:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

warnings_log <- character(0)
note_warning <- function(txt) {
  warnings_log <<- c(warnings_log, txt)
  cat("[WARNING] ", txt, "\n", sep = "")
}

# Check the settings before doing any work
if (!ANALYSIS_BY %in% c("PseudoBulk", "Clusters")) {
  stop("ANALYSIS_BY must be \"PseudoBulk\" or \"Clusters\".", call. = FALSE)
}
if (!is.numeric(THRESHOLD_PROB) || length(THRESHOLD_PROB) != 1 ||
    THRESHOLD_PROB <= 0 || THRESHOLD_PROB >= 1) {
  stop("THRESHOLD_PROB must be a single number between 0 and 1 (e.g. 0.95).", call. = FALSE)
}
if (isTRUE(CTRL_PER_TIMEPOINT) && length(NORM_CONDITIONS) == 0) {
  stop("CTRL_PER_TIMEPOINT is TRUE, so NORM_CONDITIONS must list the control condition(s),\n",
       "e.g. NORM_CONDITIONS <- c(\"NoTreat_0hr\", \"ContDepr_5min\", \"ContDepr_15min\").", call. = FALSE)
}
by_cluster <- ANALYSIS_BY == "Clusters"
use_cat    <- nzchar(THRESHOLD_CATEGORY_FOLDER)
if (!is.logical(LOG2FC_INF_AS_NA) || length(LOG2FC_INF_AS_NA) != 1 || is.na(LOG2FC_INF_AS_NA)) {
  stop("LOG2FC_INF_AS_NA must be TRUE or FALSE.", call. = FALSE)
}
if (!is.logical(FILL_EMPTY_CLUSTERS) || length(FILL_EMPTY_CLUSTERS) != 1 || is.na(FILL_EMPTY_CLUSTERS)) {
  stop("FILL_EMPTY_CLUSTERS must be TRUE or FALSE.", call. = FALSE)
}
if (FILL_EMPTY_CLUSTERS && ANALYSIS_BY != "Clusters") {
  stop("FILL_EMPTY_CLUSTERS only applies to ANALYSIS_BY = \"Clusters\".", call. = FALSE)
}
if (FILL_EMPTY_CLUSTERS && nzchar(THRESHOLD_CATEGORY_FOLDER)) {
  stop("FILL_EMPTY_CLUSTERS cannot be combined with expression categories (THRESHOLD_CATEGORY_FOLDER).", call. = FALSE)
}
if (!use_cat && length(CATEGORIES_TO_ANALYSE) > 0) {
  stop("CATEGORIES_TO_ANALYSE is set but THRESHOLD_CATEGORY_FOLDER is empty.", call. = FALSE)
}



if (!MEASURE %in% c("percent_above", "mean_expression")) stop("MEASURE must be \"percent_above\" or \"mean_expression\".", call. = FALSE)
if (!MEAN_SCALE %in% c("as_stored", "back_transformed")) stop("MEAN_SCALE must be \"as_stored\" or \"back_transformed\".", call. = FALSE)
use_mean <- identical(MEASURE, "mean_expression")


#### Choose the input folder ###################################################
if (nzchar(INPUT_FOLDER)) {
  if (!dir.exists(INPUT_FOLDER)) stop("INPUT_FOLDER not found:\n  ", INPUT_FOLDER, call. = FALSE)
  input_dir <- normalizePath(INPUT_FOLDER, winslash = "/", mustWork = TRUE)
} else {
  Original_WD <- get_working_dir("project")
  here_has_data <- file.exists(file.path(Original_WD, EXPR_ANALYSIS_NAME))
  dirs <- list.dirs(Original_WD, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[!grepl("^[_.]", basename(dirs))]
  dirs <- dirs[file.exists(file.path(dirs, EXPR_ANALYSIS_NAME))]
  dirs <- sort(dirs, decreasing = TRUE)                       # newest first
  if (here_has_data && length(dirs) == 0) {
    input_dir <- Original_WD
  } else if (length(dirs) == 0) {
    stop("No folder with ", EXPR_ANALYSIS_NAME, " was found in:\n  ", Original_WD,
         "\nChoose the right folder, or set INPUT_FOLDER at the top of this script.", call. = FALSE)
  } else {
    choices <- c(if (here_has_data) basename(Original_WD), basename(dirs))
    pick <- ask_choice("Which data folder should be analysed?", choices, preselect = choices[1])
    input_dir <- if (here_has_data && pick == basename(Original_WD)) Original_WD else dirs[basename(dirs) == pick]
    input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  }
}
cat("[INFO] Input folder:\n       ", input_dir, "\n", sep = "")

# Sequence check: new data must come from the DataCleanup scripts; published
# (older) folders have no run record and are accepted with a note.
input_record <- read_run_record(input_dir)
if (is.null(input_record)) {
  note_warning(paste0("No ", RUN_RECORD_NAME, " in the input folder: treated as data made before the run ",
                      "records existed (e.g. the published dataset). Provenance cannot be checked."))
  input_provenance <- "no run record (legacy data)"
} else {
  require_modules(input_dir, c("01_Get_Metadata", "02_FCS_to_CSV"),
                  how_to_fix = "Run the DataCleanup scripts first.")
  input_provenance <- "run record present; DataCleanup modules found"
}



#### Read the small files ######################################################
metadata <- read_metadata(input_dir)
panel    <- read_panel(input_dir)
if (anyDuplicated(panel$Fixed_Param)) stop("Fixed_Param is not unique in panel.csv.", call. = FALSE)

filenums <- data.table::fread(file.path(input_dir, FILENUMS_NAME), data.table = FALSE)[[1]]
n_cells_all <- length(filenums)
if (!all(unique(filenums) %in% metadata$File_Num)) {
  stop("filenums.csv contains File_Num values that are not in metadata.csv.", call. = FALSE)
}

cluster_file <- sprintf("cluster_R%d_assigns%s.csv", CLUSTER_ROUND, if (CLUSTERS_MERGED) "_merged" else "")
cluster_path <- file.path(input_dir, cluster_file)
if (!file.exists(cluster_path)) {
  stop(cluster_file, " was not found in the input folder.\n",
       if (CLUSTERS_MERGED) "Create it with the cluster-merge script, or set CLUSTERS_MERGED <- FALSE.\n" else "",
       "Folder: ", input_dir, call. = FALSE)
}
clusters <- data.table::fread(cluster_path, data.table = FALSE)[[1]]
if (length(clusters) != n_cells_all) {
  stop(cluster_file, " has ", length(clusters), " rows but ", FILENUMS_NAME, " has ", n_cells_all,
       ". The cluster file belongs to a different dataset.", call. = FALSE)
}



#### Expression categories (optional) ##########################################
cat_all <- NULL
category_info <- list(Used = FALSE)
if (use_cat) {
  if (!dir.exists(THRESHOLD_CATEGORY_FOLDER)) {
    stop("THRESHOLD_CATEGORY_FOLDER not found:\n  ", THRESHOLD_CATEGORY_FOLDER, call. = FALSE)
  }
  cat_dir  <- normalizePath(THRESHOLD_CATEGORY_FOLDER, winslash = "/", mustWork = TRUE)
  cat_path <- file.path(cat_dir, CATEGORY_FILE_NAME)
  if (!file.exists(cat_path)) {
    stop(CATEGORY_FILE_NAME, " was not found in THRESHOLD_CATEGORY_FOLDER.\nFolder: ", cat_dir, call. = FALSE)
  }
  cat_all <- data.table::fread(cat_path, colClasses = "character", na.strings = c("", "NA"), data.table = FALSE)[[1]]
  if (length(cat_all) != n_cells_all) {
    stop(CATEGORY_FILE_NAME, " has ", length(cat_all), " rows but ", FILENUMS_NAME, " has ", n_cells_all,
         ". The categories belong to a different dataset.", call. = FALSE)
  }
  cat_record <- read_run_record(cat_dir)
  if (is.null(cat_record) || is.null(cat_record[["02_AssignThresholdCategories"]])) {
    note_warning(paste0("THRESHOLD_CATEGORY_FOLDER has no run record from 02_AssignThresholdCategories.R: ",
                        "category provenance cannot be checked (e.g. a 13.1 output folder)."))
    cat_provenance <- "no 02_AssignThresholdCategories record"
  } else {
    cat_provenance <- "run record from 02_AssignThresholdCategories present"
    rec_input <- cat_record[["02_AssignThresholdCategories"]][["INPUT_FOLDER"]]
    if (!is.null(rec_input) && !identical(normalizePath(rec_input, winslash = "/", mustWork = FALSE), input_dir)) {
      note_warning(paste0("The categories were made from a different input folder than this analysis:\n  ",
                          rec_input, "\nCell counts match, but check that this is intended."))
    }
  }
  if (anyNA(cat_all)) {
    note_warning(paste0(sum(is.na(cat_all)), " cell(s) have no category (missing expression value); they are left out."))
  }
  cats_present <- sort(unique(cat_all[!is.na(cat_all)]))
  cats_wanted  <- if (length(CATEGORIES_TO_ANALYSE) == 0) cats_present else as.character(CATEGORIES_TO_ANALYSE)
  cats_missing <- setdiff(cats_wanted, cats_present)
  if (length(cats_missing) > 0) {
    stop("CATEGORIES_TO_ANALYSE contains categories that no cell has: ", paste(cats_missing, collapse = ", "),
         "\nCategories in the file: ", paste(cats_present, collapse = ", "), call. = FALSE)
  }
  cat("[INFO] Categories analysed: ", paste(cats_wanted, collapse = ", "), "\n", sep = "")
  category_info <- list(Used = TRUE, Folder = cat_dir, File_MD5 = unname(tools::md5sum(cat_path)),
                        Provenance = cat_provenance, Categories_In_File = as.list(cats_present),
                        Categories_Analysed = as.list(cats_wanted))
}



#### Map the metadata columns (published layout or new layout) #################
# Published metadata: Condition = "BDNF_15min" (label), Treatment = "BDNF" (base)
# New metadata (01_Get_Metadata v8): Treatment = "BDNF_15min" (label), Condition = "BDNF"
looks_like_label <- function(x) mean(grepl("_[0-9.]+(hr|min)$", as.character(x)), na.rm = TRUE)
if (is.null(LABEL_COLUMN)) {
  cand <- intersect(c("Treatment", "Condition"), colnames(metadata))
  if (length(cand) == 0) stop("metadata.csv needs a Treatment or Condition column.", call. = FALSE)
  LABEL_COLUMN_USED <- cand[which.max(vapply(cand, function(cn) looks_like_label(metadata[[cn]]), numeric(1)))]
} else {
  if (!LABEL_COLUMN %in% colnames(metadata)) stop("LABEL_COLUMN '", LABEL_COLUMN, "' is not in metadata.csv.", call. = FALSE)
  LABEL_COLUMN_USED <- LABEL_COLUMN
}
base_cand <- setdiff(intersect(c("Treatment", "Condition"), colnames(metadata)), LABEL_COLUMN_USED)
BASE_COLUMN_USED <- if (length(base_cand) > 0) base_cand[1] else NA_character_

if (!SAMPLE_SET_COLUMN %in% colnames(metadata)) {
  alt <- intersect(c("Expt", "Expt_ID"), colnames(metadata))
  if (length(alt) == 0) {
    stop("metadata.csv has no column '", SAMPLE_SET_COLUMN, "' (or Expt / Expt_ID). Set SAMPLE_SET_COLUMN.", call. = FALSE)
  }
  note_warning(paste0("SAMPLE_SET_COLUMN '", SAMPLE_SET_COLUMN, "' not found; using '", alt[1], "'."))
  SAMPLE_SET_COLUMN <- alt[1]
}
if (!"TimePoint" %in% colnames(metadata)) stop("metadata.csv has no TimePoint column.", call. = FALSE)
if (!"Replicate" %in% colnames(metadata)) stop("metadata.csv has no Replicate column.", call. = FALSE)

md <- data.frame(
  File_Num         = metadata$File_Num,
  File_Name        = metadata$File_Name,
  Sample_Condition = as.character(metadata[[LABEL_COLUMN_USED]]),
  Base_Treatment   = if (is.na(BASE_COLUMN_USED)) sub("_.*", "", as.character(metadata[[LABEL_COLUMN_USED]]))
                     else as.character(metadata[[BASE_COLUMN_USED]]),
  Sample_Set_Value = metadata[[SAMPLE_SET_COLUMN]],
  TimePoint        = suppressWarnings(as.numeric(metadata$TimePoint)),
  Replicate        = suppressWarnings(as.numeric(metadata$Replicate)),
  stringsAsFactors = FALSE)
if ("Reordered" %in% colnames(metadata)) {
  md$Reordered <- suppressWarnings(as.numeric(metadata$Reordered))
  reordered_source <- "metadata.csv column Reordered"
} else {
  ord <- order(!(md$TimePoint %in% THRESHOLD_TIMEPOINT), md$Base_Treatment, md$TimePoint, md$Replicate)
  md$Reordered <- NA_real_
  md$Reordered[ord] <- seq_along(ord)
  reordered_source <- "derived (baseline first, then treatment, timepoint, replicate): metadata.csv has no Reordered column"
  note_warning("metadata.csv has no 'Reordered' column; the plotting order was derived automatically.")
}
cat("[INFO] Sample label column: ", LABEL_COLUMN_USED, "; base treatment column: ",
    ifelse(is.na(BASE_COLUMN_USED), "(derived from the label)", BASE_COLUMN_USED), "\n", sep = "")



#### Choose samples ############################################################
in_set <- md$Sample_Set_Value %in% SAMPLE_SET
if (!any(in_set)) stop("No samples have ", SAMPLE_SET_COLUMN, " in: ", paste(SAMPLE_SET, collapse = ", "), call. = FALSE)
unknown_excl <- setdiff(CONDITIONS_TO_EXCLUDE, md$Base_Treatment)
if (length(unknown_excl) > 0) {
  note_warning(paste0("CONDITIONS_TO_EXCLUDE not found in metadata (no effect): ", paste(unknown_excl, collapse = ", ")))
}
unknown_files <- setdiff(FILES_TO_EXCLUDE, md$File_Name)
if (length(unknown_files) > 0) {
  note_warning(paste0("FILES_TO_EXCLUDE not found in metadata (no effect): ", paste(unknown_files, collapse = ", ")))
}
keep_files <- md$File_Num[in_set & !(md$Base_Treatment %in% CONDITIONS_TO_EXCLUDE)]
md_used <- md[md$File_Num %in% keep_files, , drop = FALSE]
cat("[INFO] Samples in the analysis: ", nrow(md_used), " (", SAMPLE_SET_COLUMN, " in ",
    paste(SAMPLE_SET, collapse = ","), ")\n", sep = "")



#### Read the expression matrices ##############################################
plot_markers <- panel$Fixed_Param[panel$Plot == 1]
if (length(plot_markers) == 0) stop("No channels have Plot = 1 in panel.csv.", call. = FALSE)

cat("[INFO] Reading expression matrices (this can take a minute)...\n")
expr_a <- data.table::fread(file.path(input_dir, EXPR_ANALYSIS_NAME), data.table = FALSE)
expr_o <- data.table::fread(file.path(input_dir, EXPR_OTHER_NAME), data.table = FALSE)
if (nrow(expr_a) != n_cells_all || nrow(expr_o) != n_cells_all) {
  stop("The expression matrices and ", FILENUMS_NAME, " have different numbers of rows (",
       nrow(expr_a), " / ", nrow(expr_o), " / ", n_cells_all, ").", call. = FALSE)
}
expr_all <- cbind(expr_a, expr_o)
rm(expr_a, expr_o)
missing_markers <- setdiff(plot_markers, colnames(expr_all))
if (length(missing_markers) > 0) {
  stop("These Plot = 1 channels from panel.csv are not in the expression matrices: ",
       paste(missing_markers, collapse = ", "), call. = FALSE)
}

# The published scripts back-transform EVERY plotted channel with 5 * sinh(x).
# Say so if the panel describes something different.
pl <- panel[panel$Plot == 1, ]
if (any(pl$Transform != "Asinh" | pl$Factor != ASINH_FACTOR)) {
  note_warning(paste0("panel.csv lists Plot = 1 channels that are not Asinh with factor ", ASINH_FACTOR, ": ",
                      paste(pl$Fixed_Param[pl$Transform != "Asinh" | pl$Factor != ASINH_FACTOR], collapse = ", "),
                      ". They are still back-transformed with 5 * sinh(x), as in the published scripts."))
}

keep_cell <- filenums %in% keep_files
if (length(CLUSTERS_TO_REMOVE) > 0) keep_cell <- keep_cell & !(clusters %in% CLUSTERS_TO_REMOVE)
if (use_cat) keep_cell <- keep_cell & (cat_all %in% cats_wanted)     # NA categories never match
cell_values <- as.matrix(expr_all[keep_cell, plot_markers, drop = FALSE])
if (!use_mean || MEAN_SCALE == "back_transformed") cell_values <- ASINH_FACTOR * sinh(cell_values)   # back-transform (percent mode: always)
cells <- data.table::as.data.table(cell_values)
rm(cell_values)
cells[, File_Num := filenums[keep_cell]]
cells[, Cluster  := clusters[keep_cell]]
cells[, Group    := if (by_cluster) Cluster else 0L]
cells[, ThreshCat := if (use_cat) cat_all[keep_cell] else ""]       # "" = no categories
rm(expr_all)
n_cells_used <- nrow(cells)
cat("[INFO] Cells: ", n_cells_all, " in the folder, ", n_cells_used, " analysed\n", sep = "")
if (n_cells_used == 0) stop("No cells are left after the sample and cluster choices.", call. = FALSE)
cells_per_file <- cells[, .N, by = File_Num][order(File_Num)]



#### Threshold from the baseline samples #######################################
thr_files_all <- md$File_Num[in_set & !is.na(md$TimePoint) & md$TimePoint == THRESHOLD_TIMEPOINT]
thr_files     <- intersect(thr_files_all, cells$File_Num)
if (length(thr_files) == 0) {
  stop("No baseline cells found: no analysed sample has TimePoint = ", THRESHOLD_TIMEPOINT,
       ". Check THRESHOLD_TIMEPOINT, SAMPLE_SET and CONDITIONS_TO_EXCLUDE.", call. = FALSE)
}
baseline_labels <- sort(unique(md$Sample_Condition[md$File_Num %in% thr_files]))
cat("[INFO] Baseline: ", paste(baseline_labels, collapse = ", "), " (", length(thr_files), " files)\n", sep = "")

base_cells  <- cells[File_Num %in% thr_files]
per_file_q  <- base_cells[, lapply(.SD, stats::quantile, probs = THRESHOLD_PROB, names = FALSE),
                          by = .(Group, ThreshCat, File_Num), .SDcols = plot_markers]
thr_by_group <- per_file_q[, lapply(.SD, mean), by = .(Group, ThreshCat), .SDcols = plot_markers]
data.table::setorder(thr_by_group, Group, ThreshCat)
thr_key     <- paste(thr_by_group$Group, thr_by_group$ThreshCat, sep = "|")
base_counts <- base_cells[, .(Baseline_Files = data.table::uniqueN(File_Num), Baseline_Cells = .N), by = .(Group, ThreshCat)]
groups_all  <- sort(unique(cells$Group))
all_keys    <- unique(cells[, .(Group, ThreshCat)])
no_base     <- all_keys[!(paste(Group, ThreshCat, sep = "|") %in% thr_key)]
if (nrow(no_base) > 0) {
  note_warning(paste0("No baseline cells for ", ifelse(use_cat, "cluster/category ", "cluster(s) "),
                      paste(if (use_cat) paste0(no_base$Group, "/", no_base$ThreshCat) else no_base$Group, collapse = ", "),
                      ": their percentages are NA."))
}
rm(base_cells, per_file_q)



#### Percent of cells above the threshold, per sample and cluster ##############
thr_mat <- as.matrix(thr_by_group[, plot_markers, with = FALSE])
pct <- cells[, {
  gi <- match(paste(Group[1], ThreshCat[1], sep = "|"), thr_key)
  tv <- if (is.na(gi)) rep(NA_real_, length(plot_markers)) else thr_mat[gi, ]
  out <- if (use_mean) {
    vapply(seq_along(plot_markers), function(j) mean(.SD[[j]]), numeric(1))                    # MEASURE = "mean_expression"
  } else {
    vapply(seq_along(plot_markers), function(j) mean(.SD[[j]] > tv[j]) * 100, numeric(1))      # percent above the threshold
  }
  c(list(N_Cells = .N), stats::setNames(as.list(out), plot_markers))
}, by = .(Group, ThreshCat, File_Num), .SDcols = plot_markers]
data.table::setnames(pct, "Group", "Cluster")

raw <- merge(pct, data.table::as.data.table(md_used), by = "File_Num", all.x = TRUE, sort = FALSE)
if (any(is.na(raw$File_Name))) stop("Internal check failed: results without metadata.", call. = FALSE)
id_cols <- c("File_Num", "File_Name", "Sample_Condition", "Reordered", "TimePoint", "Replicate", "Cluster", "ThreshCat", "N_Cells")
raw <- raw[, c(id_cols, plot_markers), with = FALSE]
data.table::setorder(raw, Reordered, Cluster, ThreshCat)
if (length(FILES_TO_EXCLUDE) > 0) raw <- raw[!(File_Name %in% FILES_TO_EXCLUDE)]
n_filled_rows <- 0L
if (FILL_EMPTY_CLUSTERS) {
  # as 12.1.1: complete(Cluster, File_Num) -> a row (all markers NA) for every sample x cluster with no cells
  want_files <- md_used$File_Num[!(md_used$File_Name %in% FILES_TO_EXCLUDE) & md_used$File_Num %in% cells$File_Num]
  all_pairs  <- data.table::CJ(File_Num = sort(unique(want_files)), Cluster = groups_all)
  absent     <- all_pairs[!raw[, .(File_Num, Cluster)], on = c("File_Num", "Cluster")]
  if (nrow(absent) > 0) {
    add_rows <- merge(absent, data.table::as.data.table(md_used), by = "File_Num", all.x = TRUE, sort = FALSE)
    add_rows[, ThreshCat := ""]
    add_rows[, N_Cells := 0L]
    for (m in plot_markers) data.table::set(add_rows, j = m, value = NA_real_)
    add_rows <- add_rows[, c(id_cols, plot_markers), with = FALSE]
    raw <- rbind(raw, add_rows)
    data.table::setorder(raw, Reordered, Cluster, ThreshCat)
    n_filled_rows <- nrow(add_rows)
    note_warning(paste0(n_filled_rows, " sample x cluster combination(s) have no cells; rows with NA values were added ",
                        "(FILL_EMPTY_CLUSTERS, as 12.1.1)."))
  }
}
condition_order <- unique(raw$Sample_Condition)          # plotting order of the labels



#### Fold change and log2 fold change #########################################
class_exclude <- c("H3K9ac_115In", "TrkB_s_158Gd", "yH2AX_165Ho", "CC3_173Yb", "TrkB_i_174Yb", "DNA_191Ir", "DNA_193Ir")
if (is.null(MARKERS_FOR_DERIVED)) {
  is_id <- !startsWith(plot_markers, "p") & !(plot_markers %in% class_exclude)
  derived_markers <- plot_markers[!is_id & !(plot_markers %in% c("DNA_191Ir", "DNA_193Ir"))]
  derived_rule <- "published rule (A4)"
} else {
  bad <- setdiff(MARKERS_FOR_DERIVED, plot_markers)
  if (length(bad) > 0) stop("MARKERS_FOR_DERIVED not in panel.csv Plot = 1: ", paste(bad, collapse = ", "), call. = FALSE)
  derived_markers <- MARKERS_FOR_DERIVED
  derived_rule <- "set by the user"
}

norm_conditions <- if (isTRUE(CTRL_PER_TIMEPOINT)) NORM_CONDITIONS else
  if (length(NORM_CONDITIONS) > 0) NORM_CONDITIONS else baseline_labels
not_found <- setdiff(norm_conditions, raw$Sample_Condition)
if (length(not_found) > 0) {
  stop("Control condition(s) not found among the analysed samples: ", paste(not_found, collapse = ", "), call. = FALSE)
}

key_of <- function(tp, cl, tc) paste(as.character(tp), cl, tc, sep = "|")
norm_tab <- raw[Sample_Condition %in% norm_conditions,
                c(list(NormTimepoint = mean(TimePoint, na.rm = TRUE)),
                  lapply(.SD, mean, na.rm = TRUE)),
                by = .(Sample_Condition, Cluster, ThreshCat), .SDcols = derived_markers]
if (anyDuplicated(key_of(norm_tab$NormTimepoint, norm_tab$Cluster, norm_tab$ThreshCat))) {
  stop("More than one control condition has the same TimePoint, so the fold-change control is ambiguous.\n",
       "Control conditions: ", paste(norm_conditions, collapse = ", "), call. = FALSE)
}

fc <- data.table::copy(raw[, c(setdiff(id_cols, "N_Cells"), derived_markers), with = FALSE])
fc_tp <- if (isTRUE(CTRL_PER_TIMEPOINT)) fc$TimePoint else rep(THRESHOLD_TIMEPOINT, nrow(fc))
nrow_match <- match(key_of(fc_tp, fc$Cluster, fc$ThreshCat), key_of(norm_tab$NormTimepoint, norm_tab$Cluster, norm_tab$ThreshCat))
if (anyNA(nrow_match)) {
  note_warning(paste0(sum(is.na(nrow_match)), " row(s) have no control to compare to (their fold change is NA)."))
}
for (m in derived_markers) {
  nv <- norm_tab[[m]][nrow_match]
  nv[!is.na(nv) & nv == 0] <- NA_real_              # avoid division by zero (as published)
  data.table::set(fc, j = m, value = fc[[m]] / nv)
}
log2fc <- data.table::copy(fc)
for (m in derived_markers) data.table::set(log2fc, j = m, value = log2(log2fc[[m]]))
n_neg_inf <- sum(vapply(derived_markers, function(m) sum(is.infinite(log2fc[[m]]) & log2fc[[m]] < 0), numeric(1)))
if (LOG2FC_INF_AS_NA) {
  for (m in derived_markers) {
    v <- log2fc[[m]]; v[is.infinite(v)] <- NA_real_
    data.table::set(log2fc, j = m, value = v)
  }
}

reps_avg <- function(tab, markers) {
  out <- tab[, c(list(TimePoint = mean(TimePoint, na.rm = TRUE), N_Replicates = .N),
                 lapply(.SD, mean, na.rm = TRUE)),
             by = .(Sample_Condition, Cluster, ThreshCat), .SDcols = markers]
  out[order(match(Sample_Condition, condition_order), Cluster, ThreshCat)]
}
raw_reps    <- reps_avg(raw, derived_markers)
fc_reps     <- reps_avg(fc, derived_markers)
log2fc_reps <- reps_avg(log2fc, derived_markers)

if (n_neg_inf > 0) {
  note_warning(paste0(n_neg_inf, " log2 fold-change value(s) are -Inf (", if (use_mean) "a mean of 0 in a sample" else "0 % above threshold in a sample", "). ",
                      if (LOG2FC_INF_AS_NA) "They were set to NA before averaging replicates (LOG2FC_INF_AS_NA, as 12.1.1)."
                      else "RepsAvg values for those markers are -Inf too (as in the 10.1 / 12.1 tables)."))
}

# Sample / cluster / category combinations with no cells (no row is written for them, as published)
expected_files  <- md_used$File_Num[!(md_used$File_Name %in% FILES_TO_EXCLUDE)]
expected_combos <- data.table::CJ(File_Num = sort(unique(expected_files)), Cluster = groups_all,
                                  ThreshCat = sort(unique(cells$ThreshCat)))
have_combos <- paste(raw$File_Num, raw$Cluster, raw$ThreshCat, sep = "|")
empty_combos <- expected_combos[!(paste(File_Num, Cluster, ThreshCat, sep = "|") %in% have_combos)]
if (nrow(empty_combos) > 0) {
  note_warning(paste0(nrow(empty_combos), " sample/cluster", if (use_cat) "/category" else "",
                      " combination(s) have no cells and no row (see run record): ",
                      paste(head(paste0(empty_combos$File_Num, "/", empty_combos$Cluster,
                                        if (use_cat) paste0("/", empty_combos$ThreshCat)), 6), collapse = ", "),
                      if (nrow(empty_combos) > 6) ", ..." else ""))
}



#### Save the tables ###########################################################
stage_label <- paste0(if (use_mean) "MeanExpr_" else "PercAboveThresh_", ANALYSIS_BY, "_", gsub("[^0-9]", "", paste(SAMPLE_SET, collapse = "")),
                      if (use_cat) paste0("_Cat", paste(cats_wanted, collapse = "-")) else "",
                      if (nzchar(OUTPUT_LABEL)) paste0("_", OUTPUT_LABEL) else "")
Stage_WD <- make_stage_dir(input_dir, stage_label)
prefix <- paste0("Expt", paste(SAMPLE_SET, collapse = ""), if (use_mean) "_MeanExpr_" else paste0("_PercAboveThresh_", THRESHOLD_PROB * 100, "thPerc_"))

out_tables <- list(
  RAWvalues         = raw,
  RAWvalues_RepsAvg = raw_reps,
  FC                = fc,
  FC_RepsAvg        = fc_reps,
  log2FC            = log2fc,
  log2FC_RepsAvg    = log2fc_reps)
if (!use_cat) out_tables <- lapply(out_tables, function(tb) tb[, setdiff(names(tb), "ThreshCat"), with = FALSE])
out_files <- character(0)
for (nm in names(out_tables)) {
  f <- file.path(Stage_WD, paste0(prefix, nm, ".csv"))
  data.table::fwrite(out_tables[[nm]], f)
  out_files <- c(out_files, f)
}

thr_long <- data.table::melt(thr_by_group, id.vars = c("Group", "ThreshCat"), variable.name = "Marker", value.name = "Threshold")
thr_long <- merge(thr_long, base_counts, by = c("Group", "ThreshCat"), all.x = TRUE)
data.table::setnames(thr_long, "Group", "Cluster")
data.table::setorder(thr_long, Cluster, ThreshCat, Marker)
if (!use_cat) thr_long[, ThreshCat := NULL]
if (!use_mean) {                                   # no threshold is used in mean mode
  thr_file <- file.path(Stage_WD, "Thresholds.csv")
  data.table::fwrite(thr_long, thr_file)
  out_files <- c(out_files, thr_file)
}

cat("[INFO] Tables written to:\n       ", Stage_WD, "\n", sep = "")



#### Run record ################################################################
Inputs <- list(
  Script_Version        = SCRIPT_VERSION,
  Completed_At          = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER          = input_dir,
  OUTPUT_FOLDER         = Stage_WD,
  Input_Provenance      = input_provenance,
  Settings = list(
    MEASURE = MEASURE, MEAN_SCALE = if (use_mean) MEAN_SCALE else NA_character_,
    ANALYSIS_BY = ANALYSIS_BY, CLUSTER_ROUND = CLUSTER_ROUND, CLUSTERS_MERGED = CLUSTERS_MERGED,
    Cluster_File = cluster_file,
    SAMPLE_SET_COLUMN = SAMPLE_SET_COLUMN, SAMPLE_SET = as.list(SAMPLE_SET),
    CONDITIONS_TO_EXCLUDE = as.list(CONDITIONS_TO_EXCLUDE), CLUSTERS_TO_REMOVE = as.list(CLUSTERS_TO_REMOVE),
    FILES_TO_EXCLUDE = as.list(FILES_TO_EXCLUDE),
    THRESHOLD_PROB = THRESHOLD_PROB, THRESHOLD_TIMEPOINT = THRESHOLD_TIMEPOINT,
    CTRL_PER_TIMEPOINT = CTRL_PER_TIMEPOINT, NORM_CONDITIONS_Used = as.list(norm_conditions),
    THRESHOLD_CATEGORY_FOLDER = THRESHOLD_CATEGORY_FOLDER, CATEGORIES_TO_ANALYSE = as.list(CATEGORIES_TO_ANALYSE),
    LOG2FC_INF_AS_NA = LOG2FC_INF_AS_NA, FILL_EMPTY_CLUSTERS = FILL_EMPTY_CLUSTERS,
    Quantile_Type = "R default (type 7)", Back_Transform = if (!use_mean) paste0(ASINH_FACTOR, " * sinh(x) on every Plot = 1 channel") else if (MEAN_SCALE == "back_transformed") paste0(ASINH_FACTOR, " * sinh(x), then mean") else "none (mean of the values as stored, as Fig 5F)"),
  Metadata_Mapping = list(
    Label_Column = LABEL_COLUMN_USED, Base_Treatment_Column = BASE_COLUMN_USED,
    Reordered = reordered_source),
  Baseline = list(Conditions = as.list(baseline_labels), Files = as.list(sort(thr_files))),
  Samples_Analysed       = as.list(md_used$File_Name),
  Cells_In_Folder        = n_cells_all,
  Cells_Analysed         = n_cells_used,
  Cells_Per_File         = lapply(seq_len(nrow(cells_per_file)), function(i)
    list(File_Num = cells_per_file$File_Num[i], Cells = cells_per_file$N[i])),
  Clusters_Analysed      = as.list(groups_all),
  Baseline_Counts        = lapply(seq_len(nrow(base_counts)), function(i)
    list(Cluster = base_counts$Group[i], Category = if (use_cat) base_counts$ThreshCat[i] else NA,
         Files = base_counts$Baseline_Files[i], Cells = base_counts$Baseline_Cells[i])),
  Filled_Empty_Rows      = n_filled_rows,
  Categories             = category_info,
  Empty_Combinations     = lapply(seq_len(nrow(empty_combos)), function(i)
    list(File_Num = empty_combos$File_Num[i], Cluster = empty_combos$Cluster[i],
         Category = if (use_cat) empty_combos$ThreshCat[i] else NA)),
  Markers_All            = as.list(plot_markers),
  Markers_Derived_Tables = as.list(derived_markers),
  Markers_Derived_Rule   = derived_rule,
  Files_Removed_From_Tables = as.list(intersect(FILES_TO_EXCLUDE, md$File_Name)),
  Random_Seed            = "not used in this script",
  Warnings               = as.list(warnings_log),
  Output_MD5             = file_md5(out_files),
  Session                = session_info_list(c("data.table", "jsonlite")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = input_dir)

cat("\n========================================================\n")
cat("[DONE] ", nrow(raw), " rows in RAWvalues; ", length(derived_markers), " markers in the fold-change tables\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("       Next: check against the published tables with Tools/Compare_GoldenFiles.R\n")
cat("========================================================\n")
