################################################################################
# 02_AssignThresholdCategories.R -- sort every cell into expression categories (Lo/Med/Hi ...)
################################################################################
# Author:   Jonathon Sewell (with Claude)
# Version:  2
# Created:  2026-10-23
# Edited:   2026-10-27
#
# PURPOSE
#   Gives every cell a category for one or two markers by cutting the marker's
#   expression at values YOU choose. Example: TrkB_i cut at 0.15, 1.4 and 3.7 gives
#   four levels (1 = lowest ... 4 = highest). With two markers the two levels are
#   written side by side: p75NTR level 1 and TrkB_i level 4 is category "14".
#   The categories are then used by 01_PercAboveThreshold.R (setting
#   THRESHOLD_CATEGORY_FOLDER) to ask "what percent of the TrkB_i-high cells ...".
#
#   This script replaces the DATA half of 13.1_ExtractByMarkerExpression_HiLo_v5.R.
#   The plots (histograms, biaxial plot, percentage lollipop plot) move to the
#   Visualization layer; this script writes the numbers those plots need.
#
# HOW TO USE
#   1. Set INPUT_FOLDER (the folder with expression_matrix_*.csv) and the markers
#      and cut values below. 2. Source the script. 3. Look at Category_Percentages.csv.
#
# WHAT YOU GET (in a new folder inside the input folder)
#   ThresholdCategory.csv         one row per cell, same order as the expression matrix
#   Category_Definitions.csv      the cut values and what each level means
#   Category_Percentages.csv      cells and percent of all cells in each category
#   Category_Cells_By_File.csv    cells per sample file in each category (n behind later results)
#   run_parameters.json           settings, checksums, warnings
#   Nothing from the input folder is copied (13.1 copied the whole dataset into
#   every output folder; later scripts now point at the input folder instead).
#
# CHANGE LOG FROM 13.1_ExtractByMarkerExpression_HiLo_v5.R
# A. Results
#  A1 [RESULT: none] Category rule is unchanged: a cell exactly on a cut value goes to the
#     HIGHER level (13.1 used between() on both sides, and the later level overwrote the
#     earlier one). Categories use the values in the expression matrix as they are (the
#     transformed values, not back-transformed). All cells in the folder are used.
#  A2 [BEHAVIOR] 13.1 asked for THRESHOLD.OPERATOR (">=" / "<=") but never used it; every
#     cut behaved as ">=". The setting is removed so it cannot suggest otherwise.
#  A3 [BEHAVIOR] 13.1 wrote expression_matrix_all_lo.csv / _hi.csv by looking for categories
#     called "lo" and "hi"; with numeric levels these files have always been empty (0 rows).
#     They are no longer written.
# C. Version 2 (2026-10-27), for the revision: how sensitive are the categories to the cut values?
#  C1 [RESULT: none at the defaults] New setting CUT_SHIFT (default 0). A value of 0.2 moves every
#     user-defined cut value 20 % up, -0.2 moves them 20 % down. Use it to re-make the categories with
#     slightly different cuts and ask whether a conclusion still holds. At 0 the output is identical to v1.
#  C2 New file Histogram_Data.csv: counts of cells per bin for each marker, so the histogram with the cut
#     values drawn on it can be plotted without reading the whole matrix again.
#  C3 The run record and Category_Definitions.csv state the shift used; the output folder name carries it.
# B. Efficiency / readability / reproducibility
#  B1 Only the marker columns are read (13.1 read and re-wrote every matrix).
#  B2 Category_Cells_By_File.csv and a run record are new.
#  B3 No random numbers are used.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
INPUT_FOLDER <- ""        # folder with the expression matrices ("" = choose in a pop-up)

# One marker (x) or two markers (x and y). Names must be x and y.
SELECTED_MARKERS <- c(x = "p75NTR_161Dy",
                      y = "TrkB_i_174Yb")      # delete the y line to use one marker

# Cut values for each marker, in increasing order: n cut values give n + 1 levels.
CUT_VALUES <- list(x = c(1.45),
                   y = c(0.15, 1.4, 3.7))
# Look at the marker's histogram (or another justification) before choosing cut values,
# then change this to TRUE.
CUTS_JUSTIFIED <- TRUE

# Sensitivity check: move ALL user-defined cut values by this fraction (0 = as chosen above).
# 0.2 = 20 % higher, -0.2 = 20 % lower. Leave at 0 for the published categories.
# Only for positive cut values; not used with CUTS_BY_QUANTILE.
CUT_SHIFT <- 0

# Alternative: cut at quantiles of all cells instead (used for every marker).
CUTS_BY_QUANTILE <- FALSE
CUT_QUANTILES    <- c(1/3, 2/3)   # e.g. c(1/3, 2/3) = low / middle / high thirds

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "02_AssignThresholdCategories.R"
SCRIPT_VERSION           <- 2
MODULE_BASE_NAME         <- "02_AssignThresholdCategories"
REQUIRED_HELPERS_VERSION <- 7
EXPR_ANALYSIS_NAME       <- "expression_matrix_analysis.csv"
EXPR_OTHER_NAME          <- "expression_matrix_other.csv"
FILENUMS_NAME            <- "filenums.csv"



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
cat("[START] Expression categories:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

warnings_log <- character(0)
note_warning <- function(txt) {
  warnings_log <<- c(warnings_log, txt)
  cat("[WARNING] ", txt, "\n", sep = "")
}

# Check the settings before doing any work
axes <- names(SELECTED_MARKERS)
if (is.null(axes) || !all(axes %in% c("x", "y")) || anyDuplicated(axes) || length(axes) == 0 ||
    !("x" %in% axes)) {
  stop("SELECTED_MARKERS must be named x (and optionally y), e.g. c(x = \"p75NTR_161Dy\", y = \"TrkB_i_174Yb\").",
       call. = FALSE)
}
if (anyDuplicated(unname(SELECTED_MARKERS))) stop("SELECTED_MARKERS lists the same marker twice.", call. = FALSE)
if (!isTRUE(CUTS_JUSTIFIED)) {
  stop("Look at the histogram of each marker (or use another justification) before choosing cut values,\n",
       "then set CUTS_JUSTIFIED <- TRUE at the top of this script.", call. = FALSE)
}
if (!is.numeric(CUT_SHIFT) || length(CUT_SHIFT) != 1 || is.na(CUT_SHIFT) || CUT_SHIFT <= -1) {
  stop("CUT_SHIFT must be one number greater than -1 (0 = no shift, 0.2 = 20 % up, -0.2 = 20 % down).", call. = FALSE)
}
if (CUT_SHIFT != 0 && CUTS_BY_QUANTILE) {
  stop("CUT_SHIFT applies to user-defined CUT_VALUES, not to quantile cuts. Set CUT_SHIFT <- 0.", call. = FALSE)
}
if (CUTS_BY_QUANTILE) {
  if (!is.numeric(CUT_QUANTILES) || length(CUT_QUANTILES) == 0 || any(CUT_QUANTILES <= 0 | CUT_QUANTILES >= 1) ||
      is.unsorted(CUT_QUANTILES, strictly = TRUE)) {
    stop("CUT_QUANTILES must be increasing numbers between 0 and 1, e.g. c(1/3, 2/3).", call. = FALSE)
  }
} else {
  for (ax in axes) {
    cv <- CUT_VALUES[[ax]]
    if (is.null(cv) || !is.numeric(cv) || length(cv) == 0 || anyNA(cv) || is.unsorted(cv, strictly = TRUE)) {
      stop("CUT_VALUES$", ax, " must be one or more numbers in increasing order.", call. = FALSE)
    }
    if (CUT_SHIFT != 0 && any(cv <= 0)) {
      stop("CUT_SHIFT needs positive cut values; CUT_VALUES$", ax, " has zero or negative values.", call. = FALSE)
    }
  }
}



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
    pick <- ask_choice("Which data folder should be used?", choices, preselect = choices[1])
    input_dir <- if (here_has_data && pick == basename(Original_WD)) Original_WD else dirs[basename(dirs) == pick]
    input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  }
}
cat("[INFO] Input folder:\n       ", input_dir, "\n", sep = "")

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



#### Read the marker columns ###################################################
# Find which matrix holds each marker, and read only those columns
hdr_a <- colnames(data.table::fread(file.path(input_dir, EXPR_ANALYSIS_NAME), nrows = 0))
hdr_o <- colnames(data.table::fread(file.path(input_dir, EXPR_OTHER_NAME), nrows = 0))
find_marker_file <- function(mk) {
  if (mk %in% hdr_a) EXPR_ANALYSIS_NAME else if (mk %in% hdr_o) EXPR_OTHER_NAME else NA_character_
}
where <- vapply(SELECTED_MARKERS, find_marker_file, character(1))
if (anyNA(where)) {
  stop("These markers were not found in the expression matrices: ",
       paste(SELECTED_MARKERS[is.na(where)], collapse = ", "),
       "\nThe columns are named like Tuj1_89Y (see panel.csv, column Fixed_Param).", call. = FALSE)
}
vals <- lapply(axes, function(ax) {
  mk <- SELECTED_MARKERS[[ax]]
  d <- data.table::fread(file.path(input_dir, where[[ax]]), select = mk, data.table = FALSE)
  as.numeric(d[[1]])
})
names(vals) <- axes
n_cells <- length(vals[[1]])
if (any(vapply(vals, length, numeric(1)) != n_cells)) stop("The marker columns have different lengths.", call. = FALSE)

filenums_path <- file.path(input_dir, FILENUMS_NAME)
filenums <- if (file.exists(filenums_path)) data.table::fread(filenums_path, data.table = FALSE)[[1]] else NULL
if (!is.null(filenums) && length(filenums) != n_cells) {
  stop(FILENUMS_NAME, " has ", length(filenums), " rows but the expression matrices have ", n_cells, ".", call. = FALSE)
}
cat("[INFO] Cells: ", n_cells, "\n", sep = "")
for (ax in axes) {
  if (anyNA(vals[[ax]])) {
    note_warning(paste0(sum(is.na(vals[[ax]])), " cell(s) have a missing value for ", SELECTED_MARKERS[[ax]],
                        "; their category is NA."))
  }
}



#### Decide the cut values and assign levels ###################################
cuts <- lapply(axes, function(ax) {
  if (CUTS_BY_QUANTILE) {
    unname(stats::quantile(vals[[ax]], probs = CUT_QUANTILES, na.rm = TRUE, type = 7))
  } else {
    as.numeric(CUT_VALUES[[ax]]) * (1 + CUT_SHIFT)
  }
})
names(cuts) <- axes

levels_by_axis <- lapply(axes, function(ax) {
  # level = 1 + the number of cut values that are <= the cell's value (a cell ON a cut goes up)
  findInterval(vals[[ax]], cuts[[ax]], left.open = FALSE) + 1L
})
names(levels_by_axis) <- axes

for (ax in axes) {
  r <- range(vals[[ax]], na.rm = TRUE)
  out_of_range <- cuts[[ax]][cuts[[ax]] <= r[1] | cuts[[ax]] > r[2]]
  if (length(out_of_range) > 0) {
    note_warning(paste0("Cut value(s) ", paste(out_of_range, collapse = ", "), " for ", SELECTED_MARKERS[[ax]],
                        " are outside the data range (", signif(r[1], 4), " to ", signif(r[2], 4),
                        "); at least one level is empty."))
  }
  if ((length(cuts[[ax]]) + 1L) > 9L) {
    stop("More than 9 levels for ", SELECTED_MARKERS[[ax]], ": the written category (e.g. \"14\") would be ambiguous.",
         call. = FALSE)
  }
}

category <- if (length(axes) == 1) {
  as.character(levels_by_axis[[1]])
} else {
  ifelse(is.na(levels_by_axis[["x"]]) | is.na(levels_by_axis[["y"]]), NA_character_,
         paste0(levels_by_axis[["x"]], levels_by_axis[["y"]]))
}
cat_table <- table(category, useNA = "no")
cat("[INFO] Categories (cells): ", paste0(names(cat_table), " = ", as.integer(cat_table), collapse = "; "), "\n", sep = "")
empty_possible <- setdiff(
  if (length(axes) == 1) as.character(seq_len(length(cuts[[1]]) + 1L)) else
    as.vector(outer(seq_len(length(cuts[["x"]]) + 1L), seq_len(length(cuts[["y"]]) + 1L), paste0)),
  names(cat_table))
if (length(empty_possible) > 0) {
  note_warning(paste0("No cells in category: ", paste(sort(empty_possible), collapse = ", "),
                      ". Later results have no rows for these categories."))
}



#### Summary tables ############################################################
# Definitions: what each level means
defs <- do.call(rbind, lapply(axes, function(ax) {
  cv <- cuts[[ax]]
  data.frame(Axis = ax, Marker = SELECTED_MARKERS[[ax]], Level = seq_len(length(cv) + 1L),
             Lower_Cut = c(NA_real_, cv), Upper_Cut = c(cv, NA_real_),
             Rule = c(paste0("< ", cv[1]),
                      if (length(cv) > 1) paste0(">= ", cv[-length(cv)], " and < ", cv[-1]),
                      paste0(">= ", cv[length(cv)])),
             Cells = as.integer(table(factor(levels_by_axis[[ax]], levels = seq_len(length(cv) + 1L)))),
             Percent_Of_Cells = as.numeric(table(factor(levels_by_axis[[ax]], levels = seq_len(length(cv) + 1L)))) /
               n_cells * 100,
             stringsAsFactors = FALSE, row.names = NULL)
}))
defs$Cuts_Chosen_By <- if (CUTS_BY_QUANTILE) paste0("quantiles of all cells: ", paste(signif(CUT_QUANTILES, 4), collapse = ", ")) else
  if (CUT_SHIFT != 0) paste0("user-defined values shifted by ", CUT_SHIFT * 100, " %") else "user-defined values"

perc <- data.frame(ThresholdCategory = names(cat_table), Cells = as.integer(cat_table),
                   Percent_Of_Cells = as.integer(cat_table) / n_cells * 100, stringsAsFactors = FALSE)

by_file <- NULL
if (!is.null(filenums)) {
  by_file <- as.data.frame(table(File_Num = filenums, ThresholdCategory = category), stringsAsFactors = FALSE)
  by_file <- by_file[by_file$Freq > 0, , drop = FALSE]
  by_file$File_Num <- as.integer(by_file$File_Num)
  names(by_file)[names(by_file) == "Freq"] <- "Cells"
  by_file <- by_file[order(by_file$File_Num, by_file$ThresholdCategory), , drop = FALSE]
}



#### Save ######################################################################
label <- paste0("ThreshCat_", paste(vapply(axes, function(ax) {
  cv <- cuts[[ax]]
  paste0(SELECTED_MARKERS[[ax]], "-", paste0("Gr", round(cv, 2), collapse = "-"))
}, character(1)), collapse = "_"))
if (CUT_SHIFT != 0) label <- paste0(label, "_Shift", sprintf("%+g", CUT_SHIFT * 100), "pct")
Stage_WD <- make_stage_dir(input_dir, label)

# Histogram counts (100 equal bins over each marker's range) for plotting the cuts on the distribution
hist_data <- do.call(rbind, lapply(axes, function(ax) {
  v <- vals[[ax]][is.finite(vals[[ax]])]
  br <- seq(min(v), max(v), length.out = 101)
  cnt <- graphics::hist(v, breaks = br, plot = FALSE, include.lowest = TRUE)$counts
  data.frame(Axis = ax, Marker = SELECTED_MARKERS[[ax]], Bin_Lower = br[-101], Bin_Upper = br[-1],
             Cells = as.integer(cnt), stringsAsFactors = FALSE)
}))

out_files <- character(0)
f1 <- file.path(Stage_WD, "ThresholdCategory.csv")
data.table::fwrite(data.frame(ThresholdCategory = category), f1);              out_files <- c(out_files, f1)
f2 <- file.path(Stage_WD, "Category_Definitions.csv")
utils::write.csv(defs, f2, row.names = FALSE);                                 out_files <- c(out_files, f2)
f3 <- file.path(Stage_WD, "Category_Percentages.csv")
utils::write.csv(perc, f3, row.names = FALSE);                                 out_files <- c(out_files, f3)
f5 <- file.path(Stage_WD, "Histogram_Data.csv")
utils::write.csv(hist_data, f5, row.names = FALSE);                            out_files <- c(out_files, f5)
if (!is.null(by_file)) {
  f4 <- file.path(Stage_WD, "Category_Cells_By_File.csv")
  utils::write.csv(by_file, f4, row.names = FALSE);                            out_files <- c(out_files, f4)
}
cat("[INFO] Files written to:\n       ", Stage_WD, "\n", sep = "")



#### Run record ################################################################
Inputs <- list(
  Script_Version   = SCRIPT_VERSION,
  Completed_At     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER     = input_dir,
  OUTPUT_FOLDER    = Stage_WD,
  Input_Provenance = input_provenance,
  Settings = list(
    SELECTED_MARKERS = as.list(SELECTED_MARKERS),
    CUTS_BY_QUANTILE = CUTS_BY_QUANTILE,
    CUT_QUANTILES    = if (CUTS_BY_QUANTILE) as.list(CUT_QUANTILES) else list(),
    Cut_Values_Used  = lapply(cuts, as.list),
    CUTS_JUSTIFIED   = CUTS_JUSTIFIED,
    CUT_SHIFT        = CUT_SHIFT,
    Rule             = "level = 1 + number of cut values <= the cell's value (a cell on a cut goes to the higher level)"),
  Cells            = n_cells,
  Marker_Files     = as.list(where),
  Category_Cells   = as.list(stats::setNames(as.integer(cat_table), names(cat_table))),
  Random_Seed      = "not used in this script",
  Warnings         = as.list(warnings_log),
  Output_MD5       = file_md5(out_files),
  Session          = session_info_list(c("data.table", "jsonlite")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = input_dir)

cat("\n========================================================\n")
cat("[DONE] ", n_cells, " cells sorted into ", length(cat_table), " categories\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("       Next: in 01_PercAboveThreshold.R set THRESHOLD_CATEGORY_FOLDER to:\n       ", Stage_WD, "\n", sep = "")
cat("========================================================\n")
