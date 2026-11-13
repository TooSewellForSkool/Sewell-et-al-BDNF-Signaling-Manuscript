################################################################################
# 03_Statistics.R
################################################################################
# Layer:    Analysis -- statistics (no plots; the Visualization layer reads the tables written here)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-25
# Edited:   2026-10-25
# Version:  2   (1 = the PUBLISHED statistics only (golden-checked, runs P1-P4; saved in _superseded).
#                2 = v1 + design summary with litter composition + assumption diagnostics.
#                    The p-values are the same as v1.
#                    Planned: 3 = suggested test, ANOVA / post-hoc tests and multiple-testing adjustment.)
#
# Purpose:
# * Take one table made by 01_PercAboveThreshold.R (FC, log2FC or RAWvalues; one row per sample
#   and cluster) and compare every group at every time point with a REFERENCE group
#   (published: NoTreat at 0 hr), separately for every marker.
# * Write the p-values as tables that the plotting scripts read, instead of calculating
#   them again inside the plotting code.
#
# This version does exactly what the published line-graph scripts did
# (LineGraphAnalysis_..._(stats table).R, function MakeStats_ConditionByTime):
#     for every marker, and every (time point, group):   t.test(group values, reference values)
#     with R's default settings, unpaired, NO adjustment for multiple testing.
#
# ==============================================================================
# CHANGE LOG: published scripts -> v1
# ==============================================================================
# Convention:  [RESULT]   = can change a p-value
#              [BEHAVIOR] = changes how the script behaves, not the numbers
#
# A. Needs your attention (nothing here is applied silently)
#  A1 [RESULT: none expected] Same calculation as the published scripts. The golden check
#     (Tools/Compare_GoldenFiles.R, run P1 in README_Analysis.txt) must reproduce the published
#     Supplementary Table 3 before this is used for anything new.
#  A2 [WORDING, not a result] The published figure legends / Supplementary Tables say
#     "unpaired Student t-test". The scripts call R's t.test(..., paired = FALSE) with the
#     default var.equal = FALSE, which is WELCH's t-test (unequal variances). TEST = "welch" is
#     therefore the published method; "student" is offered but was NOT what was run.
#  A3 [RESULT, declared] No adjustment for multiple testing (published). Not changed in v1.
#  A4 [RESULT, declared] Published behaviour that is kept as it was, and that you may want to change on purpose:
#       (a) With COMPARE_BY = "Cluster" or "ThreshCat" the reference values are POOLED over
#           all clusters (categories) in the table (REFERENCE_SCOPE = "pooled"). "matched" uses only
#           the reference values of the same cluster (category).
#       (b) The reference group is also tested against itself (time point 0). When a group's mean
#           equals the reference mean (fold-change tables are normalised so the 0 hr mean is 1) the
#           p-value is exactly 1.
#       (c) A group with fewer than 2 values gets p = NA.
#       (d) Samples are treated as independent. Replicate numbers are NOT litters (see
#           Inventory/Decisions_Log.md, 2026-10-25), so no pairing or blocking is done.
#  A5 [BEHAVIOR] The published code took the table, the clusters/categories and the treatments to
#     compare from several settings spread over a long plotting script. They are the settings below.
#     The script stops (instead of pooling silently) when the selection would mix clusters,
#     categories or treatments inside one compared group; the published panels never did.
#  A6 [BEHAVIOR] New columns in the output: N_Group, N_Reference, Mean_Group, Mean_Reference, Test,
#     so the sample size and the test are stored next to every p-value.
#
# D. Bugfix 2026-10-28 (v2.1; Script_Version stays 2; no p-value, diagnostic or output changes for any data that ran before)
#  D1 [BUGFIX] Diagnostics (RUN_DIAGNOSTICS) built the group labels of a family with paste0("t", time points). When a
#     treatment had NO rows in a cluster/category, paste0("t", numeric(0)) returned "t" (one label too many) and the
#     diagnostics stopped with a length error. Such families now get an empty label vector and are handled as before.
#     Found while fixing the same line in 04_TimecourseStatistics (v1.2); only small cluster/category tables are affected.
#
# C. Changes in v2 (design summary + diagnostics; the p-values do NOT change)
#  C1 [BEHAVIOR] Statistics_design.csv now also lists the LITTERS behind every group (N_Litters, Litters,
#     Max_Samples_Per_Litter), read from Litter_Map.csv (next to this script). Litter_Map.csv was built from
#     the 'Date Collected' column of the sample-set PDF and covers the 45 samples of the published tables;
#     CHECK IT against your records. Statistics_litters.csv lists the litter of every sample used.
#     Samples from the same litter are NOT independent; the tests still treat them as independent
#     (published, and as the reviewer assumed). The warning below tells you when that applies.
#  C2 [BEHAVIOR] Statistics_diagnostics.csv: for every marker and every "family" (one treatment over its time
#     points together with the reference, in one cluster/category) it reports
#       * Shapiro-Wilk test on the residuals (value minus its group mean, all groups together),
#       * Brown-Forsythe test (equal spread across groups; Levene's test around the median),
#       * variance ratio (largest / smallest group variance), skewness of the residuals,
#       * Spearman correlation between group mean and group SD (does the spread grow with the mean?),
#       * the same Shapiro / Brown-Forsythe checks on log2(values) when STAT_TABLE = "FC"
#         (helps to decide between FC and log2FC for the tests).
#     NOTHING IS DECIDED FROM THESE: they are reported next to the tests, not used to choose a test.
#     With 2-4 samples per group these tests have very little power, a 'not significant' result does
#     not show that the data are normal. Read them as a screen, together with the design table.
#  C3 [BEHAVIOR] New settings RUN_DIAGNOSTICS and LITTER_MAP. The run record lists the diagnostic summary.

# B. Efficiency / readability / reproducibility
#  B1 One run record (run_parameters.json) with all settings, the input table (checksum),
#     the number of comparisons, warnings, output checksums and the session.
#  B2 Three tables: Statistics_long.csv (one row per marker x group x time point),
#     Statistics_pvalues_wide.csv (the layout of the Supplementary Tables), Statistics_design.csv (n per group).
#  B3 No random numbers are used, so no seed is needed.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
INPUT_FOLDER   <- ""            # a folder made by 01_PercAboveThreshold.R ("" = choose in a pop-up)

STAT_TABLE     <- "RAWvalues"          # which table is tested: "FC", "log2FC" or "RAWvalues"
                                # (the *_RepsAvg tables are averages and are never tested)
TEST           <- "welch"       # "welch" (published; R's t.test default), "student", or "wilcox"

# What is compared with what
COMPARE_BY          <- "Treatment"   # the groups that are compared with the reference:
                                     #   "Treatment" (one cluster/category at a time), "Cluster" or "ThreshCat"
REFERENCE_TREATMENT <- "NoTreat"     # the reference group = this treatment ...
REFERENCE_TIMEPOINT <- 0             # ... at this time point (hours)
REFERENCE_SCOPE     <- "pooled"      # only for COMPARE_BY "Cluster"/"ThreshCat": "pooled" (published) or "matched"

# Which part of the table is used (c() = everything in the table)
TREATMENTS     <- c("NoTreat", "CompMed")           # e.g. c("NoTreat", "ContDepr", "CompMed", "BDNF")  (must include REFERENCE_TREATMENT)
CLUSTERS       <- c()           # e.g. c(1, 2, 4, 7, 13)
CATEGORIES     <- c()           # e.g. c("22")   (only for tables made with expression categories)
MARKERS        <- c("pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La",
                    "pPLK1_113In", "yH2AX_165Ho", "TrkB_s_158Gd",
                    "pATM_159Tb", "pPLCG-2_144Nd", "pSrc_154Sm",
                    "H3K9ac_115In", "p-p38_142Nd", "CC3_173Yb",
                    "pSTAT3_153Eu", "pSTAT5_147Sm", "p-p53_151Eu",
                    "p-cJUN_162Dy", "pNFKB_166Er", "pTrkB_170Er", 
                    "p75NTR_161Dy", "TrkB_i_174Yb")           # e.g. c("pERK_171Yb", "pAKt_168Er", "TrkB_s_158Gd")   c() = every marker in the table

# c("pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La",
#   "pPLK1_113In", "yH2AX_165Ho", "TrkB_s_158Gd",
#   "pATM_159Tb", "pPLCG-2_144Nd", "pSrc_154Sm",
#   "H3K9ac_115In", "p-p38_142Nd", "CC3_173Yb",
#   "pSTAT3_153Eu", "pSTAT5_147Sm", "p-p53_151Eu",
#   "p-cJUN_162Dy", "pNFKB_166Er", "pTrkB_170Er", 
#   "p75NTR_161Dy", "TrkB_i_174Yb")



OUTPUT_LABEL   <- "FigS5_onlySignalingMarkers"            # optional text added to the output folder name

RUN_DIAGNOSTICS <- TRUE         # TRUE: also write Statistics_diagnostics.csv (normality / equal-spread checks)
LITTER_MAP      <- ""           # "" = Litter_Map.csv in the ScriptsToUse/02_Analysis folder; or a full file path
                                #      (columns: Sample_Condition, Replicate, Litter). "none" = skip the litter columns.

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "03_Statistics.R"
SCRIPT_VERSION           <- 2
MODULE_BASE_NAME         <- "03_Statistics"
REQUIRED_HELPERS_VERSION <- 7
UPSTREAM_MODULE          <- "01_PercAboveThreshold"
ID_COLUMNS <- c("File_Num", "File_Name", "Sample_Condition", "Reordered", "TimePoint",
                "Replicate", "Cluster", "ThreshCat", "N_Cells", "N_Replicates")



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
cat("[START] Statistics against a reference group:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

warnings_log <- character(0)
note_warning <- function(txt) {
  warnings_log <<- c(warnings_log, txt)
  cat("[WARNING] ", txt, "\n", sep = "")
}

# Check the settings before doing any work
one_of <- function(value, allowed, name) {
  if (!is.character(value) || length(value) != 1 || !value %in% allowed) {
    stop(name, " must be one of: ", paste0("\"", allowed, "\"", collapse = ", "), ".", call. = FALSE)
  }
}
one_of(STAT_TABLE,       c("FC", "log2FC", "RAWvalues"), "STAT_TABLE")
one_of(TEST,             c("welch", "student", "wilcox"), "TEST")
one_of(COMPARE_BY,       c("Treatment", "Cluster", "ThreshCat"), "COMPARE_BY")
one_of(REFERENCE_SCOPE,  c("pooled", "matched"), "REFERENCE_SCOPE")
if (!is.logical(RUN_DIAGNOSTICS) || length(RUN_DIAGNOSTICS) != 1 || is.na(RUN_DIAGNOSTICS)) {
  stop("RUN_DIAGNOSTICS must be TRUE or FALSE.", call. = FALSE)
}
if (!is.numeric(REFERENCE_TIMEPOINT) || length(REFERENCE_TIMEPOINT) != 1) {
  stop("REFERENCE_TIMEPOINT must be a single number (hours), e.g. 0.", call. = FALSE)
}
if (TEST == "student") {
  note_warning(paste0("TEST = \"student\": the published statistics used Welch's t-test (TEST = \"welch\"), ",
                      "not Student's. Results will differ slightly from the published p-values."))
}
if (length(TREATMENTS) > 0 && !REFERENCE_TREATMENT %in% TREATMENTS) {
  stop("TREATMENTS must include the reference treatment (", REFERENCE_TREATMENT, ").", call. = FALSE)
}

table_pattern <- switch(STAT_TABLE, FC = "_FC\\.csv$", log2FC = "_log2FC\\.csv$", RAWvalues = "_RAWvalues\\.csv$")



#### Choose the input folder ###################################################
folder_has_table <- function(d) length(list.files(d, pattern = table_pattern)) > 0
if (nzchar(INPUT_FOLDER)) {
  if (!dir.exists(INPUT_FOLDER)) stop("INPUT_FOLDER not found:\n  ", INPUT_FOLDER, call. = FALSE)
  input_dir <- normalizePath(INPUT_FOLDER, winslash = "/", mustWork = TRUE)
} else {
  Original_WD <- get_working_dir("project")
  here_has_table <- folder_has_table(Original_WD)
  dirs <- list.dirs(Original_WD, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[!grepl("^[_.]", basename(dirs))]
  dirs <- dirs[vapply(dirs, folder_has_table, logical(1))]
  dirs <- sort(dirs, decreasing = TRUE)                       # newest first
  if (here_has_table && length(dirs) == 0) {
    input_dir <- Original_WD
  } else if (length(dirs) == 0) {
    stop("No folder with a ", STAT_TABLE, " table (made by ", UPSTREAM_MODULE, ".R) was found in:\n  ",
         Original_WD, "\nChoose the right folder, or set INPUT_FOLDER at the top of this script.", call. = FALSE)
  } else {
    choices <- c(if (here_has_table) basename(Original_WD), basename(dirs))
    pick <- ask_choice("Which table folder should be analysed?", choices, preselect = choices[1])
    input_dir <- if (here_has_table && pick == basename(Original_WD)) Original_WD else dirs[basename(dirs) == pick]
    input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  }
}
cat("[INFO] Input folder:\n       ", input_dir, "\n", sep = "")

# Sequence check: the table must come from 01_PercAboveThreshold.R
require_modules(input_dir, UPSTREAM_MODULE, how_to_fix = "Run 01_PercAboveThreshold.R first.")
input_record <- read_run_record(input_dir)

table_files <- list.files(input_dir, pattern = table_pattern, full.names = TRUE)
if (length(table_files) != 1) {
  stop("Expected exactly one ", STAT_TABLE, " table in the input folder, found ", length(table_files), ":\n  ",
       if (length(table_files) > 0) paste(basename(table_files), collapse = "\n  ") else "(none)",
       "\nFolder: ", input_dir, call. = FALSE)
}
table_path <- table_files[1]
cat("[INFO] Table: ", basename(table_path), "\n", sep = "")



#### Read the table and select the part to test ################################
tb <- data.table::fread(table_path, na.strings = c("", "NA"), data.table = TRUE)
if (names(tb)[1] %in% c("", "V1") && ncol(tb) > 1 &&
    identical(as.integer(tb[[1]]), seq_len(nrow(tb)))) tb[, (1) := NULL]   # row-number column
for (cn in c("Sample_Condition", "TimePoint", "Cluster")) {
  if (!cn %in% names(tb)) {
    stop("The table has no column '", cn, "'. Was it made by 01_PercAboveThreshold.R?\nFile: ", table_path, call. = FALSE)
  }
}
if (!"ThreshCat" %in% names(tb)) tb[, ThreshCat := ""]
tb[, ThreshCat := as.character(ThreshCat)]
tb[is.na(ThreshCat), ThreshCat := ""]
tb[, Treatment := sub("_[^_]+$", "", Sample_Condition)]
n_rows_table <- nrow(tb)

if (length(TREATMENTS) > 0) {
  missing_t <- setdiff(TREATMENTS, unique(tb$Treatment))
  if (length(missing_t) > 0) {
    stop("TREATMENTS not found in the table: ", paste(missing_t, collapse = ", "),
         "\nTreatments in the table: ", paste(sort(unique(tb$Treatment)), collapse = ", "), call. = FALSE)
  }
  tb <- tb[Treatment %in% TREATMENTS]
}
if (length(CLUSTERS) > 0) {
  missing_c <- setdiff(as.character(CLUSTERS), as.character(unique(tb$Cluster)))
  if (length(missing_c) > 0) {
    stop("CLUSTERS not found in the table: ", paste(missing_c, collapse = ", "),
         "\nClusters in the table: ", paste(sort(unique(tb$Cluster)), collapse = ", "), call. = FALSE)
  }
  tb <- tb[as.character(Cluster) %in% as.character(CLUSTERS)]
}
if (length(CATEGORIES) > 0) {
  if (all(tb$ThreshCat == "")) stop("CATEGORIES is set but the table has no expression categories.", call. = FALSE)
  missing_k <- setdiff(as.character(CATEGORIES), unique(tb$ThreshCat))
  if (length(missing_k) > 0) {
    stop("CATEGORIES not found in the table: ", paste(missing_k, collapse = ", "),
         "\nCategories in the table: ", paste(sort(unique(tb$ThreshCat)), collapse = ", "), call. = FALSE)
  }
  tb <- tb[ThreshCat %in% as.character(CATEGORIES)]
}
if (COMPARE_BY == "ThreshCat" && all(tb$ThreshCat == "")) {
  stop("COMPARE_BY is \"ThreshCat\" but the table has no expression categories.", call. = FALSE)
}

# Every compared group must be ONE thing (the published panels never mixed clusters/categories/treatments)
n_clusters <- data.table::uniqueN(tb$Cluster)
n_cats     <- data.table::uniqueN(tb$ThreshCat)
if (COMPARE_BY != "Cluster" && n_clusters > 1) {
  stop("The selection holds ", n_clusters, " clusters (", paste(sort(unique(tb$Cluster)), collapse = ", "), ") ",
       "but COMPARE_BY is \"", COMPARE_BY, "\".\nSet CLUSTERS to ONE cluster (run once per cluster), ",
       "or use COMPARE_BY = \"Cluster\".", call. = FALSE)
}
if (COMPARE_BY != "ThreshCat" && n_cats > 1) {
  stop("The selection holds ", n_cats, " categories (", paste(sort(unique(tb$ThreshCat)), collapse = ", "), ") ",
       "but COMPARE_BY is \"", COMPARE_BY, "\".\nSet CATEGORIES to ONE category (run once per category), ",
       "or use COMPARE_BY = \"ThreshCat\".", call. = FALSE)
}
non_ref_treatments <- setdiff(unique(tb$Treatment), REFERENCE_TREATMENT)
if (COMPARE_BY != "Treatment" && length(non_ref_treatments) > 1) {
  stop("COMPARE_BY is \"", COMPARE_BY, "\" but the selection holds ", length(non_ref_treatments),
       " treatments besides the reference (", paste(non_ref_treatments, collapse = ", "), ").\n",
       "Their samples would be mixed inside one group. Set TREATMENTS to the reference plus ONE treatment.", call. = FALSE)
}

tb[, Group := as.character(get(COMPARE_BY))]
is_ref <- tb$Treatment == REFERENCE_TREATMENT & abs(tb$TimePoint - REFERENCE_TIMEPOINT) < 1e-9
ref_tb <- tb[is_ref]
if (nrow(ref_tb) < 2) {
  stop("The reference group (", REFERENCE_TREATMENT, " at ", REFERENCE_TIMEPOINT, " hr) has ", nrow(ref_tb),
       " row(s) in the selection; at least 2 are needed.", call. = FALSE)
}

marker_cols <- setdiff(names(tb), c(ID_COLUMNS, "Treatment", "Group"))
markers <- if (length(MARKERS) == 0) marker_cols else MARKERS
if (length(setdiff(markers, marker_cols)) > 0) {
  stop("MARKERS not found in the table: ", paste(setdiff(markers, marker_cols), collapse = ", "), call. = FALSE)
}
cat("[INFO] ", nrow(tb), " of ", n_rows_table, " table rows selected; ", length(markers), " markers; ",
    "reference: ", REFERENCE_TREATMENT, " at ", REFERENCE_TIMEPOINT, " hr (", nrow(ref_tb), " rows)\n", sep = "")

nonfinite <- sum(vapply(markers, function(m) sum(!is.finite(tb[[m]]) & !is.na(tb[[m]])), numeric(1)))
if (nonfinite > 0) {
  note_warning(paste0(nonfinite, " non-finite value(s) (Inf/-Inf) in the tested columns. As in the published scripts, ",
                      "they are passed to the test as they are, which usually gives p = NA for that comparison."))
}



#### The test (as published) ###################################################
run_test <- function(x, y) {
  p <- tryCatch({
    if (length(x) < 2) {
      NA_real_
    } else if (TEST == "welch") {
      stats::t.test(x, y, paired = FALSE)$p.value
    } else if (TEST == "student") {
      stats::t.test(x, y, paired = FALSE, var.equal = TRUE)$p.value
    } else {
      suppressWarnings(stats::wilcox.test(x, y, paired = FALSE)$p.value)
    }
  }, error = function(e) NA_real_)
  as.numeric(p)
}
stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}

combos <- unique(tb[, .(TimePoint, Group)])
data.table::setorder(combos, TimePoint)
combos[, Group_Order := suppressWarnings(as.numeric(Group))]
if (anyNA(combos$Group_Order)) {
  # text groups (treatments): reference treatment first, then alphabetical
  combos[, Group_Order := match(Group, c(REFERENCE_TREATMENT, sort(setdiff(unique(Group), REFERENCE_TREATMENT))))]
}
data.table::setorder(combos, TimePoint, Group_Order)
grp_rows <- lapply(seq_len(nrow(combos)), function(i)
  which(abs(tb$TimePoint - combos$TimePoint[i]) < 1e-9 & tb$Group == combos$Group[i]))
use_matched <- COMPARE_BY != "Treatment" && REFERENCE_SCOPE == "matched"

res <- data.table::rbindlist(lapply(markers, function(m) {
  data.table::rbindlist(lapply(seq_len(nrow(combos)), function(i) {
    rows <- grp_rows[[i]]
    g    <- combos$Group[i]
    x    <- tb[[m]][rows]
    ref_idx <- if (use_matched) which(ref_tb$Group == g) else seq_len(nrow(ref_tb))
    y    <- ref_tb[[m]][ref_idx]
    data.table::data.table(
      Marker         = m,
      Compare_By     = COMPARE_BY,
      Group          = g,
      Treatment      = paste(unique(tb$Treatment[rows]), collapse = "+"),
      Cluster        = paste(unique(tb$Cluster[rows]), collapse = "+"),
      ThreshCat      = paste(unique(tb$ThreshCat[rows]), collapse = "+"),
      TimePoint      = combos$TimePoint[i],
      N_Group        = sum(!is.na(x)),
      N_Reference    = sum(!is.na(y)),
      Mean_Group     = if (any(!is.na(x))) mean(x, na.rm = TRUE) else NA_real_,
      Mean_Reference = if (any(!is.na(y))) mean(y, na.rm = TRUE) else NA_real_,
      Test           = TEST,
      p_value        = run_test(x, y))
  }))
}))
res[, Significance := stars(p_value)]
n_na <- sum(is.na(res$p_value))
if (n_na > 0) {
  note_warning(paste0(n_na, " of ", nrow(res), " comparisons have p = NA (a group or the reference has fewer than ",
                      "2 usable values, or the test could not be calculated)."))
}
small <- unique(res[N_Group < 3, .(Group, TimePoint, N_Group)])
if (nrow(small) > 0) {
  note_warning(paste0(nrow(small), " group/time point combination(s) have fewer than 3 samples (smallest n = ",
                      min(small$N_Group), "); with so few samples a p-value is very uncertain."))
}

# Published-layout wide table: one row per (time point, group), one column per marker
wide <- data.table::dcast(res, TimePoint + Group ~ Marker, value.var = "p_value")
wide <- wide[match(paste(combos$TimePoint, combos$Group), paste(wide$TimePoint, wide$Group))]
data.table::setcolorder(wide, c("TimePoint", "Group", markers))

# Litters behind the samples (Litter_Map.csv)
tb[, Key := paste(Cluster, ThreshCat, sep = "|")]
litter_file <- if (identical(LITTER_MAP, "none")) NA_character_ else if (nzchar(LITTER_MAP)) LITTER_MAP else
  file.path(dirname(HELPERS_DIR), "02_Analysis", "Litter_Map.csv")
has_litter <- FALSE
if (!is.na(litter_file)) {
  if (!file.exists(litter_file)) {
    note_warning(paste0("Litter map not found (", litter_file, "): the design table has no litter columns. ",
                        "Set LITTER_MAP to the file, or to \"none\" to silence this."))
  } else if (!"Replicate" %in% names(tb)) {
    note_warning("The table has no 'Replicate' column, so litters cannot be added to the design table.")
  } else {
    litter_map <- data.table::fread(litter_file, data.table = TRUE, colClasses = "character")
    if (!all(c("Sample_Condition", "Replicate", "Litter") %in% names(litter_map))) {
      stop("The litter map needs the columns Sample_Condition, Replicate and Litter.\nFile: ", litter_file, call. = FALSE)
    }
    if (anyDuplicated(litter_map[, .(Sample_Condition, Replicate)]) > 0) {
      stop("The litter map has the same Sample_Condition + Replicate more than once.\nFile: ", litter_file, call. = FALSE)
    }
    tb[, Litter := litter_map$Litter[match(paste(Sample_Condition, Replicate),
                                          paste(litter_map$Sample_Condition, litter_map$Replicate))]]
    n_unmapped <- data.table::uniqueN(tb[is.na(Litter), .(Sample_Condition, Replicate)])
    if (n_unmapped > 0) {
      note_warning(paste0(n_unmapped, " sample(s) are not in the litter map (litter shown as NA): ",
                          paste(utils::head(unique(tb[is.na(Litter), paste0(Sample_Condition, "_rep", Replicate)]), 5), collapse = ", "),
                          if (n_unmapped > 5) ", ..." else ""))
    }
    has_litter <- TRUE
  }
}

# Sample sizes behind the tests
lit_summary <- function(l) {
  l <- l[!is.na(l)]
  list(N_Litters = length(unique(l)), Litters = paste(sort(unique(l)), collapse = ","),
       Max_Samples_Per_Litter = if (length(l)) max(table(l)) else NA_integer_)
}
design <- tb[, c(list(N_Samples = .N, Replicates = paste(sort(unique(Replicate)), collapse = ",")),
                 if (has_litter) lit_summary(Litter)),
             by = .(Group, Treatment, TimePoint)]
data.table::setorder(design, TimePoint, Group)
litters_used <- if (has_litter) unique(tb[, .(Treatment, TimePoint, Sample_Condition, Replicate, Litter)]) else NULL
if (has_litter) {
  data.table::setorder(litters_used, Treatment, TimePoint, Replicate)
  shared <- design[!is.na(Max_Samples_Per_Litter) & Max_Samples_Per_Litter > 1]
  if (nrow(shared) > 0) {
    note_warning(paste0(nrow(shared), " group/time point combination(s) hold two or more samples from the SAME litter ",
                        "(e.g. ", shared$Treatment[1], " ", shared$TimePoint[1], " hr). The tests treat all samples as independent; ",
                        "see Statistics_design.csv."))
  }
}

# Assumption diagnostics (reported, never used to choose or skip a test)
skewness <- function(e) { m <- mean(e); d <- e - m; v <- mean(d^2); if (v == 0) NA_real_ else mean(d^3) / v^1.5 }
diagnose <- function(x, grp) {
  ok <- is.finite(x); x <- x[ok]; grp <- as.character(grp[ok])
  out <- list(N_Total = length(x), N_Groups = length(unique(grp)),
              Min_N = if (length(x)) min(table(grp)) else NA_integer_,
              Max_N = if (length(x)) max(table(grp)) else NA_integer_,
              Shapiro_W = NA_real_, Shapiro_p = NA_real_, Skewness = NA_real_,
              BrownForsythe_p = NA_real_, Variance_Ratio = NA_real_, Mean_SD_Spearman = NA_real_)
  if (length(x) < 3) return(out)
  n_by <- table(grp); keep <- names(n_by)[n_by >= 2]
  xs <- x[grp %in% keep]; gs <- grp[grp %in% keep]
  if (length(keep) >= 1 && length(xs) >= 3) {
    e <- xs - stats::ave(xs, gs)
    if (stats::sd(e) > 0) {
      sw <- tryCatch(stats::shapiro.test(e), error = function(err) NULL)
      if (!is.null(sw)) { out$Shapiro_W <- unname(sw$statistic); out$Shapiro_p <- sw$p.value }
      out$Skewness <- skewness(e)
    }
  }
  if (length(keep) >= 2) {
    v <- tapply(xs, gs, stats::var)
    if (max(v) > 0) out$Variance_Ratio <- if (min(v) > 0) max(v) / min(v) else Inf
    z <- abs(xs - stats::ave(xs, gs, FUN = stats::median))
    if (stats::var(z) > 0) {
      out$BrownForsythe_p <- tryCatch(stats::anova(stats::lm(z ~ factor(gs)))[["Pr(>F)"]][1],
                                      error = function(err) NA_real_)
    }
  }
  if (length(keep) >= 4) {
    mu <- tapply(xs, gs, mean); sdv <- tapply(xs, gs, stats::sd)
    if (stats::sd(mu) > 0 && stats::sd(sdv) > 0) {
      out$Mean_SD_Spearman <- suppressWarnings(stats::cor(mu, sdv, method = "spearman"))
    }
  }
  out
}

diag_tbl <- NULL
if (RUN_DIAGNOSTICS) {
  non_ref_tr <- setdiff(unique(tb$Treatment), REFERENCE_TREATMENT)
  keys <- unique(tb$Key)
  fam <- data.table::CJ(Treatment = non_ref_tr, Key = keys, sorted = FALSE)
  fam_rows <- lapply(seq_len(nrow(fam)), function(i) {
    r1 <- which(tb$Treatment == fam$Treatment[i] & tb$Key == fam$Key[i])
    r0 <- which(tb$Treatment == REFERENCE_TREATMENT &
                abs(tb$TimePoint - REFERENCE_TIMEPOINT) < 1e-9 & tb$Key == fam$Key[i])
    list(rows = c(r1, r0),
         grp  = c(if (length(r1)) paste0("t", tb$TimePoint[r1]) else character(0),   # [v2.1] paste0("t", numeric(0)) is "t", not character(0)
                  rep("ref", length(r0))))
  })
  fam_ok <- vapply(fam_rows, function(f) length(f$rows) >= 3, logical(1))
  want_log2 <- STAT_TABLE == "FC"
  diag_tbl <- data.table::rbindlist(lapply(markers, function(m) {
    data.table::rbindlist(lapply(which(fam_ok), function(i) {
      f <- fam_rows[[i]]; x <- tb[[m]][f$rows]
      d <- diagnose(x, f$grp)
      d2 <- if (want_log2) diagnose(suppressWarnings(log2(x)), f$grp) else NULL
      r1 <- f$rows[1]
      data.table::data.table(
        Marker = m, Treatment = fam$Treatment[i], Cluster = tb$Cluster[r1], ThreshCat = tb$ThreshCat[r1],
        data.table::as.data.table(d),
        Shapiro_p_log2        = if (want_log2) d2$Shapiro_p else NA_real_,
        BrownForsythe_p_log2  = if (want_log2) d2$BrownForsythe_p else NA_real_,
        Variance_Ratio_log2   = if (want_log2) d2$Variance_Ratio else NA_real_,
        N_Log2_Dropped        = if (want_log2) d$N_Total - d2$N_Total else NA_integer_)
    }))
  }))
  if (nrow(diag_tbl) == 0) {
    note_warning("No family had enough values for diagnostics (need a non-reference treatment plus the reference).")
    diag_tbl <- NULL
  } else {
    diag_tbl[, Flag_Normality := !is.na(Shapiro_p) & Shapiro_p < 0.05]
    diag_tbl[, Flag_EqualSpread := !is.na(BrownForsythe_p) & BrownForsythe_p < 0.05]
    diag_tbl[, Note := ifelse(!is.na(Min_N) & Min_N < 3, "smallest group n < 3: variance and normality checks have almost no power",
                       ifelse(!is.na(N_Total) & N_Total < 12, "fewer than 12 values: normality check has low power", ""))]
    n_fam <- nrow(diag_tbl)
    if (n_fam > 0) {
      cat("[INFO] Diagnostics: ", n_fam, " marker x family checks. Shapiro-Wilk p < 0.05: ",
          sum(diag_tbl$Flag_Normality), " (", sprintf("%.1f", 100 * mean(diag_tbl$Flag_Normality[!is.na(diag_tbl$Shapiro_p)])),
          "% of those testable; about 5% is expected by chance alone). Brown-Forsythe p < 0.05: ",
          sum(diag_tbl$Flag_EqualSpread), " (", sprintf("%.1f", 100 * mean(diag_tbl$Flag_EqualSpread[!is.na(diag_tbl$BrownForsythe_p)])),
          "%).\n", sep = "")
      if (STAT_TABLE == "FC") {
        cat("[INFO] On log2 scale: Shapiro p < 0.05 in ",
            sprintf("%.1f", 100 * mean(diag_tbl$Shapiro_p_log2 < 0.05, na.rm = TRUE)), "% (FC scale: ",
            sprintf("%.1f", 100 * mean(diag_tbl$Shapiro_p < 0.05, na.rm = TRUE)), "%); Brown-Forsythe p < 0.05 in ",
            sprintf("%.1f", 100 * mean(diag_tbl$BrownForsythe_p_log2 < 0.05, na.rm = TRUE)), "% (FC scale: ",
            sprintf("%.1f", 100 * mean(diag_tbl$BrownForsythe_p < 0.05, na.rm = TRUE)), "%).\n", sep = "")
      }
    }
  }
}



#### Save the tables ###########################################################
stage_label <- paste0("Statistics_", STAT_TABLE, "_", TEST, "_by", COMPARE_BY,
                      if (length(CLUSTERS) > 0) paste0("_Cl", paste(CLUSTERS, collapse = ",")) else "",
                      if (length(CATEGORIES) > 0) paste0("_Cat", paste(CATEGORIES, collapse = ",")) else "",
                      if (nzchar(OUTPUT_LABEL)) paste0("_", OUTPUT_LABEL) else "")
Stage_WD <- make_stage_dir(input_dir, stage_label)

out_files <- character(0)
f_long   <- file.path(Stage_WD, "Statistics_long.csv");        data.table::fwrite(res,    f_long)
f_wide   <- file.path(Stage_WD, "Statistics_pvalues_wide.csv"); data.table::fwrite(wide,   f_wide)
f_design <- file.path(Stage_WD, "Statistics_design.csv");      data.table::fwrite(design, f_design)
out_files <- c(f_long, f_wide, f_design)
if (has_litter) {
  f_lit <- file.path(Stage_WD, "Statistics_litters.csv"); data.table::fwrite(litters_used, f_lit)
  out_files <- c(out_files, f_lit)
}
if (!is.null(diag_tbl) && nrow(diag_tbl) > 0) {
  f_diag <- file.path(Stage_WD, "Statistics_diagnostics.csv"); data.table::fwrite(diag_tbl, f_diag)
  out_files <- c(out_files, f_diag)
}
cat("[INFO] Tables written to:\n       ", Stage_WD, "\n", sep = "")



#### Run record ################################################################
Inputs <- list(
  Script_Version = SCRIPT_VERSION,
  Completed_At   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER   = input_dir,
  OUTPUT_FOLDER  = Stage_WD,
  Settings = list(
    STAT_TABLE = STAT_TABLE, TEST = TEST, COMPARE_BY = COMPARE_BY,
    REFERENCE_TREATMENT = REFERENCE_TREATMENT, REFERENCE_TIMEPOINT = REFERENCE_TIMEPOINT,
    REFERENCE_SCOPE = REFERENCE_SCOPE,
    RUN_DIAGNOSTICS = RUN_DIAGNOSTICS, LITTER_MAP = LITTER_MAP,
    TREATMENTS = as.list(TREATMENTS), CLUSTERS = as.list(CLUSTERS), CATEGORIES = as.list(CATEGORIES),
    MARKERS = as.list(MARKERS),
    Adjustment = "none (as published)",
    Test_Detail = switch(TEST,
      welch   = "stats::t.test(x, y, paired = FALSE)  [Welch; the published method]",
      student = "stats::t.test(x, y, paired = FALSE, var.equal = TRUE)",
      wilcox  = "stats::wilcox.test(x, y, paired = FALSE)")),
  Input_Table            = basename(table_path),
  Input_Table_MD5        = unname(tools::md5sum(table_path)),
  Upstream_Module        = UPSTREAM_MODULE,
  Upstream_Script_Version = input_record[[UPSTREAM_MODULE]][["Script_Version"]],
  Table_Rows_Total       = n_rows_table,
  Table_Rows_Used        = nrow(tb),
  Reference_Rows         = nrow(ref_tb),
  Markers_Tested         = as.list(markers),
  Groups_Tested          = lapply(seq_len(nrow(combos)), function(i)
    list(TimePoint = combos$TimePoint[i], Group = combos$Group[i], N = length(grp_rows[[i]]))),
  Comparisons            = nrow(res),
  Comparisons_NA         = n_na,
  Litter_Map_File        = if (has_litter) litter_file else "not used",
  Litter_Map_MD5         = if (has_litter) unname(tools::md5sum(litter_file)) else "not used",
  Diagnostics            = if (!is.null(diag_tbl) && nrow(diag_tbl) > 0) list(
    Families_Checked = nrow(diag_tbl),
    Shapiro_p_below_0.05 = sum(diag_tbl$Flag_Normality),
    BrownForsythe_p_below_0.05 = sum(diag_tbl$Flag_EqualSpread)) else "not run",
  Random_Seed            = "not used in this script",
  Warnings               = as.list(warnings_log),
  Output_MD5             = file_md5(out_files),
  Session                = session_info_list(c("data.table", "jsonlite")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = input_dir)

cat("\n========================================================\n")
cat("[DONE] ", nrow(res), " comparisons (", length(markers), " markers x ", nrow(combos), " groups/time points)\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("       Next: check against the published tables with Tools/Compare_GoldenFiles.R (run P1)\n")
cat("========================================================\n")
