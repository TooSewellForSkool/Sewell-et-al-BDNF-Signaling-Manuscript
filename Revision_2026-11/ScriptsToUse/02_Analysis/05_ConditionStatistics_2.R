################################################################################
# 05_ConditionStatistics.R
################################################################################
# Layer:    Analysis -- statistics for the revision (no plots; the Visualization layer reads the tables)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-26
# Edited:   2026-10-26
# Version:  2   (NEW; not part of the published analysis -- written for the eLife revision, RP 111513)
#
# Purpose:
# * Compare the CONDITIONS of one experiment with each other at one time point (the bar graphs:
#   Cont. Depr., BDNF test, BDNF + K252a test), for every marker, within every cluster / category.
# * Answer the reviewer's statistics comment: an omnibus test (one-way ANOVA or Welch's ANOVA) across the
#   three conditions, then ALL THREE pairwise comparisons with a correction for the number of comparisons,
#   and a second correction (Benjamini-Hochberg) reported next to it.
# * The published bar graphs showed unadjusted Welch t-test p-values for the same three pairs
#   (checked: 41 of 41 published panels reproduce). Those numbers stay in the output as p_unadjusted.
#
# Where it fits:   01_PercAboveThreshold.R  (or the merged-category table)  ->  05_ConditionStatistics.R
#   INPUT_FOLDER is a folder that holds the fold-change table (file name contains "_FC", not "RepsAvg").
#   One row per sample x cluster/category, one column per marker (as written by 01_PercAboveThreshold.R).
#
# ==============================================================================
# WHAT IS COMPARED  ("family" = one marker in one cluster/category)
# ==============================================================================
#   Groups: the conditions listed in CONDITIONS (default: ContDepr, BDNFTest, BDNF+K252aTest).
#   Pairs, in this fixed order:  1 vs 2,  1 vs 3,  2 vs 3   (ContDepr-BDNFTest, ContDepr-BDNF+K252aTest,
#   BDNFTest-BDNF+K252aTest; the same order as the published bars, bottom to top).
#   Samples are treated as independent. Litters are not modelled.
#
# ==============================================================================
# THE P-VALUE COLUMNS (nothing is dropped; read the column names)
# ==============================================================================
#   p_unadjusted      pairwise Welch t-test (R default, unpaired), no adjustment = the published p-values
#   p_Holm            Holm adjustment of p_unadjusted over the 3 pairs of the family
#   p_primary         PRIMARY result, adjusted over the 3 pairs of the family. With the default settings
#                     (POSTHOC "welch_holm") it equals p_Holm, for EVERY marker (decision 2026-10-26).
#   p_Tukey           Tukey HSD from the classic one-way ANOVA (assumes equal spread), always calculated when
#                     3 groups have 2+ samples, for the supplementary tables. Omnibus_p_classic = classic ANOVA p.
#   q_BH_all_markers  ALTERNATIVE: Benjamini-Hochberg over p_unadjusted of ALL markers x 3 pairs of the
#                     same cluster/category. Its size depends on how many markers were analysed.
#   Omnibus_p / q_omnibus_BH_markers   omnibus p and its BH-adjusted version across the markers.
#
# ==============================================================================
# THE SUGGESTION (decision tree; it never blocks anything)
# ==============================================================================
#   D1  Independent groups, one factor, all pairs  -> one-way ANOVA + all-pairs post-hoc (Tukey).
#   D2  Equal spread?  If the smallest group has fewer than 3 samples, OR the largest/smallest variance ratio
#       is above VAR_RATIO_LIMIT, OR Brown-Forsythe p < 0.05 -> "welch" (Welch's ANOVA + Welch t-tests with
#       Holm). Otherwise -> "classic" (classic ANOVA + Tukey).
#   D3  Normality (Shapiro-Wilk on residuals) is REPORTED but does not change the test (n = 2-4 per group).
#   Override with OMNIBUS and POSTHOC below ("auto" follows the suggestion).
#
# ==============================================================================
# CHANGE LOG: new script (no earlier version)
# ==============================================================================
#  A1 [RESULT, NEW] Different tests from the published ones on purpose (omnibus test, adjustment). The
#     published p-values are kept as p_unadjusted.
#  A2 [BEHAVIOR] A group with fewer than 2 usable values is left out of the models and its pairs get NA.
#     With fewer than 3 usable groups there is no omnibus test and no Tukey (the pair falls back to Holm).
#  A4 [RESULT, v2] Default changed from the spread-based suggestion to ONE method for every marker:
#     Welch's ANOVA + Welch t-tests with Holm (OMNIBUS "welch", POSTHOC "welch_holm"). The suggestion
#     (Suggested_Mode) is still written to Conditions_families.csv; Tukey and the classic ANOVA are reported
#     in the columns p_Tukey and Omnibus_p_classic. Set both settings to "auto" to get the v1 behaviour.
#  A3 [BEHAVIOR] Whether stars are shown at low n is decided by the plotting script (MIN_N_FOR_STARS);
#     here every p-value that can be calculated is reported, with the group sizes next to it.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
INPUT_FOLDER   <- ""            # folder holding the fold-change table ("" = choose in a pop-up)
TABLE_FILE     <- ""            # "" = the single "_FC" table in the folder (asked if there are several)

CONDITIONS     <- c("ContDepr", "BDNFTest", "BDNF+K252aTest")  # order sets the pairs (see above)
CLUSTERS       <- character(0)  # e.g. c("1","2","3"); empty = all clusters/categories in the table
MARKERS        <- character(0)  # e.g. c("pAKt_168Er","pCreb_176Lu"); empty = every marker column

OMNIBUS        <- "welch"       # "welch" (Welch's ANOVA; default), "classic" (ANOVA) or "auto" (follow the suggestion)
POSTHOC        <- "welch_holm"  # "welch_holm" (default), "tukey" (needs equal variances) or "auto"
VAR_RATIO_LIMIT <- 3            # D2: largest/smallest group variance above this -> suggest "welch"
OUTPUT_LABEL   <- ""            # optional text added to the output folder name

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "05_ConditionStatistics.R"
SCRIPT_VERSION           <- 2
MODULE_BASE_NAME         <- "05_ConditionStatistics"
REQUIRED_HELPERS_VERSION <- 7
ID_COLUMNS <- c("File_Num", "File_Name", "Sample_Condition", "Reordered", "TimePoint", "TimePoints",
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
cat("[START] Compare conditions (ANOVA + all pairs) for every marker and cluster/category:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

warnings_log <- character(0)
note_warning <- function(txt) {
  warnings_log <<- c(warnings_log, txt)
  cat("[WARNING] ", txt, "\n", sep = "")
}
one_of <- function(value, allowed, name) {
  if (!is.character(value) || length(value) != 1 || !value %in% allowed) {
    stop(name, " must be one of: ", paste0("\"", allowed, "\"", collapse = ", "), ".", call. = FALSE)
  }
}
one_of(OMNIBUS, c("auto", "classic", "welch"), "OMNIBUS")
one_of(POSTHOC, c("auto", "tukey", "welch_holm"), "POSTHOC")
if (!is.numeric(VAR_RATIO_LIMIT) || length(VAR_RATIO_LIMIT) != 1 || VAR_RATIO_LIMIT <= 1) {
  stop("VAR_RATIO_LIMIT must be one number above 1 (e.g. 3).", call. = FALSE)
}
if (!is.character(CONDITIONS) || length(CONDITIONS) < 2 || anyDuplicated(CONDITIONS)) {
  stop("CONDITIONS must hold at least two different condition names.", call. = FALSE)
}
if (OMNIBUS != "auto" && POSTHOC != "auto" && ((OMNIBUS == "classic") != (POSTHOC == "tukey"))) {
  note_warning(paste0("OMNIBUS = \"", OMNIBUS, "\" with POSTHOC = \"", POSTHOC, "\" mixes an equal-variance ",
                      "and an unequal-variance method. That is allowed, but check that it is what you want."))
}



#### Choose the input folder and the table ######################################
fc_tables_in <- function(d) {
  f <- list.files(d, pattern = "\\.csv$", full.names = FALSE)
  f[grepl("_FC", f) & !grepl("RepsAvg|log2FC", f)]
}
if (nzchar(INPUT_FOLDER)) {
  if (!dir.exists(INPUT_FOLDER)) stop("INPUT_FOLDER not found:\n  ", INPUT_FOLDER, call. = FALSE)
  input_dir <- normalizePath(INPUT_FOLDER, winslash = "/", mustWork = TRUE)
} else {
  Original_WD <- get_working_dir("project")
  dirs <- list.dirs(Original_WD, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[!grepl("/[_.][^/]*", sub(Original_WD, "", dirs, fixed = TRUE))]
  dirs <- sort(dirs[vapply(dirs, function(d) length(fc_tables_in(d)) > 0, logical(1))], decreasing = TRUE)
  if (length(fc_tables_in(Original_WD)) > 0) dirs <- c(Original_WD, dirs)
  dirs <- unique(dirs)
  if (length(dirs) == 0) stop("No folder with a fold-change table (*_FC.csv) was found in:\n  ", Original_WD,
                              "\nRun 01_PercAboveThreshold.R first, or set INPUT_FOLDER.", call. = FALSE)
  labels <- basename(dirs)
  pick <- ask_choice("Which folder holds the table to analyse?", labels, preselect = labels[1])
  input_dir <- normalizePath(dirs[match(pick, labels)], winslash = "/", mustWork = TRUE)
}
cat("[INFO] Input folder:\n       ", input_dir, "\n", sep = "")
tabs <- fc_tables_in(input_dir)
if (nzchar(TABLE_FILE)) {
  if (!file.exists(file.path(input_dir, TABLE_FILE))) stop("TABLE_FILE not found in the input folder: ", TABLE_FILE, call. = FALSE)
  table_name <- TABLE_FILE
} else if (length(tabs) == 1) {
  table_name <- tabs
} else if (length(tabs) > 1) {
  table_name <- ask_choice("Which table should be analysed?", tabs, preselect = tabs[1])
} else stop("No *_FC.csv table in the input folder.", call. = FALSE)
table_path <- file.path(input_dir, table_name)
prior_record <- read_run_record(input_dir)
if (is.null(prior_record)) {
  note_warning("No run_parameters.json in the input folder: the table is treated as data made before run records existed. Provenance cannot be checked.")
}
cat("[INFO] Table: ", table_name, "\n", sep = "")



#### Read the table #############################################################
tb <- data.table::fread(table_path, na.strings = c("", "NA"), data.table = TRUE)
if (names(tb)[1] %in% c("", "V1") && ncol(tb) > 1 &&
    identical(as.integer(tb[[1]]), seq_len(nrow(tb)))) tb[, (1) := NULL]
if (!"Sample_Condition" %in% names(tb)) stop("The table has no column 'Sample_Condition'.", call. = FALSE)
if (!"Cluster" %in% names(tb)) tb[, Cluster := NA_character_]
if (!"ThreshCat" %in% names(tb)) tb[, ThreshCat := ""]
tb[, ThreshCat := as.character(ThreshCat)]
tb[is.na(ThreshCat), ThreshCat := ""]
tb[, Cluster := as.character(Cluster)]
tb[, Condition := sub("_[^_]+$", "", Sample_Condition)]
tb <- tb[Condition %in% CONDITIONS]
missing_cond <- setdiff(CONDITIONS, unique(tb$Condition))
if (length(missing_cond) > 0) {
  stop("These CONDITIONS were not found in the table: ", paste(missing_cond, collapse = ", "),
       "\nConditions present (Sample_Condition without its last '_' part): ",
       paste(unique(sub("_[^_]+$", "", data.table::fread(table_path, select = "Sample_Condition")[[1]])), collapse = ", "),
       call. = FALSE)
}
if (length(CLUSTERS) > 0) tb <- tb[Cluster %in% as.character(CLUSTERS)]
if (nrow(tb) == 0) stop("No rows left after selecting CONDITIONS and CLUSTERS.", call. = FALSE)
tb[, Key := paste(Cluster, ThreshCat, sep = "|")]
marker_cols <- setdiff(names(tb), c(ID_COLUMNS, "Condition", "Key"))
marker_cols <- marker_cols[vapply(marker_cols, function(m) is.numeric(tb[[m]]), logical(1))]
markers <- if (length(MARKERS) == 0) marker_cols else MARKERS
if (length(setdiff(markers, marker_cols)) > 0) {
  stop("MARKERS not found in the table: ", paste(setdiff(markers, marker_cols), collapse = ", "), call. = FALSE)
}
cat("[INFO] ", nrow(tb), " table rows; ", length(markers), " markers; ", data.table::uniqueN(tb$Key),
    " cluster/category groups; conditions: ", paste(CONDITIONS, collapse = ", "), "\n", sep = "")



#### The tests ##################################################################
stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}
fmt <- function(v, digits = 3) ifelse(is.na(v), "NA", formatC(v, digits = digits, format = "g"))
pair_idx <- utils::combn(length(CONDITIONS), 2)       # columns = pairs: (1,2),(1,3),(2,3) for three groups

# Checks for one family (x = values, g = factor of conditions)
diagnose <- function(x, g) {
  n_by <- as.integer(table(g)); keep <- n_by >= 2
  v <- tapply(x, g, stats::var)[keep]
  vr <- if (length(v) >= 2 && all(is.finite(v))) {
    if (min(v) > 0) max(v) / min(v) else Inf
  } else NA_real_
  bf <- tryCatch({
    xs <- x[g %in% levels(g)[keep]]; gs <- droplevels(g[g %in% levels(g)[keep]])
    med <- tapply(xs, gs, stats::median)
    summary(stats::aov(abs(xs - med[as.character(gs)]) ~ gs))[[1]][["Pr(>F)"]][1]
  }, error = function(e) NA_real_)
  sw <- tryCatch({
    xs <- x[g %in% levels(g)[keep]]; gs <- droplevels(g[g %in% levels(g)[keep]])
    r <- stats::resid(stats::lm(xs ~ gs))
    if (length(r) >= 3 && length(unique(r)) > 1) stats::shapiro.test(r)$p.value else NA_real_
  }, error = function(e) NA_real_)
  list(Min_N = min(n_by), Variance_Ratio = vr, BrownForsythe_p = bf, Shapiro_p = sw)
}

suggest_mode <- function(d) {
  reasons <- character(0)
  if (is.na(d$Min_N) || d$Min_N < 3) reasons <- c(reasons, paste0("smallest group n = ", d$Min_N, " (spread cannot be judged)"))
  if (!is.na(d$Variance_Ratio) && d$Variance_Ratio > VAR_RATIO_LIMIT) {
    reasons <- c(reasons, paste0("variance ratio ", fmt(d$Variance_Ratio), " > ", VAR_RATIO_LIMIT))
  }
  if (!is.na(d$BrownForsythe_p) && d$BrownForsythe_p < 0.05) {
    reasons <- c(reasons, paste0("Brown-Forsythe p = ", fmt(d$BrownForsythe_p)))
  }
  mode <- if (length(reasons) > 0) "welch" else "classic"
  if (mode == "classic") reasons <- "spread looks similar (ratio and Brown-Forsythe within limits, n >= 3)"
  if (!is.na(d$Shapiro_p) && d$Shapiro_p < 0.05) {
    reasons <- c(reasons, paste0("note: Shapiro-Wilk p = ", fmt(d$Shapiro_p), " (not used to change the test)"))
  }
  list(mode = mode, reason = paste(reasons, collapse = "; "))
}

# One family. x = values; g = factor over CONDITIONS
analyse_family <- function(x, g, omnibus, posthoc) {
  ok <- is.finite(x); x <- x[ok]; g <- droplevels_keep(g[ok])
  n_by <- as.integer(table(g)); names(n_by) <- levels(g)
  mean_by <- tapply(x, g, mean); mean_by[n_by == 0] <- NA
  res <- data.frame(Group1 = CONDITIONS[pair_idx[1, ]], Group2 = CONDITIONS[pair_idx[2, ]],
                    N_Group1 = as.integer(unname(n_by[CONDITIONS[pair_idx[1, ]]])), N_Group2 = as.integer(unname(n_by[CONDITIONS[pair_idx[2, ]]])),
                    Mean_Group1 = as.numeric(mean_by[CONDITIONS[pair_idx[1, ]]]),
                    Mean_Group2 = as.numeric(mean_by[CONDITIONS[pair_idx[2, ]]]),
                    p_unadjusted = NA_real_, p_Holm = NA_real_, p_primary = NA_real_, p_Tukey = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(res))) {
    if (res$N_Group1[i] >= 2 && res$N_Group2[i] >= 2) {
      res$p_unadjusted[i] <- tryCatch(stats::t.test(x[g == res$Group1[i]], x[g == res$Group2[i]], paired = FALSE)$p.value,
                                      error = function(e) NA_real_)
    }
  }
  res$p_Holm <- stats::p.adjust(res$p_unadjusted, method = "holm")
  res$p_primary <- res$p_Holm
  used <- CONDITIONS[n_by[CONDITIONS] >= 2]
  omni_p <- NA_real_; omni_classic <- NA_real_
  if (length(used) >= 3) {
    sel <- g %in% used
    d <- data.frame(x = x[sel], g = factor(as.character(g[sel]), levels = used))
    omni_p <- tryCatch({
      if (omnibus == "classic") summary(stats::aov(x ~ g, data = d))[[1]][["Pr(>F)"]][1]
      else stats::oneway.test(x ~ g, data = d, var.equal = FALSE)$p.value
    }, error = function(e) NA_real_)
    omni_classic <- tryCatch(summary(stats::aov(x ~ g, data = d))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
    tk <- tryCatch(stats::TukeyHSD(stats::aov(x ~ g, data = d))$g, error = function(e) NULL)
    if (!is.null(tk)) {
      for (i in seq_len(nrow(res))) {
        nm <- paste0(res$Group2[i], "-", res$Group1[i])
        if (nm %in% rownames(tk)) res$p_Tukey[i] <- tk[nm, "p adj"]
      }
    }
    if (posthoc == "tukey") res$p_primary <- res$p_Tukey
  }
  list(omnibus_p = as.numeric(omni_p), omnibus_classic = as.numeric(omni_classic), pairs = res, n_by = n_by)
}
droplevels_keep <- function(g) factor(as.character(g), levels = CONDITIONS)

long_list <- list(); fam_list <- list(); gs_list <- list(); k <- 0L
keys <- unique(tb$Key)
for (key in keys) {
  rows <- which(tb$Key == key)
  clus <- tb$Cluster[rows[1]]; thr <- tb$ThreshCat[rows[1]]
  g <- factor(tb$Condition[rows], levels = CONDITIONS)
  for (m in markers) {
    x <- tb[[m]][rows]
    ok <- is.finite(x)
    if (sum(ok) < 4) next
    dg <- diagnose(x[ok], droplevels_keep(g[ok]))
    sg <- suggest_mode(dg)
    omni_used <- if (OMNIBUS == "auto") sg$mode else OMNIBUS
    post_used <- if (POSTHOC == "auto") (if (sg$mode == "classic") "tukey" else "welch_holm") else POSTHOC
    a <- analyse_family(x, g, omni_used, post_used)
    n_used <- sum(a$n_by >= 2)
    k <- k + 1L
    pr <- a$pairs
    if (n_used < 3) post_shown <- "welch_holm (fewer than 3 usable groups)" else post_shown <- post_used
    long_list[[k]] <- data.table::data.table(
      Marker = m, Cluster = clus, ThreshCat = thr, Pair = seq_len(nrow(pr)), pr[, c("Group1", "Group2")],
      N_Group1 = pr$N_Group1, N_Group2 = pr$N_Group2, Mean_Group1 = pr$Mean_Group1, Mean_Group2 = pr$Mean_Group2,
      Omnibus_Test = if (n_used >= 3) omni_used else "none (fewer than 3 groups)", Omnibus_p = a$omnibus_p,
      Omnibus_p_classic = a$omnibus_classic, PostHoc = post_shown, p_unadjusted = pr$p_unadjusted, p_Holm = pr$p_Holm,
      p_primary = pr$p_primary, p_Tukey = pr$p_Tukey)
    fam_list[[k]] <- data.table::data.table(
      Marker = m, Cluster = clus, ThreshCat = thr, N_Total = sum(ok), Min_N = dg$Min_N,
      Variance_Ratio = dg$Variance_Ratio, BrownForsythe_p = dg$BrownForsythe_p, Shapiro_p = dg$Shapiro_p,
      Suggested_Mode = sg$mode, Suggestion_Reason = sg$reason,
      Omnibus_Test_Used = if (n_used >= 3) omni_used else "none (fewer than 3 groups)",
      PostHoc_Used = post_shown, Omnibus_p = a$omnibus_p)
    xs <- split(x[ok], droplevels_keep(g[ok]))
    gs_list[[k]] <- data.table::data.table(
      Marker = m, Cluster = clus, ThreshCat = thr, Condition = CONDITIONS,
      N = vapply(CONDITIONS, function(cn) length(xs[[cn]]), integer(1)),
      Mean = vapply(CONDITIONS, function(cn) if (length(xs[[cn]])) mean(xs[[cn]]) else NA_real_, numeric(1)),
      SD = vapply(CONDITIONS, function(cn) if (length(xs[[cn]]) >= 2) stats::sd(xs[[cn]]) else NA_real_, numeric(1)))
  }
}
if (length(long_list) == 0) stop("No marker/cluster combination had enough values (at least 4 in total) to test.", call. = FALSE)
long <- data.table::rbindlist(long_list)
fams <- data.table::rbindlist(fam_list)
gstat <- data.table::rbindlist(gs_list)

# Alternative family: BH over all markers x pairs of a cluster/category, on p_unadjusted
long[, q_BH_all_markers := stats::p.adjust(p_unadjusted, method = "BH"), by = .(Cluster, ThreshCat)]
long[, Family_Size_BH := sum(!is.na(p_unadjusted)), by = .(Cluster, ThreshCat)]
fams[, q_omnibus_BH_markers := stats::p.adjust(Omnibus_p, method = "BH"), by = .(Cluster, ThreshCat)]
long[, Significance := stars(p_primary)]
long[, Scale := "FC"]
data.table::setorder(long, Cluster, ThreshCat, Marker, Pair)
data.table::setorder(fams, Cluster, ThreshCat, Marker)
n_na <- sum(is.na(long$p_primary))
if (n_na > 0) {
  note_warning(paste0(n_na, " of ", nrow(long), " comparisons have p = NA (a group has fewer than 2 values)."))
}
if (any(long$N_Group1 < 3 | long$N_Group2 < 3, na.rm = TRUE)) {
  note_warning(paste0(sum(long$N_Group1 < 3 | long$N_Group2 < 3, na.rm = TRUE),
                      " comparison(s) involve a group with fewer than 3 samples; interpret with care."))
}
n_sig_unadj   <- sum(long$p_unadjusted < 0.05, na.rm = TRUE)
n_sig_primary <- sum(long$p_primary < 0.05, na.rm = TRUE)
n_sig_bh      <- sum(long$q_BH_all_markers < 0.05, na.rm = TRUE)
cat("[INFO] Comparisons with p (or q) < 0.05:  unadjusted ", n_sig_unadj, ";  primary (3 pairs per marker) ",
    n_sig_primary, ";  BH over all markers ", n_sig_bh, "  (of ", sum(!is.na(long$p_unadjusted)), ")\n", sep = "")
cat("[INFO] Suggested mode: ", paste(names(table(fams$Suggested_Mode)), table(fams$Suggested_Mode), collapse = ", "),
    " (families = marker x cluster/category)\n", sep = "")



#### Save the tables ###########################################################
stage_label <- paste0("Conditions_", OMNIBUS, "-", POSTHOC, if (nzchar(OUTPUT_LABEL)) paste0("_", OUTPUT_LABEL) else "")
Stage_WD <- make_stage_dir(input_dir, stage_label)
f1 <- file.path(Stage_WD, "Conditions_long.csv");       data.table::fwrite(long,  f1)
f2 <- file.path(Stage_WD, "Conditions_families.csv");   data.table::fwrite(fams,  f2)
f3 <- file.path(Stage_WD, "Conditions_groupstats.csv"); data.table::fwrite(gstat, f3)
f4 <- file.path(Stage_WD, "Conditions_input_values.csv")
data.table::fwrite(tb[, c("Sample_Condition", "Condition", "Cluster", "ThreshCat",
                          intersect(c("Replicate", "TimePoint", "TimePoints"), names(tb)), markers), with = FALSE], f4)
out_files <- c(f1, f2, f3, f4)
cat("[INFO] Tables written to:\n       ", Stage_WD, "\n", sep = "")



#### Run record ################################################################
Inputs <- list(
  Script_Version = SCRIPT_VERSION,
  Completed_At   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER   = input_dir,
  OUTPUT_FOLDER  = Stage_WD,
  Settings = list(
    CONDITIONS = as.list(CONDITIONS), CLUSTERS = as.list(as.character(CLUSTERS)), MARKERS = as.list(MARKERS),
    OMNIBUS = OMNIBUS, POSTHOC = POSTHOC, VAR_RATIO_LIMIT = VAR_RATIO_LIMIT,
    Pairs = as.list(paste(CONDITIONS[pair_idx[1, ]], "vs", CONDITIONS[pair_idx[2, ]])),
    Primary_Family = "the 3 pairs of one marker in one cluster/category (default: Welch t-tests + Holm for every marker; Tukey in p_Tukey)",
    Alternative_Family = "Benjamini-Hochberg over all markers x pairs of one cluster/category, on p_unadjusted"),
  Input_Table            = table_name,
  Input_Table_MD5        = unname(tools::md5sum(table_path)),
  Input_Provenance       = if (is.null(prior_record)) "no run record (legacy data)" else "run record present",
  Families_Analysed      = nrow(fams),
  Comparisons            = nrow(long),
  Comparisons_NA         = n_na,
  Significant_Tukey      = sum(long$p_Tukey < 0.05, na.rm = TRUE),
  Significant_Unadjusted = n_sig_unadj,
  Significant_Primary    = n_sig_primary,
  Significant_BH_All_Markers = n_sig_bh,
  Suggested_Mode_Counts  = as.list(table(fams$Suggested_Mode)),
  Random_Seed            = "not used (no random numbers in this script)",
  Warnings               = as.list(warnings_log),
  Output_MD5             = file_md5(out_files),
  Session                = session_info_list(c("data.table", "jsonlite")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = input_dir)

cat("\n========================================================\n")
cat("[DONE] ", nrow(long), " comparisons in ", nrow(fams), " marker x cluster/category families\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("       Main result: Conditions_long.csv (column p_primary); suggestions: Conditions_families.csv\n")
cat("========================================================\n")
