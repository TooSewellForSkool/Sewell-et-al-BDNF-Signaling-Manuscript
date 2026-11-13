################################################################################
# 04_TimecourseStatistics.R
################################################################################
# Layer:    Analysis -- statistics for the revision (no plots; the Visualization layer reads the tables)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-25
# Edited:   2026-10-28
# Version:  1.2   (NEW; not part of the published analysis -- written for the eLife revision, RP 111513)
#
# Purpose:
# * Answer the reviewer's statistics comment: compare the time points of ONE treatment for ONE marker with
#   a one-way ANOVA (omnibus test), then compare every time point with the reference (NoTreat 0 hr) with a
#   post-hoc test that accounts for the number of comparisons, and report a second adjustment next to it.
# * Suggest the test from the design and the assumption checks, and let you override the suggestion.
#
# Where it fits:   01_PercAboveThreshold.R  ->  03_Statistics.R (published tests + design + diagnostics)
#                  ->  04_TimecourseStatistics.R (this script).
#   INPUT_FOLDER is the folder that 03_Statistics.R (version 2 or newer, RUN_DIAGNOSTICS = TRUE) wrote.
#   This script re-reads the SAME table with the SAME selection (treatments, clusters, markers) that 03
#   recorded, so the two scripts always describe the same data.
#
# ==============================================================================
# WHAT IS COMPARED  ("family" = one treatment, one cluster/category)
# ==============================================================================
#   Within a family, the groups are the reference (NoTreat at 0 hr) and the time points of the treatment.
#   Each marker is analysed on its own.  For each marker and family:
#     1. Omnibus test: is there any difference among the groups?   (classic ANOVA or Welch's ANOVA)
#     2. Post-hoc: every time point against the reference.           (Dunnett, or Welch t-tests + Holm)
#   The reference is the same cluster's/category's NoTreat 0 hr (never pooled over clusters).
#
# ==============================================================================
# THE P-VALUE COLUMNS (nothing is dropped; read the column names)
# ==============================================================================
#   p_unadjusted      pairwise Welch t-test time point vs reference, no adjustment
#                     (the same calculation as the published tables, for comparison)
#   p_primary         PRIMARY result. Adjusted within the family "this treatment x this marker"
#                     (Dunnett's adjustment, or Holm applied to p_unadjusted).
#   q_BH_all_markers  ALTERNATIVE family, reported alongside: Benjamini-Hochberg false discovery rate over
#                     p_unadjusted of ALL markers x time points of the same treatment (and cluster/category).
#                     Use it to ask "how many of the many marker results could be chance?".
#   Omnibus_p / q_omnibus_BH_markers   omnibus p-value and its BH-adjusted version across markers.
#
# ==============================================================================
# THE SUGGESTION (a decision tree; it never blocks anything)
# ==============================================================================
#   Design first, thenf checks that 03_Statistics.R wrote (Statistics_diagnostics.csv):
#   D1  Independent groups, one factor (time), every time point against ONE reference
#       -> one-way ANOVA + many-to-one post-hoc (Dunnett). Samples are treated as independent
#          (see Statistics_litters.csv for the litters; litters are not modelled here).
#   D2  Equal spread across groups?  If the smallest group has fewer than 3 samples (so spread cannot be
#       judged), OR the largest/smallest variance ratio is above VAR_RATIO_LIMIT, OR Brown-Forsythe p < 0.05
#       -> "welch" (Welch's ANOVA + pairwise Welch t-tests with Holm). Otherwise -> "classic"
#       (classic ANOVA + Dunnett, which assume equal variances).
#   D3  Normality (Shapiro-Wilk on residuals) is REPORTED in the reason text but does not change the test:
#       with 2-4 samples per group a rank-based test cannot reach small p-values and normality cannot be
#       verified either way.
#   D4  Scale: with STAT_TABLE "FC" the log2 diagnostics are compared with the FC diagnostics
#       (Suggested_Scale). It is a suggestion only; re-run 03 and 04 with STAT_TABLE "log2FC" to use it.
#   Override with OMNIBUS and POSTHOC below ("auto" follows the suggestion).
#
# ==============================================================================
# CHANGE LOG: new script (no earlier version)
# ==============================================================================
#  A1 [RESULT, NEW] Different tests from the published ones on purpose (ANOVA, adjustment). The published
#     tests stay available in 03_Statistics.R; nothing published is changed.
#  A2 [RESULT, declared] Dunnett p-values come from a multivariate-t integration that uses random numbers
#     (package multcomp). The seed is set (SEED below) and recorded, so a re-run gives the same numbers; a
#     different seed changes the 4th-5th decimal place.
#  A5 [BUGFIX, v1.1, no change to any result that ran before] A marker whose values in a group of a family were ALL
#     missing or infinite (typical for small expression categories, e.g. TrkB-negative / p75NTR-high) stopped the
#     script with "arguments imply differing number of rows: 0, 1". Such a group now stays in the tables with n = 0
#     and NA p-values. Output for data that ran before is identical.
#  A6 [BUGFIX, v1.2, no change to any result that ran before] A treatment with NO rows in a cluster/category (e.g. BDNF in a
#     small expression category) was counted as having one time point, because paste0("t", numeric(0)) returns "t".
#     The group labels were then one longer than the values and the script stopped with "arguments must have same
#     length". Such families are now skipped, with the warning that was already written for them.
#  A3 [BEHAVIOR] A group with fewer than 2 usable values is left out of the models and gets NA p-values.
#     With fewer than 3 groups there is no omnibus test (NA).
#  A4 [DECISION NEEDED LATER] Litters are not modelled. A mixed model with litter as a random effect is a
#     possible long-term upgrade; with 3 litters per time point it would be weakly estimated.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
INPUT_FOLDER   <- ""            # a folder made by 03_Statistics.R (version 2+) ("" = choose in a pop-up)

OMNIBUS        <- "welch"        # "auto" (follow the suggestion), "classic" (ANOVA) or "welch" (Welch's ANOVA)
POSTHOC        <- "welch_holm"        # "auto", "dunnett" (needs equal variances) or "welch_holm"
VAR_RATIO_LIMIT <- 3            # D2: largest/smallest group variance above this -> suggest "welch"
SEED           <- 1             # for Dunnett's random-number integration (recorded in the run record)
OUTPUT_LABEL   <- ""            # optional text added to the output folder name

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "04_TimecourseStatistics.R"
SCRIPT_VERSION           <- "1.2"
MODULE_BASE_NAME         <- "04_TimecourseStatistics"
REQUIRED_HELPERS_VERSION <- 7
UPSTREAM_MODULE          <- "03_Statistics"
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
need_multcomp <- POSTHOC %in% c("auto", "dunnett")
pkg_versions <- dc_load_packages(cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite",
                                          if (need_multcomp) "multcomp"))

cat("\n========================================================\n")
cat("[START] ANOVA across time points and post-hoc tests against the reference:\n",
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
one_of(POSTHOC, c("auto", "dunnett", "welch_holm"), "POSTHOC")
if (!is.numeric(VAR_RATIO_LIMIT) || length(VAR_RATIO_LIMIT) != 1 || VAR_RATIO_LIMIT <= 1) {
  stop("VAR_RATIO_LIMIT must be one number above 1 (e.g. 3).", call. = FALSE)
}
if (!is.numeric(SEED) || length(SEED) != 1) stop("SEED must be one number.", call. = FALSE)
if (OMNIBUS != "auto" && POSTHOC != "auto" &&
    ((OMNIBUS == "classic") != (POSTHOC == "dunnett"))) {
  note_warning(paste0("OMNIBUS = \"", OMNIBUS, "\" with POSTHOC = \"", POSTHOC, "\" mixes an equal-variance ",
                      "and an unequal-variance method. That is allowed, but check that it is what you want."))
}



#### Choose the input folder (a 03_Statistics.R output folder) #################
is_stat_folder <- function(d) file.exists(file.path(d, "Statistics_diagnostics.csv")) &&
                              file.exists(file.path(d, RUN_RECORD_NAME))
if (nzchar(INPUT_FOLDER)) {
  if (!dir.exists(INPUT_FOLDER)) stop("INPUT_FOLDER not found:\n  ", INPUT_FOLDER, call. = FALSE)
  input_dir <- normalizePath(INPUT_FOLDER, winslash = "/", mustWork = TRUE)
} else {
  Original_WD <- get_working_dir("project")
  dirs <- list.dirs(Original_WD, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[!grepl("/[_.][^/]*", sub(Original_WD, "", dirs, fixed = TRUE))]
  dirs <- sort(dirs[vapply(dirs, is_stat_folder, logical(1))], decreasing = TRUE)
  if (is_stat_folder(Original_WD)) dirs <- c(Original_WD, dirs)
  dirs <- unique(dirs)
  if (length(dirs) == 0) {
    stop("No 03_Statistics.R output folder (with Statistics_diagnostics.csv) was found in:\n  ", Original_WD,
         "\nRun 03_Statistics.R (version 2 or newer, RUN_DIAGNOSTICS = TRUE) first, or set INPUT_FOLDER.", call. = FALSE)
  }
  labels <- basename(dirs)
  pick <- ask_choice("Which Statistics folder should be analysed?", labels, preselect = labels[1])
  input_dir <- normalizePath(dirs[match(pick, labels)], winslash = "/", mustWork = TRUE)
}
cat("[INFO] Input folder:\n       ", input_dir, "\n", sep = "")

record <- require_modules(input_dir, UPSTREAM_MODULE, how_to_fix = "Run 03_Statistics.R first.")
s3 <- record[[UPSTREAM_MODULE]]
if (is.null(s3$Script_Version) || as.numeric(s3$Script_Version) < 2) {
  stop("The Statistics folder was made by 03_Statistics.R version ", s3$Script_Version,
       ". Version 2 or newer (with the diagnostics) is needed.", call. = FALSE)
}
f_diag <- file.path(input_dir, "Statistics_diagnostics.csv")
if (!file.exists(f_diag)) {
  stop("Statistics_diagnostics.csv is missing. Run 03_Statistics.R with RUN_DIAGNOSTICS = TRUE.", call. = FALSE)
}

# The table behind the Statistics folder: it sits in the folder the Statistics folder was made in
as_chr <- function(v) if (is.null(v) || length(v) == 0) character(0) else as.character(unlist(v))
table_path <- file.path(dirname(input_dir), s3$Input_Table)
if (!file.exists(table_path)) {
  stop("The table that 03_Statistics.R used was not found:\n  ", table_path,
       "\nThe Statistics folder must stay inside the folder of the table (made by 01_PercAboveThreshold.R).", call. = FALSE)
}
if (!identical(unname(as.character(tools::md5sum(table_path))), as.character(s3$Input_Table_MD5))) {
  stop("The table has changed since 03_Statistics.R used it (checksum differs):\n  ", table_path,
       "\nRun 03_Statistics.R again.", call. = FALSE)
}
STAT_TABLE          <- s3$Settings$STAT_TABLE
REFERENCE_TREATMENT <- s3$Settings$REFERENCE_TREATMENT
REFERENCE_TIMEPOINT <- as.numeric(s3$Settings$REFERENCE_TIMEPOINT)
TREATMENTS <- as_chr(s3$Settings$TREATMENTS)
CLUSTERS   <- as_chr(s3$Settings$CLUSTERS)
CATEGORIES <- as_chr(s3$Settings$CATEGORIES)
MARKERS    <- as_chr(s3$Settings$MARKERS)
cat("[INFO] Table: ", basename(table_path), "  (", STAT_TABLE, ")\n", sep = "")



#### Read the table and apply the same selection as 03_Statistics.R ############
tb <- data.table::fread(table_path, na.strings = c("", "NA"), data.table = TRUE)
if (names(tb)[1] %in% c("", "V1") && ncol(tb) > 1 &&
    identical(as.integer(tb[[1]]), seq_len(nrow(tb)))) tb[, (1) := NULL]
for (cn in c("Sample_Condition", "TimePoint", "Cluster")) {
  if (!cn %in% names(tb)) stop("The table has no column '", cn, "'.", call. = FALSE)
}
if (!"ThreshCat" %in% names(tb)) tb[, ThreshCat := ""]
tb[, ThreshCat := as.character(ThreshCat)]
tb[is.na(ThreshCat), ThreshCat := ""]
tb[, Treatment := sub("_[^_]+$", "", Sample_Condition)]
if (length(TREATMENTS) > 0) tb <- tb[Treatment %in% TREATMENTS]
if (length(CLUSTERS)   > 0) tb <- tb[as.character(Cluster) %in% CLUSTERS]
if (length(CATEGORIES) > 0) tb <- tb[ThreshCat %in% CATEGORIES]
tb[, Key := paste(Cluster, ThreshCat, sep = "|")]
marker_cols <- setdiff(names(tb), c(ID_COLUMNS, "Treatment", "Key"))
markers <- if (length(MARKERS) == 0) marker_cols else MARKERS
if (length(setdiff(markers, marker_cols)) > 0) {
  stop("MARKERS not found in the table: ", paste(setdiff(markers, marker_cols), collapse = ", "), call. = FALSE)
}
is_ref_row <- tb$Treatment == REFERENCE_TREATMENT & abs(tb$TimePoint - REFERENCE_TIMEPOINT) < 1e-9
treatments_tested <- setdiff(unique(tb$Treatment), REFERENCE_TREATMENT)
if (length(treatments_tested) == 0) {
  stop("The selection holds only the reference treatment (", REFERENCE_TREATMENT, "); there is nothing to compare.", call. = FALSE)
}
cat("[INFO] ", nrow(tb), " table rows; ", length(markers), " markers; treatments tested: ",
    paste(treatments_tested, collapse = ", "), "; reference: ", REFERENCE_TREATMENT, " at ",
    REFERENCE_TIMEPOINT, " hr\n", sep = "")

# Checks written by 03_Statistics.R (one row per marker x family)
dg <- data.table::fread(f_diag, na.strings = c("", "NA"), data.table = TRUE)
dg[, ThreshCat := ifelse(is.na(ThreshCat), "", as.character(ThreshCat))]
dg[, Key := paste(Cluster, ThreshCat, sep = "|")]
dg_lookup <- stats::setNames(seq_len(nrow(dg)), paste(dg$Marker, dg$Treatment, dg$Key))
if (!all(c("Variance_Ratio", "BrownForsythe_p", "Shapiro_p", "Min_N") %in% names(dg))) {
  stop("Statistics_diagnostics.csv has no variance ratio / Brown-Forsythe / Shapiro columns.", call. = FALSE)
}



#### The tests ##################################################################
stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}
fmt <- function(v, digits = 3) ifelse(is.na(v), "NA", formatC(v, digits = digits, format = "g"))

# The suggestion for one family (decision tree D2/D3; see the header)
suggest_mode <- function(d) {
  reasons <- character(0)
  if (is.na(d$Min_N) || d$Min_N < 3) {
    reasons <- c(reasons, paste0("smallest group n = ", d$Min_N, " (spread cannot be judged)"))
  }
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

# One family: x = values, g = group labels ("ref" or "t<hours>")
analyse_family <- function(x, g, omnibus, posthoc) {
  lv <- unique(g)                                   # [v1.1] keep groups whose values are all missing (n = 0)
  ok <- is.finite(x); x <- x[ok]; g <- factor(g[ok], levels = lv)
  n_by <- table(g)
  out_groups <- setdiff(lv, "ref")
  out_groups <- out_groups[order(as.numeric(sub("^t", "", out_groups)))]
  res <- data.frame(Group = out_groups,
                    N_Group = as.integer(n_by[out_groups]),
                    Mean_Group = as.numeric(tapply(x, g, mean)[out_groups]),
                    N_Reference = if ("ref" %in% names(n_by)) as.integer(n_by[["ref"]]) else 0L,
                    Mean_Reference = if ("ref" %in% names(n_by)) mean(x[g == "ref"]) else NA_real_,
                    p_unadjusted = NA_real_, p_primary = NA_real_, stringsAsFactors = FALSE)
  omni_p <- NA_real_
  if (length(out_groups) == 0) return(list(omnibus_p = omni_p, groups = res))
  if (res$N_Reference[1] < 2) return(list(omnibus_p = omni_p, groups = res))
  can_test <- res$N_Group >= 2
  for (i in which(can_test)) {
    res$p_unadjusted[i] <- tryCatch(stats::t.test(x[g == res$Group[i]], x[g == "ref"], paired = FALSE)$p.value,
                                    error = function(e) NA_real_)
  }
  keep <- c("ref", res$Group[can_test])
  if (length(keep) < 2) return(list(omnibus_p = omni_p, groups = res))
  sel <- g %in% keep
  d <- data.frame(x = x[sel], g = factor(g[sel], levels = keep))
  if (length(keep) >= 3) {
    omni_p <- tryCatch({
      if (omnibus == "classic") summary(stats::aov(x ~ g, data = d))[[1]][["Pr(>F)"]][1]
      else stats::oneway.test(x ~ g, data = d, var.equal = FALSE)$p.value
    }, error = function(e) NA_real_)
  }
  if (posthoc == "dunnett") {
    pv <- tryCatch({
      fit <- stats::aov(x ~ g, data = d)
      mc  <- multcomp::glht(fit, linfct = multcomp::mcp(g = "Dunnett"))
      sm  <- summary(mc)
      p   <- as.numeric(sm$test$pvalues)
      names(p) <- sub(" - ref$", "", names(sm$test$coefficients))
      p
    }, error = function(e) NULL)
    if (!is.null(pv)) res$p_primary[match(names(pv), res$Group)] <- as.numeric(pv)
  } else {
    res$p_primary[can_test] <- stats::p.adjust(res$p_unadjusted[can_test], method = "holm")
  }
  list(omnibus_p = as.numeric(omni_p), groups = res)
}

set.seed(SEED)
fam <- data.table::CJ(Treatment = treatments_tested, Key = unique(tb$Key), sorted = FALSE)
fam_info <- lapply(seq_len(nrow(fam)), function(i) {
  r1 <- which(tb$Treatment == fam$Treatment[i] & tb$Key == fam$Key[i])
  r0 <- which(is_ref_row & tb$Key == fam$Key[i])
  list(rows = c(r1, r0),
       grp = c(if (length(r1)) paste0("t", tb$TimePoint[r1]) else character(0), rep("ref", length(r0))),   # [v1.2] paste0("t", numeric(0)) is "t", not character(0)
       cluster = if (length(c(r1, r0))) tb$Cluster[c(r1, r0)[1]] else NA,
       thresh  = if (length(c(r1, r0))) tb$ThreshCat[c(r1, r0)[1]] else "")
})
fam_ok <- vapply(fam_info, function(f) sum(f$grp != "ref") >= 1 && sum(f$grp == "ref") >= 2, logical(1))
if (!any(fam_ok)) stop("No treatment/cluster combination has both time points and a reference with 2+ samples.", call. = FALSE)
if (any(!fam_ok)) {
  note_warning(paste0(sum(!fam_ok), " treatment/cluster combination(s) were skipped: no time point values or fewer than 2 reference samples."))
}

long_list <- list(); fam_list <- list(); k <- 0L; n_missing_diag <- 0L
for (m in markers) {
  for (i in which(fam_ok)) {
    f <- fam_info[[i]]
    key_name <- paste(m, fam$Treatment[i], fam$Key[i])
    d_row <- if (key_name %in% names(dg_lookup)) as.list(dg[dg_lookup[[key_name]]]) else NULL
    if (is.null(d_row)) {
      n_missing_diag <- n_missing_diag + 1L
      d_row <- list(Min_N = NA_integer_, Variance_Ratio = NA_real_, BrownForsythe_p = NA_real_,
                    Shapiro_p = NA_real_, Variance_Ratio_log2 = NA_real_, N_Total = NA_integer_)
    }
    sg <- suggest_mode(d_row)
    omni_used <- if (OMNIBUS == "auto") sg$mode else OMNIBUS
    post_used <- if (POSTHOC == "auto") (if (sg$mode == "classic") "dunnett" else "welch_holm") else POSTHOC
    a <- analyse_family(tb[[m]][f$rows], f$grp, omni_used, post_used)
    n_groups_used <- sum(a$groups$N_Group >= 2) + 1L
    k <- k + 1L
    gdf <- a$groups
    long_list[[k]] <- data.table::data.table(
      Marker = m, Treatment = fam$Treatment[i], Cluster = f$cluster, ThreshCat = f$thresh,
      TimePoint = as.numeric(sub("^t", "", gdf$Group)),
      N_Group = gdf$N_Group, N_Reference = gdf$N_Reference,
      Mean_Group = gdf$Mean_Group, Mean_Reference = gdf$Mean_Reference,
      Omnibus_Test = if (n_groups_used >= 3) omni_used else "none (fewer than 3 groups)",
      Omnibus_p = a$omnibus_p, PostHoc = post_used,
      p_unadjusted = gdf$p_unadjusted, p_primary = gdf$p_primary)
    fam_list[[k]] <- data.table::data.table(
      Marker = m, Treatment = fam$Treatment[i], Cluster = f$cluster, ThreshCat = f$thresh,
      N_Total = d_row$N_Total, Min_N = d_row$Min_N,
      Variance_Ratio = d_row$Variance_Ratio, BrownForsythe_p = d_row$BrownForsythe_p, Shapiro_p = d_row$Shapiro_p,
      Suggested_Mode = sg$mode, Suggestion_Reason = sg$reason,
      Suggested_Scale = if (STAT_TABLE == "FC" && !is.null(d_row$Variance_Ratio_log2) &&
                            !is.na(d_row$Variance_Ratio_log2) && !is.na(d_row$Variance_Ratio) &&
                            d_row$Variance_Ratio_log2 < d_row$Variance_Ratio) "log2FC" else STAT_TABLE,
      Omnibus_Test_Used = if (n_groups_used >= 3) omni_used else "none (fewer than 3 groups)",
      PostHoc_Used = post_used, Omnibus_p = a$omnibus_p)
  }
}
if (n_missing_diag > 0) {
  note_warning(paste0(n_missing_diag, " marker/family combination(s) have no row in Statistics_diagnostics.csv; ",
                      "for those the suggestion falls back to \"welch\"."))
}
long <- data.table::rbindlist(long_list)
fams <- data.table::rbindlist(fam_list)

# Alternative family: Benjamini-Hochberg over all markers x time points of a treatment (and cluster/category)
long[, q_BH_all_markers := stats::p.adjust(p_unadjusted, method = "BH"), by = .(Treatment, Cluster, ThreshCat)]
fams[, q_omnibus_BH_markers := stats::p.adjust(Omnibus_p, method = "BH"), by = .(Treatment, Cluster, ThreshCat)]
long[, Significance := stars(p_primary)]
long[, Scale := STAT_TABLE]
data.table::setorder(long, Treatment, Cluster, ThreshCat, TimePoint, Marker)
n_na <- sum(is.na(long$p_primary))
if (n_na > 0) {
  note_warning(paste0(n_na, " of ", nrow(long), " post-hoc comparisons have p = NA (a group has fewer than 2 values, ",
                      "or the test could not be calculated)."))
}
if (any(long$Omnibus_Test == "none (fewer than 3 groups)")) {
  note_warning("Some families have fewer than 3 usable groups: no omnibus test there (NA).")
}

# Wide layout for the plotting scripts: rows = treatment x cluster/category x time point, columns = markers
wide <- data.table::dcast(long, Treatment + Cluster + ThreshCat + TimePoint ~ Marker, value.var = "p_primary")
data.table::setcolorder(wide, c("Treatment", "Cluster", "ThreshCat", "TimePoint", markers))

n_sig_unadj   <- sum(long$p_unadjusted < 0.05, na.rm = TRUE)
n_sig_primary <- sum(long$p_primary < 0.05, na.rm = TRUE)
n_sig_bh      <- sum(long$q_BH_all_markers < 0.05, na.rm = TRUE)
cat("[INFO] Comparisons with p (or q) < 0.05:  unadjusted ", n_sig_unadj, ";  primary (within treatment x marker) ",
    n_sig_primary, ";  BH over all markers ", n_sig_bh, "  (of ", sum(!is.na(long$p_unadjusted)), ")\n", sep = "")
cat("[INFO] Suggested mode: ", paste(names(table(fams$Suggested_Mode)), table(fams$Suggested_Mode), collapse = ", "),
    " (families = marker x treatment x cluster/category)\n", sep = "")



#### Save the tables ###########################################################
stage_label <- paste0("Timecourse_", STAT_TABLE, "_", OMNIBUS, "-", POSTHOC,
                      if (nzchar(OUTPUT_LABEL)) paste0("_", OUTPUT_LABEL) else "")
Stage_WD <- make_stage_dir(input_dir, stage_label)
f1 <- file.path(Stage_WD, "Timecourse_long.csv");           data.table::fwrite(long,  f1)
f2 <- file.path(Stage_WD, "Timecourse_families.csv");       data.table::fwrite(fams,  f2)
f3 <- file.path(Stage_WD, "Timecourse_pvalues_wide.csv");   data.table::fwrite(wide,  f3)
out_files <- c(f1, f2, f3)
cat("[INFO] Tables written to:\n       ", Stage_WD, "\n", sep = "")



#### Run record ################################################################
Inputs <- list(
  Script_Version = SCRIPT_VERSION,
  Completed_At   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER   = input_dir,
  OUTPUT_FOLDER  = Stage_WD,
  Settings = list(
    OMNIBUS = OMNIBUS, POSTHOC = POSTHOC, VAR_RATIO_LIMIT = VAR_RATIO_LIMIT, SEED = SEED,
    STAT_TABLE = STAT_TABLE, REFERENCE_TREATMENT = REFERENCE_TREATMENT, REFERENCE_TIMEPOINT = REFERENCE_TIMEPOINT,
    TREATMENTS = as.list(TREATMENTS), CLUSTERS = as.list(CLUSTERS), CATEGORIES = as.list(CATEGORIES),
    MARKERS = as.list(MARKERS),
    Primary_Family = "one treatment x one marker (x one cluster/category): time points against the reference",
    Alternative_Family = "Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted"),
  Input_Table            = basename(table_path),
  Input_Table_MD5        = unname(tools::md5sum(table_path)),
  Statistics_Folder      = basename(input_dir),
  Statistics_Script_Version = s3$Script_Version,
  Families_Analysed      = nrow(fams),
  Comparisons            = nrow(long),
  Comparisons_NA         = n_na,
  Significant_Unadjusted = n_sig_unadj,
  Significant_Primary    = n_sig_primary,
  Significant_BH_All_Markers = n_sig_bh,
  Suggested_Mode_Counts  = as.list(table(fams$Suggested_Mode)),
  Random_Seed            = SEED,
  Warnings               = as.list(warnings_log),
  Output_MD5             = file_md5(out_files),
  Session                = session_info_list(c("data.table", "jsonlite", if (need_multcomp) "multcomp")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = input_dir)

cat("\n========================================================\n")
cat("[DONE] ", nrow(long), " comparisons in ", nrow(fams), " marker x treatment families\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("       Main result: Timecourse_long.csv (column p_primary); suggestions: Timecourse_families.csv\n")
cat("========================================================\n")
