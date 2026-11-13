################################################################################
# Compare_GoldenFiles.R -- check that a NEW script gives the same results as the OLD one
################################################################################
# Author:   Jonathon Sewell (with Claude)
# Version:  3.5   (2.0 = pairs, composite keys; 2.1 = per-pair folders, declared dropped/extra columns;
#                  3.0 = ONE published-figures root folder + checksum check of every OLD file,
#                  settings check of every NEW run folder, so a wrong folder is caught, not guessed;
#                  3.1 = 13 named runs: A-C (tables), D95/D75/D50 (12.1.1 percentiles),
#                  T1-T4 (threshold categories), K1-K3 (analysis within categories);
#                  (v3.4: P3 tolerance 5e-4)
#                   3.2 = run P1 (03_Statistics.R vs published Supplementary Table 3), per-pair tolerance,
#                  settings checks can look at an earlier module's section of the run record;
#                  3.3 = rows only in NEW allowed for values not in OLD (P1: more markers than published),
#                  run E (01 v3, 20 clusters, 12.1.1 Clusters tables), runs P2-P4 (Supplementary Table 4);
#                  3.5 = Script_Version in the settings check means 'this version or newer' (was: the text must contain the old
#                  version number, so a rebuild with 01 v4 was rejected although only the settings matter))
# Created:  2026-10-21
# Edited:   2026-10-29 (v3.5)
#
# PURPOSE
#   Run the OLD (published) result and the NEW script on the same data, then use
#   this tool to show that the output files match. Differences are listed in plain
#   language and marked as EXPECTED (you declared them) or UNEXPECTED. The tool
#   also writes a one-line-per-comparison record (CSV) with file checksums, which
#   you can keep with the project or send for review.
#
# CHANGES FROM v1 (the settings block changed; v1 is in _superseded)
#   * Several comparisons in one run; each is one "pair" (see PAIRS below)
#   * Rows are matched on one OR SEVERAL key columns (e.g. File_Num + Cluster)
#   * Columns can be renamed on the OLD side (published tables used other names)
#   * Columns can be ignored, declared as expected to differ, or declared as new
#   * Reports the number of values compared and the largest difference
#   * Writes golden_record_<time>.csv (pair, files, MD5s, rows, values compared,
#     largest difference, counts of unexpected/expected differences, PASS/FAIL)
#   * Negative control (default on): changes one value in a COPY of the first pair's
#     NEW table and confirms the tool notices. If it does not, the check is not
#     trustworthy and the run is marked FAIL.
#
# HOW TO USE
#   1. Example A (the golden checks): fill in GOLDEN_ROOT and, in NEW_RUNS, the folder of each run
#      you have made. Leave a run's folder "" to skip that run.
#      Other uses (Example B): fill in OLD_FOLDER and NEW_FOLDER, then edit PAIRS.
#   2. Click "Source". Read the report, then keep (or send) the two output files.
#
# WHAT IS COMPARED
#   Values are read as text, then compared as numbers when both sides are numbers
#   (so "8" equals 8) and as text otherwise. "NaN" is treated as NA (missing) if
#   NAN_AS_NA is TRUE. A first column with an empty heading and values 1..n (the
#   row-number column written by write.csv) is dropped automatically.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
# --- Golden check for 02_Analysis/01_PercAboveThreshold.R (Example A) ---------------
# GOLDEN_ROOT: the published-figures folder that was scanned for the file manifest, i.e. the
#   folder named "2026-02-23 Figures_v12,13" (it contains "Figure 2 - PseudoBulk" etc.).
#   Every published (OLD) file below is looked up inside it by its path in the manifest and
#   checked by MD5 checksum, so the wrong file or folder is reported, never compared.
GOLDEN_ROOT <- "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-02-23 Figures_v12,13"
# GOLDEN_STATS_FOLDER: the folder "GoldenData" that comes with these scripts (Tools/GoldenData). It holds the
#   published p-values (made from Supplementary_Tables.xlsx) that run P1 is compared with. Its files are
#   checked by MD5 checksum too.
GOLDEN_STATS_FOLDER <- "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/ScriptsToUse/Tools/GoldenData"
# NEW_RUNS: one folder per run you have made (the folder the script created inside the input folder,
#   named <date_time>_<something>). Leave "" for runs you have not made. README_Analysis.txt lists
#   exactly which script and settings each run needs; the tool checks them and says if a folder is wrong.
NEW_RUNS <- list(
  A   = "",   # 01_PercAboveThreshold:  PseudoBulk, no exclusions                         (45 samples)
  B   = "",   # 01_PercAboveThreshold:  PseudoBulk, CONDITIONS_TO_EXCLUDE c("CompMed")    (32 samples)
  C   = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-43-28_PercAboveThresh_Clusters_1234_CLmerged",   # 01_PercAboveThreshold:  Clusters, CLUSTERS_MERGED TRUE                    (45 samples x clusters 1-4)
  D95 = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-39-20_PercAboveThresh_PseudoBulk_1234_PB",   # 01 (v2): as A, LOG2FC_INF_AS_NA TRUE, THRESHOLD_PROB 0.95  -> 12.1.1 Fig 3D tables
  D75 = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-39-46_PercAboveThresh_PseudoBulk_1234_PB75",   # 01 (v2): as D95 with THRESHOLD_PROB 0.75                    -> SuppFig 3
  D50 = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-40-32_PercAboveThresh_PseudoBulk_1234_PB50",   # 01 (v2): as D95 with THRESHOLD_PROB 0.5                     -> SuppFig 3
  T1  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-48-58_ThreshCat_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7",   # 02_AssignThresholdCategories: TrkB_i_174Yb, cuts 0.15, 1.4, 3.7
  T2  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-49-24_ThreshCat_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95",   # 02_AssignThresholdCategories: p75NTR_161Dy, cuts 0.45, 1.45, 3.95
  T3  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-50-12_ThreshCat_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7",   # 02_AssignThresholdCategories: p75NTR_161Dy 0.45, 1.45, 3.95 + TrkB_i_174Yb 0.15, 1.4, 3.7
  T4  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-50-33_ThreshCat_p75NTR_161Dy-Gr1.45_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7",   # 02_AssignThresholdCategories: p75NTR_161Dy 1.45 + TrkB_i_174Yb 0.15, 1.4, 3.7
  K1  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-52-00_PercAboveThresh_PseudoBulk_1234_Cat1-2-3-4_T1_TrkB",   # 01 (v2): categories from T1, all categories, PseudoBulk           -> Fig 5 TrkB_i
  K2  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-53-30_PercAboveThresh_PseudoBulk_1234_Cat11-41-14-44_T3_6AB",   # 01 (v2): categories from T3, CATEGORIES_TO_ANALYSE 11,41,14,44    -> Fig 6A
  K3  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-52-54_PercAboveThresh_PseudoBulk_1234_Cat1-2-3-4_T2_p75",   # 01 (v2): categories from T2, all categories, PseudoBulk           -> p75NTR (Fig 5 folder 02-27)
  P1  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-39-20_PercAboveThresh_PseudoBulk_1234_PB/2026-10-29_16-57-16_Statistics_RAWvalues_welch_byTreatment_Fig2C",   # 03_Statistics (v1): input = the run A folder, RAWvalues, Welch, by Treatment -> Supplementary Table 3
  E   = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_16-34-51_PercAboveThresh_Clusters_1234_CL20",   # 01 (v3): Clusters, CLUSTER_ROUND 1, CLUSTERS_MERGED FALSE (20 clusters), LOG2FC_INF_AS_NA TRUE,
              #          FILL_EMPTY_CLUSTERS TRUE                                          -> 12.1.1 Fig 3 cluster tables (900 rows)
  P2  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_16-34-51_PercAboveThresh_Clusters_1234_CL20/2026-10-29_17-31-09_Statistics_FC_welch_byCluster_Cl1-2-4-7-13_Table_S4A",   # 03_Statistics: input = the run E folder, FC, Welch, by Cluster, NoTreat + BDNF, clusters 1,2,4,7,13 -> Table S4A (Fig 3C)
  P3  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_16-34-51_PercAboveThresh_Clusters_1234_CL20/2026-10-29_17-31-58_Statistics_FC_welch_byCluster_Cl1-3-15-17_Table_S4B",   # 03_Statistics: input = the run E folder, FC, Welch, by Cluster, NoTreat + CompMed, clusters 1,3,15,17 -> Table S4B (Fig 3D)
  P4  = "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_16-34-51_PercAboveThresh_Clusters_1234_CL20/2026-10-29_17-32-33_Statistics_FC_welch_byCluster_Cl9-17-18_Table_S4C"    # 03_Statistics: input = the run E folder, FC, Welch, by Cluster, NoTreat + BDNF, clusters 9,17,18 -> Table S4C (Fig 3E)
)

# --- Other uses (Example B) ----------------------------------------------------------
OLD_FOLDER <- ""     # folder with the OLD (published) files ("" = ask)
NEW_FOLDER <- ""     # folder with the NEW script's files ("" = ask)
REPORT_FOLDER <- ""  # where to save the report and record ("" = the first NEW run folder, or NEW_FOLDER)

NUMERIC_TOLERANCE <- 1e-9   # numbers closer than this are equal (rounding only); a pair can set its own (see pair())
NAN_AS_NA         <- TRUE   # TRUE: "NaN" and missing are the same
MAX_EXAMPLES      <- 5      # example differences shown per column
NEGATIVE_CONTROL  <- TRUE   # prove that the tool notices a change (see above)

# Each comparison is one pair(...). Arguments:
#   label       short name shown in the report
#   old, new    file names (relative to OLD_FOLDER / NEW_FOLDER) or full paths
#   key         column(s) that identify a row, e.g. "File_Name" or c("File_Num", "Cluster")
#               (NA = compare row by row, in order)
#   rename_old  OLD columns to rename before comparing: c(NewName = "OldName")
#   ignore      columns not compared (listed in the report)
#   expected    columns whose differences you EXPECT (listed, but not a failure)
#   new_columns columns that only exist in NEW on purpose (listed, not a failure)
#   extra_new_ok TRUE = any column that only exists in NEW is allowed (listed, not a failure)
#   dropped_columns columns that only exist in OLD on purpose (listed, not a failure)
#   old_folder, new_folder  folders for this pair ("" = OLD_FOLDER / NEW_FOLDER above)
#   require_folders TRUE = a blank folder means "skip this pair" instead of asking for a folder
#   tolerance   this pair's own numeric tolerance (NA = NUMERIC_TOLERANCE)
#   extra_new_rows_by  column name: rows that exist only in NEW are allowed when their value in this column does
#               not appear anywhere in OLD (e.g. "Marker": the new table tests more markers than were published;
#               every published marker must still match completely). "" = no extra rows allowed.
#   old_root    "" = OLD file is inside GOLDEN_ROOT; "stats" = inside GOLDEN_STATS_FOLDER
#   old_rel, old_md5, run   golden mode: OLD file path inside GOLDEN_ROOT, its expected MD5, and which
#               NEW_RUN_x folder holds the NEW file
pair <- function(label, old = "", new, key = NA, rename_old = character(0), ignore = character(0),
                 expected = character(0), new_columns = character(0), extra_new_ok = FALSE,
                 dropped_columns = character(0), old_folder = "", new_folder = "", require_folders = FALSE,
                 old_rel = "", old_md5 = "", run = "", tolerance = NA, old_root = "", extra_new_rows_by = "") {
  if (!nzchar(old) && nzchar(old_rel)) old <- basename(old_rel)
  list(label = label, old = old, new = new, key = key, rename_old = rename_old,
       ignore = ignore, expected = expected, new_columns = new_columns, extra_new_ok = extra_new_ok,
       dropped_columns = dropped_columns, old_folder = old_folder, new_folder = new_folder,
       require_folders = require_folders, old_rel = old_rel, old_md5 = old_md5, run = run,
       tolerance = tolerance, old_root = old_root, extra_new_rows_by = extra_new_rows_by)
}

# --- Example A: golden definitions (the paths and MD5s come from manifest_2026-10-22_17-20-24.csv;
#     nothing to edit below this line until Example B) -------------------------------------------
# What each run must have used. eq = exact value, has / hasnt = the text appears / does not appear.
# (checked against the run_parameters.json inside the run folder; FOLDER_NAME = the folder's own name)
# module = (optional) check another module's section of the run record, e.g. the earlier script's settings
eq    <- function(key, value, module = NULL)  list(key = key, op = "eq",    value = value, module = module)
has   <- function(key, value, module = NULL)  list(key = key, op = "has",   value = value, module = module)
hasnt <- function(key, value, module = NULL)  list(key = key, op = "hasnt", value = value, module = module)
atleast <- function(key, value, module = NULL) list(key = key, op = "atleast", value = value, module = module)  # number >= value
M01 <- "01_PercAboveThreshold"; M02 <- "02_AssignThresholdCategories"; M03 <- "03_Statistics"
CAT_TRKB <- "ThreshCat_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7"
CAT_P75  <- "ThreshCat_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95"
CAT_BOTH <- "ThreshCat_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7"
CAT_P75_1CUT <- "ThreshCat_p75NTR_161Dy-Gr1.45_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7"
d_common <- function(p) list(eq("ANALYSIS_BY", "PseudoBulk"), eq("CONDITIONS_TO_EXCLUDE", ""),
                             eq("THRESHOLD_PROB", p), eq("LOG2FC_INF_AS_NA", "true"),
                             eq("THRESHOLD_CATEGORY_FOLDER", ""), atleast("Script_Version", "2"))
# what runs P2-P4 (03_Statistics on the 20-cluster FC table of run E) must have used
p_cluster_spec <- function(treatments, clusters) list(
  eq("STAT_TABLE", "FC"), eq("TEST", "welch"), eq("COMPARE_BY", "Cluster"), eq("REFERENCE_TREATMENT", "NoTreat"),
  eq("REFERENCE_TIMEPOINT", "0"), eq("REFERENCE_SCOPE", "pooled"), eq("TREATMENTS", treatments),
  eq("CLUSTERS", clusters), eq("CATEGORIES", ""),
  eq("ANALYSIS_BY", "Clusters", module = M01), eq("CLUSTER_ROUND", "1", module = M01), eq("CLUSTERS_MERGED", "false", module = M01),
  eq("FILL_EMPTY_CLUSTERS", "true", module = M01), atleast("Script_Version", "3", module = M01))
RUN_SPEC <- list(
  A = list(module = M01, expect = list(eq("ANALYSIS_BY", "PseudoBulk"), eq("CONDITIONS_TO_EXCLUDE", ""))),
  B = list(module = M01, expect = list(eq("ANALYSIS_BY", "PseudoBulk"), eq("CONDITIONS_TO_EXCLUDE", "CompMed"))),
  C = list(module = M01, expect = list(eq("ANALYSIS_BY", "Clusters"), eq("CONDITIONS_TO_EXCLUDE", ""), eq("CLUSTERS_MERGED", "true"))),
  D95 = list(module = M01, expect = d_common("0.95")),
  D75 = list(module = M01, expect = d_common("0.75")),
  D50 = list(module = M01, expect = d_common("0.5")),
  T1 = list(module = M02, expect = list(has("FOLDER_NAME", CAT_TRKB), eq("CUTS_BY_QUANTILE", "false"))),
  T2 = list(module = M02, expect = list(has("FOLDER_NAME", CAT_P75), hasnt("FOLDER_NAME", "TrkB"), eq("CUTS_BY_QUANTILE", "false"))),
  T3 = list(module = M02, expect = list(has("FOLDER_NAME", CAT_BOTH), eq("CUTS_BY_QUANTILE", "false"))),
  T4 = list(module = M02, expect = list(has("FOLDER_NAME", CAT_P75_1CUT), eq("CUTS_BY_QUANTILE", "false"))),
  K1 = list(module = M01, expect = list(eq("ANALYSIS_BY", "PseudoBulk"), eq("CONDITIONS_TO_EXCLUDE", ""), eq("LOG2FC_INF_AS_NA", "true"),
                                        eq("THRESHOLD_PROB", "0.95"), has("THRESHOLD_CATEGORY_FOLDER", CAT_TRKB), eq("CATEGORIES_TO_ANALYSE", ""))),
  K2 = list(module = M01, expect = list(eq("ANALYSIS_BY", "PseudoBulk"), eq("CONDITIONS_TO_EXCLUDE", ""), eq("LOG2FC_INF_AS_NA", "true"),
                                        eq("THRESHOLD_PROB", "0.95"), has("THRESHOLD_CATEGORY_FOLDER", CAT_BOTH), eq("CATEGORIES_TO_ANALYSE", "11,41,14,44"))),
  K3 = list(module = M01, expect = list(eq("ANALYSIS_BY", "PseudoBulk"), eq("CONDITIONS_TO_EXCLUDE", ""), eq("LOG2FC_INF_AS_NA", "true"),
                                        eq("THRESHOLD_PROB", "0.95"), has("THRESHOLD_CATEGORY_FOLDER", CAT_P75), hasnt("THRESHOLD_CATEGORY_FOLDER", "TrkB"),
                                        eq("CATEGORIES_TO_ANALYSE", ""))),
  P1 = list(module = M03, expect = list(eq("STAT_TABLE", "RAWvalues"), eq("TEST", "welch"), eq("COMPARE_BY", "Treatment"),
                                        eq("REFERENCE_TREATMENT", "NoTreat"), eq("REFERENCE_TIMEPOINT", "0"),
                                        eq("TREATMENTS", "NoTreat,ContDepr,BDNF"), eq("CLUSTERS", ""), eq("CATEGORIES", ""),
                                        eq("MARKERS", ""),
                                        eq("ANALYSIS_BY", "PseudoBulk", module = M01), eq("CONDITIONS_TO_EXCLUDE", "", module = M01),
                                        eq("THRESHOLD_PROB", "0.95", module = M01), atleast("Script_Version", "2", module = M01))),
  E = list(module = M01, expect = list(eq("ANALYSIS_BY", "Clusters"), eq("CLUSTER_ROUND", "1"), eq("CLUSTERS_MERGED", "false"),
                                       eq("CONDITIONS_TO_EXCLUDE", ""), eq("CLUSTERS_TO_REMOVE", ""), eq("THRESHOLD_PROB", "0.95"),
                                       eq("LOG2FC_INF_AS_NA", "true"), eq("FILL_EMPTY_CLUSTERS", "true"),
                                       eq("THRESHOLD_CATEGORY_FOLDER", ""), atleast("Script_Version", "3"))),
  P2 = list(module = M03, expect = p_cluster_spec("NoTreat,BDNF", "1,2,4,7,13")),
  P3 = list(module = M03, expect = p_cluster_spec("NoTreat,CompMed", "1,3,15,17")),
  P4 = list(module = M03, expect = p_cluster_spec("NoTreat,BDNF", "9,17,18")))
RUN_FOLDERS <- NEW_RUNS

N <- "Expt1234_PercAboveThresh_95thPerc_"     # start of the NEW file names (runs A-C)
# 48-column "RAWvalues" table (all markers) of 10.1 / 12.1: rows matched on File_Num (or File_Num + Cluster)
gp_raw48 <- function(label, rel, md5, run, byc)
  pair(label, old_rel = rel, old_md5 = md5, run = run, new = paste0(N, "RAWvalues.csv"),
       key = if (byc) c("File_Num", "Cluster") else "File_Num",
       rename_old = c(TimePoint = "TimePoints"), ignore = if (byc) character(0) else "Cluster",
       new_columns = "N_Cells", require_folders = TRUE)
# NormTo0hr tables of 10.1 / 12.1: published ones list only some markers, so extra columns in the new tables are allowed.
gp <- function(label, rel, md5, run, new, key, byc)
  pair(label, old_rel = rel, old_md5 = md5, run = run, new = paste0(N, new, ".csv"), key = key,
       rename_old = c(TimePoint = "TimePoints"), ignore = if (byc) character(0) else "Cluster",
       extra_new_ok = TRUE, dropped_columns = "Replicate", require_folders = TRUE)
# Same for the 12.1.1 tables (label column is called Condition; Cluster is 0 for PseudoBulk; ThreshCat in category runs)
g121_raw <- function(label, rel, md5, run, prefix, key)
  pair(label, old_rel = rel, old_md5 = md5, run = run, new = paste0(prefix, "RAWvalues.csv"), key = key,
       rename_old = c(Sample_Condition = "Condition"), new_columns = "N_Cells", require_folders = TRUE)
g121 <- function(label, rel, md5, run, new, prefix, key, cat)
  pair(label, old_rel = rel, old_md5 = md5, run = run, new = paste0(prefix, new, ".csv"), key = key,
       rename_old = c(Sample_Condition = "Condition"), extra_new_ok = TRUE, dropped_columns = "Replicate",
       require_folders = TRUE)
# ThresholdCategory.csv: one column, one row per cell, compared row by row
gcat <- function(label, rel, md5, run)
  pair(label, old_rel = rel, old_md5 = md5, run = run, new = "ThresholdCategory.csv", key = NA,
       rename_old = c(ThresholdCategory = 'thresholding.ct.df.category[, "ThresholdCategory"]'),
       require_folders = TRUE)
# Published p-values (Supplementary Table 3, made into long form; see GoldenData/README_GoldenData.txt)
gstat <- function(label, file, md5, run, tolerance = 1e-7)
  pair(label, old_rel = file, old_md5 = md5, old_root = "stats", run = run, new = "Statistics_long.csv",
       key = c("Marker", "TimePoint", "Group"), extra_new_ok = TRUE, tolerance = tolerance, extra_new_rows_by = "Marker",
       require_folders = TRUE)
PAIRS <- list(
  gp_raw48("A: 10.1 RAWvalues, 45 samples (Fig SuppFig5 table)",
     "Figure 2 - PseudoBulk/10.1-2026-03-03_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4__InitialDataset_Everything/Expt1,2,3,4_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_RAWvalues.csv",
     "b2681972e60a8af82c6beb71e9112ddf", run = "A", byc = FALSE),
  gp("A: 10.1 45 samples RAWvalues",
     "Figure 2 - PseudoBulk/10.1-2026-03-03_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4__InitialDataset_Everything/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "e950bf80c8e8ddf7856cc0f49931c9fe", run = "A", new = "RAWvalues", key = c("Sample_Condition", "Replicate"), byc = FALSE),
  gp("A: 10.1 45 samples RAWvalues_RepsAvg",
     "Figure 2 - PseudoBulk/10.1-2026-03-03_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4__InitialDataset_Everything/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "c7c8b73446c4237405e965829778a92f", run = "A", new = "RAWvalues_RepsAvg", key = "Sample_Condition", byc = FALSE),
  gp("A: 10.1 45 samples FC",
     "Figure 2 - PseudoBulk/10.1-2026-03-03_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4__InitialDataset_Everything/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC.csv",
     "f6b71d1e72945703095358fcec0573db", run = "A", new = "FC", key = c("Sample_Condition", "Replicate"), byc = FALSE),
  gp("A: 10.1 45 samples FC_RepsAvg",
     "Figure 2 - PseudoBulk/10.1-2026-03-03_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4__InitialDataset_Everything/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "4aaf5943cf9eac73b814a2c3b62b74ae", run = "A", new = "FC_RepsAvg", key = "Sample_Condition", byc = FALSE),
  gp("A: 10.1 45 samples log2FC",
     "Figure 2 - PseudoBulk/10.1-2026-03-03_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4__InitialDataset_Everything/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "cf640ae1276322d88fd0f5d78f982f63", run = "A", new = "log2FC", key = c("Sample_Condition", "Replicate"), byc = FALSE),
  gp("A: 10.1 45 samples log2FC_RepsAvg",
     "Figure 2 - PseudoBulk/10.1-2026-03-03_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4__InitialDataset_Everything/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "ed4cf4ecb93bb1384f5fe310a38c5c0d", run = "A", new = "log2FC_RepsAvg", key = "Sample_Condition", byc = FALSE),
  gp_raw48("B: 10.1 RAWvalues, 32 samples (Fig 3C)",
     "Figure 2 - PseudoBulk/10.1-2026-02-24_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4_BDNF_InitialDataset_Everything (Fig3C)/Expt1,2,3,4_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_RAWvalues.csv",
     "55e96a882665f2ac40e71306cd4c2a11", run = "B", byc = FALSE),
  gp("B: Fig 3C RAWvalues",
     "Figure 2 - PseudoBulk/10.1-2026-02-24_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4_BDNF_InitialDataset_Everything (Fig3C)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "573dbeb1e639d0b9867a6c449a53a44e", run = "B", new = "RAWvalues", key = c("Sample_Condition", "Replicate"), byc = FALSE),
  gp("B: Fig 3C RAWvalues_RepsAvg",
     "Figure 2 - PseudoBulk/10.1-2026-02-24_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4_BDNF_InitialDataset_Everything (Fig3C)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "5e8b3b9775458178ab4b1030f0c10ce2", run = "B", new = "RAWvalues_RepsAvg", key = "Sample_Condition", byc = FALSE),
  gp("B: Fig 3C FC",
     "Figure 2 - PseudoBulk/10.1-2026-02-24_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4_BDNF_InitialDataset_Everything (Fig3C)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC.csv",
     "77c77ded08e9af8b594e5612e39fa65f", run = "B", new = "FC", key = c("Sample_Condition", "Replicate"), byc = FALSE),
  gp("B: Fig 3C FC_RepsAvg",
     "Figure 2 - PseudoBulk/10.1-2026-02-24_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4_BDNF_InitialDataset_Everything (Fig3C)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "9e91615d3728b7c1626054efff34e19f", run = "B", new = "FC_RepsAvg", key = "Sample_Condition", byc = FALSE),
  gp("B: Fig 3C log2FC",
     "Figure 2 - PseudoBulk/10.1-2026-02-24_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4_BDNF_InitialDataset_Everything (Fig3C)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "c718c84d9dd2401b8dd1d4b02ac2ce99", run = "B", new = "log2FC", key = c("Sample_Condition", "Replicate"), byc = FALSE),
  gp("B: Fig 3C log2FC_RepsAvg",
     "Figure 2 - PseudoBulk/10.1-2026-02-24_PseudoBulk_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4_BDNF_InitialDataset_Everything (Fig3C)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "8e378708ecc9bbccf04990f8f3327873", run = "B", new = "log2FC_RepsAvg", key = "Sample_Condition", byc = FALSE),
  gp_raw48("C: 12.1 RAWvalues, by cluster (SuppFig6)",
     "Figure 2 - PseudoBulk/08-2026-02-24_MergeClusterAnalysis (SuppFig6)/12.1-2026-02-24_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_RAWvalues.csv",
     "ac96ce8e089607d89c27facd4ea1b5a5", run = "C", byc = TRUE),
  gp("C: 12.1 RAWvalues",
     "Figure 2 - PseudoBulk/08-2026-02-24_MergeClusterAnalysis (SuppFig6)/12.1-2026-02-24_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "6acc698c7246c9c3a8b7edc6b6f69cbb", run = "C", new = "RAWvalues", key = c("Sample_Condition", "Replicate", "Cluster"), byc = TRUE),
  gp("C: 12.1 RAWvalues_RepsAvg",
     "Figure 2 - PseudoBulk/08-2026-02-24_MergeClusterAnalysis (SuppFig6)/12.1-2026-02-24_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "d2800ac3a1a79ee69c55cd17010efbd8", run = "C", new = "RAWvalues_RepsAvg", key = c("Sample_Condition", "Cluster"), byc = TRUE),
  gp("C: 12.1 FC",
     "Figure 2 - PseudoBulk/08-2026-02-24_MergeClusterAnalysis (SuppFig6)/12.1-2026-02-24_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC.csv",
     "ee49b81d75e0f55c55a08d0c413d5e33", run = "C", new = "FC", key = c("Sample_Condition", "Replicate", "Cluster"), byc = TRUE),
  gp("C: 12.1 FC_RepsAvg",
     "Figure 2 - PseudoBulk/08-2026-02-24_MergeClusterAnalysis (SuppFig6)/12.1-2026-02-24_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "c5b34cd8447a0d7d680898dcfa0bcd97", run = "C", new = "FC_RepsAvg", key = c("Sample_Condition", "Cluster"), byc = TRUE),
  gp("C: 12.1 log2FC",
     "Figure 2 - PseudoBulk/08-2026-02-24_MergeClusterAnalysis (SuppFig6)/12.1-2026-02-24_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "78a12c331cfd2133fe2324ee727942bc", run = "C", new = "log2FC", key = c("Sample_Condition", "Replicate", "Cluster"), byc = TRUE),
  gp("C: 12.1 log2FC_RepsAvg",
     "Figure 2 - PseudoBulk/08-2026-02-24_MergeClusterAnalysis (SuppFig6)/12.1-2026-02-24_PercAboveThresh_95thPercOfNoTreat_Expt1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "026266b11dbb666cf1fcf43e1b5b22ea", run = "C", new = "log2FC_RepsAvg", key = c("Sample_Condition", "Cluster"), byc = TRUE),
  g121_raw("D95: 12.1.1 95th percentile RAWvalues, 45 samples",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (Fig3D, SuppFig3)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_RAWvalues.csv",
     "b1a0813b77c8ef40dd04afc868a11188", run = "D95", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = "File_Num"),
  g121("D95: 12.1.1 95th percentile RAWvalues",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (Fig3D, SuppFig3)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "c22fe264e9863093b92be56d4a598488", run = "D95", new = "RAWvalues", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D95: 12.1.1 95th percentile RAWvalues_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (Fig3D, SuppFig3)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "434daf820ffe896f82e50da69fcfb9bd", run = "D95", new = "RAWvalues_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = "Sample_Condition", cat = FALSE),
  g121("D95: 12.1.1 95th percentile FC",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (Fig3D, SuppFig3)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC.csv",
     "d6771e78ceb4dae20e118c409e751a7e", run = "D95", new = "FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D95: 12.1.1 95th percentile FC_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (Fig3D, SuppFig3)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "1fb4e13753df00aa2efcb8ba590ff2e9", run = "D95", new = "FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = "Sample_Condition", cat = FALSE),
  g121("D95: 12.1.1 95th percentile log2FC",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (Fig3D, SuppFig3)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "86a9f7938dec36456148718cea091c5e", run = "D95", new = "log2FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D95: 12.1.1 95th percentile log2FC_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (Fig3D, SuppFig3)/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "1328d4ddb5e44e7a3fbfc5e106860cd2", run = "D95", new = "log2FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = "Sample_Condition", cat = FALSE),
  g121_raw("D75: 12.1.1 75th percentile RAWvalues, 45 samples",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_75thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_75thPercOfNoTreat_RAWvalues.csv",
     "a7e4fe6b82e463e13f8d457361ebabf4", run = "D75", prefix = "Expt1234_PercAboveThresh_75thPerc_", key = "File_Num"),
  g121("D75: 12.1.1 75th percentile RAWvalues",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_75thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_75thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "58973dc46d26696093b4324a0ce911d5", run = "D75", new = "RAWvalues", prefix = "Expt1234_PercAboveThresh_75thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D75: 12.1.1 75th percentile RAWvalues_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_75thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_75thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "ee7adb05df36c7e65e8774cd51306028", run = "D75", new = "RAWvalues_RepsAvg", prefix = "Expt1234_PercAboveThresh_75thPerc_", key = "Sample_Condition", cat = FALSE),
  g121("D75: 12.1.1 75th percentile FC",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_75thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_75thPercOfNoTreat_NormTo0hr_FC.csv",
     "d05897457d9c8aaefb2da02521cca002", run = "D75", new = "FC", prefix = "Expt1234_PercAboveThresh_75thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D75: 12.1.1 75th percentile FC_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_75thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_75thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "4a1c6278b769f2769d6d5a7078889a89", run = "D75", new = "FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_75thPerc_", key = "Sample_Condition", cat = FALSE),
  g121("D75: 12.1.1 75th percentile log2FC",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_75thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_75thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "ab417e0f16ddbb40c265c58aebbb7557", run = "D75", new = "log2FC", prefix = "Expt1234_PercAboveThresh_75thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D75: 12.1.1 75th percentile log2FC_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_75thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_75thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "87454fcac0f9572ba2ba066c195a6c13", run = "D75", new = "log2FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_75thPerc_", key = "Sample_Condition", cat = FALSE),
  g121_raw("D50: 12.1.1 50th percentile RAWvalues, 45 samples",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_50thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_50thPercOfNoTreat_RAWvalues.csv",
     "3558dfcceaf3989f119e3b78bb2252d9", run = "D50", prefix = "Expt1234_PercAboveThresh_50thPerc_", key = "File_Num"),
  g121("D50: 12.1.1 50th percentile RAWvalues",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_50thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_50thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "617bd55c625473189320d10348c40d01", run = "D50", new = "RAWvalues", prefix = "Expt1234_PercAboveThresh_50thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D50: 12.1.1 50th percentile RAWvalues_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_50thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_50thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "bbb1a673037162f64996077d4c71b11f", run = "D50", new = "RAWvalues_RepsAvg", prefix = "Expt1234_PercAboveThresh_50thPerc_", key = "Sample_Condition", cat = FALSE),
  g121("D50: 12.1.1 50th percentile FC",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_50thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_50thPercOfNoTreat_NormTo0hr_FC.csv",
     "d353e4e97bc30da06f5485821aa96060", run = "D50", new = "FC", prefix = "Expt1234_PercAboveThresh_50thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D50: 12.1.1 50th percentile FC_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_50thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_50thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "133231725d595b09a314832e8739d020", run = "D50", new = "FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_50thPerc_", key = "Sample_Condition", cat = FALSE),
  g121("D50: 12.1.1 50th percentile log2FC",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_50thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_50thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "5cee9df55946dfe052222ad57e2eab7d", run = "D50", new = "log2FC", prefix = "Expt1234_PercAboveThresh_50thPerc_", key = c("Sample_Condition", "Replicate"), cat = FALSE),
  g121("D50: 12.1.1 50th percentile log2FC_RepsAvg",
     "Figure 2 - PseudoBulk/12.1.1-2026-02-25_DotPlot_PseudoBulk_50thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers (SuppFig3)/Expt1,2,3,4_PercAboveThresh_50thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "46da27d92d5c79ac68cfa20619b8bd7e", run = "D50", new = "log2FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_50thPerc_", key = "Sample_Condition", cat = FALSE),
  gcat("T1: ThresholdCategory.csv, TrkB_i only (cuts 0.15, 1.4, 3.7)",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7/ThresholdCategory.csv",
     "1548f018290dbb852d559ea6682dbb8e", run = "T1"),
  gcat("T2: ThresholdCategory.csv, p75NTR only (cuts 0.45, 1.45, 3.95)",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95/ThresholdCategory.csv",
     "936f8509f01d54c6fa7e25e93d5dbd38", run = "T2"),
  gcat("T3: ThresholdCategory.csv, p75NTR + TrkB_i (Fig 6A)",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (Fig6A, SuppFig11)/ThresholdCategory.csv",
     "224fb10cd915925407edc0ff9e1edb86", run = "T3"),
  gcat("T4: ThresholdCategory.csv, p75NTR cut 1.45 + TrkB_i (SuppFig12)",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr1.45_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (SuppFig12)/ThresholdCategory.csv",
     "2037a5d4d317c009c082d652bc6d7b2d", run = "T4"),
  g121_raw("K1: Fig 5 TrkB_i categories 1-4, RAWvalues (all markers)",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_RAWvalues-1,2,3,4.csv",
     "2d32e394967ea65ee04475e0c0e7301e", run = "K1", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("File_Num", "ThreshCat")),
  g121("K1: Fig 5 TrkB_i categories 1-4, RAWvalues",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "d52049e40cbc71cc9cb15f351fe4ff3d", run = "K1", new = "RAWvalues", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K1: Fig 5 TrkB_i categories 1-4, RAWvalues_RepsAvg",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "2674a814e97f0a508fb89ddca8cd481f", run = "K1", new = "RAWvalues_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  g121("K1: Fig 5 TrkB_i categories 1-4, FC",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC.csv",
     "b58de4d6aa27a09966b8f3ab878f33a5", run = "K1", new = "FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K1: Fig 5 TrkB_i categories 1-4, FC_RepsAvg",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "2df6b3fe286f1aa25e9f8437ee741854", run = "K1", new = "FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  g121("K1: Fig 5 TrkB_i categories 1-4, log2FC",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "dbed5d8d2fba26c16226f04f1733e8c0", run = "K1", new = "log2FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K1: Fig 5 TrkB_i categories 1-4, log2FC_RepsAvg",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "35f4d61c733f531be1eea72e1629b9de", run = "K1", new = "log2FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  g121_raw("K2: Fig 6A categories 11,41,14,44, RAWvalues (all markers)",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (Fig6A, SuppFig11)/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-11,41,14,44/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_RAWvalues-11,41,14,44.csv",
     "915c0dd3350d6956b5b399e22a121eb8", run = "K2", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("File_Num", "ThreshCat")),
  g121("K2: Fig 6A categories 11,41,14,44, RAWvalues",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (Fig6A, SuppFig11)/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-11,41,14,44/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "95665f18c1f78debc72539552e3398a7", run = "K2", new = "RAWvalues", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K2: Fig 6A categories 11,41,14,44, RAWvalues_RepsAvg",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (Fig6A, SuppFig11)/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-11,41,14,44/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "63426137b6a0b7c551a1cd6d346e7d2c", run = "K2", new = "RAWvalues_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  g121("K2: Fig 6A categories 11,41,14,44, FC",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (Fig6A, SuppFig11)/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-11,41,14,44/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC.csv",
     "b0c33008647cfa8ca94ed015df5d454f", run = "K2", new = "FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K2: Fig 6A categories 11,41,14,44, FC_RepsAvg",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (Fig6A, SuppFig11)/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-11,41,14,44/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "2348d5476039d2814a5d3d2a91ce3d96", run = "K2", new = "FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  g121("K2: Fig 6A categories 11,41,14,44, log2FC",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (Fig6A, SuppFig11)/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-11,41,14,44/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "d925baf70d5b8602f6612af381d6a9c9", run = "K2", new = "log2FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K2: Fig 6A categories 11,41,14,44, log2FC_RepsAvg",
     "Figure 6 - Co-expression/13.1-2026-02-28_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7 (Fig6A, SuppFig11)/12.1.1-2026-02-28_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-11,41,14,44/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "1be2f29b8af75219b36849512550b5c1", run = "K2", new = "log2FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  g121_raw("K3: p75NTR categories 1-4 (02-27 folder), RAWvalues (all markers)",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95/12.1.1-2026-02-27_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_RAWvalues-1,2,3,4.csv",
     "58e42fd12ae7770ee404abc4d806dabd", run = "K3", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("File_Num", "ThreshCat")),
  g121("K3: p75NTR categories 1-4 (02-27 folder), RAWvalues",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95/12.1.1-2026-02-27_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "af90f5cfef82d43fae3932493f927c04", run = "K3", new = "RAWvalues", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K3: p75NTR categories 1-4 (02-27 folder), RAWvalues_RepsAvg",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95/12.1.1-2026-02-27_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "804c7dc1143c4039bed563bd6f054598", run = "K3", new = "RAWvalues_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  g121("K3: p75NTR categories 1-4 (02-27 folder), FC",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95/12.1.1-2026-02-27_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC.csv",
     "ab274980f415be44db57b4829c34aa58", run = "K3", new = "FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K3: p75NTR categories 1-4 (02-27 folder), FC_RepsAvg",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95/12.1.1-2026-02-27_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "19d0c59635e72a263d5015868e94974e", run = "K3", new = "FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  g121("K3: p75NTR categories 1-4 (02-27 folder), log2FC",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95/12.1.1-2026-02-27_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "665bc7ef361db89979d696f33fce05b5", run = "K3", new = "log2FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "ThreshCat"), cat = TRUE),
  g121("K3: p75NTR categories 1-4 (02-27 folder), log2FC_RepsAvg",
     "Figure 5 - TrkB focus/13.1-2026-02-27_ExtractByMarkerExpression_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95/12.1.1-2026-02-27_DotPlot_PseudoBulk_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers-1,2,3,4/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "1f220f9a8f40925278334f6fd5aff76d", run = "K3", new = "log2FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "ThreshCat"), cat = TRUE),
  gstat("P1: Supplementary Table 3 p-values (Fig 2C / SuppFig 4B), 20 markers x 11 groups",
     "PublishedStats_TableS3_long.csv", "bb48eadcf5a620cc9578a2a20546f4f7", run = "P1"),
  g121_raw("E: 12.1.1 Clusters RAWvalues, all 40 markers, 45 samples x 20 clusters (900 rows)",
     "Figure 3 - Cell type analysis/12.1.1-2026-03-01_DotPlot_Clusters_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_RAWvalues.csv",
     "921147404b4562f254a255f1bcf5ec8f", run = "E", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("File_Num", "Cluster")),
  g121("E: 12.1.1 Clusters RAWvalues (22 markers)",
     "Figure 3 - Cell type analysis/12.1.1-2026-03-01_DotPlot_Clusters_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues.csv",
     "b711b4e10d3543fb332cbac81885687d", run = "E", new = "RAWvalues", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "Cluster"), cat = FALSE),
  g121("E: 12.1.1 Clusters RAWvalues_RepsAvg",
     "Figure 3 - Cell type analysis/12.1.1-2026-03-01_DotPlot_Clusters_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_RAWvalues_RepsAvg.csv",
     "077959b1e1e7813cd8e2f3b4252986a3", run = "E", new = "RAWvalues_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Cluster"), cat = FALSE),
  g121("E: 12.1.1 Clusters FC",
     "Figure 3 - Cell type analysis/12.1.1-2026-03-01_DotPlot_Clusters_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC.csv",
     "2205bfedfe5f876c72a1a4b4f68b2138", run = "E", new = "FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "Cluster"), cat = FALSE),
  g121("E: 12.1.1 Clusters FC_RepsAvg",
     "Figure 3 - Cell type analysis/12.1.1-2026-03-01_DotPlot_Clusters_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_FC_RepsAvg.csv",
     "7227f9f60f96f14450ed034b45dd66e6", run = "E", new = "FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Cluster"), cat = FALSE),
  g121("E: 12.1.1 Clusters log2FC",
     "Figure 3 - Cell type analysis/12.1.1-2026-03-01_DotPlot_Clusters_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC.csv",
     "ee32a89c39b1e1dc599eb9619b7232b0", run = "E", new = "log2FC", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Replicate", "Cluster"), cat = FALSE),
  g121("E: 12.1.1 Clusters log2FC_RepsAvg",
     "Figure 3 - Cell type analysis/12.1.1-2026-03-01_DotPlot_Clusters_95thPercOfNoTreat_Expt1,2,3,4_SizeByFracAbund_ColorByFoldChange_InclAllClusts_SelectMarkers/Expt1,2,3,4_PercAboveThresh_95thPercOfNoTreat_NormTo0hr_log2FC_RepsAvg.csv",
     "82ff9a3aab2e9a161a22aa0de85fc969", run = "E", new = "log2FC_RepsAvg", prefix = "Expt1234_PercAboveThresh_95thPerc_", key = c("Sample_Condition", "Cluster"), cat = FALSE),
  gstat("P2: Supplementary Table 4A p-values (Fig 3C, NoTreat + BDNF, clusters 1,2,4,7,13), 5 markers x 5 clusters x 6 time points",
     "PublishedStats_TableS4A_long.csv", "9033afabf7d33fd6abfe61a965ca95f2", run = "P2"),
  gstat("P3: Supplementary Table 4B p-values (Fig 3D, NoTreat + CompMed, clusters 1,3,15,17), 4 markers x 4 clusters x 6 time points",
     "PublishedStats_TableS4B_long.csv", "4af228c2555111427f1e518716d20f4b", run = "P3",
     tolerance = 5e-4),  # 2 published cells (0.0889, 0.258) carry only 3 significant figures; the other 94 are 7-digit
  gstat("P4: Supplementary Table 4C p-values (Fig 3E, NoTreat + BDNF, clusters 9,17,18), 5 markers x 3 clusters x 6 time points",
     "PublishedStats_TableS4C_long.csv", "8f95860df645cbd9b352abf093339e5c", run = "P4")
)
# --- Example B: 01_Get_Metadata (old files in OLD_FOLDER, new ones in NEW_FOLDER)
# PAIRS <- list(pair("metadata", "metadata.csv", "metadata.csv", key = "File_Name"),
#               pair("panel",    "panel.csv",    "panel.csv"))



#### Set up ####################################################################
choose_folder <- function(title) {
  if (!requireNamespace("svDialogs", quietly = TRUE)) {
    stop("Please set OLD_FOLDER and NEW_FOLDER at the top of this script, or install ",
         "the 'svDialogs' package: install.packages(\"svDialogs\")", call. = FALSE)
  }
  res <- svDialogs::dlg_dir(default = getwd(), title = title)$res
  if (length(res) == 0 || !nzchar(res)) stop("No folder selected.", call. = FALSE)
  res
}
needs_old <- any(vapply(PAIRS, function(p) !p$require_folders && !nzchar(p$old_folder), logical(1)))
needs_new <- any(vapply(PAIRS, function(p) !p$require_folders && !nzchar(p$new_folder), logical(1)))
if (needs_old && !nzchar(OLD_FOLDER)) OLD_FOLDER <- choose_folder("Select the folder with the OLD (published) files")
if (needs_new && !nzchar(NEW_FOLDER)) NEW_FOLDER <- choose_folder("Select the folder with the NEW script's files")
if (nzchar(OLD_FOLDER) && !dir.exists(OLD_FOLDER)) stop("OLD_FOLDER not found: ", OLD_FOLDER, call. = FALSE)
if (nzchar(NEW_FOLDER) && !dir.exists(NEW_FOLDER)) stop("NEW_FOLDER not found: ", NEW_FOLDER, call. = FALSE)
uses_golden <- any(vapply(PAIRS, function(p) nzchar(p$old_rel) && !identical(p$old_root, "stats") && nzchar(RUN_FOLDERS[[p$run]]), logical(1)))
if (uses_golden && !nzchar(GOLDEN_ROOT))
  stop("Please set GOLDEN_ROOT (the \"2026-02-23 Figures_v12,13\" folder) at the top of this script.", call. = FALSE)
if (nzchar(GOLDEN_ROOT) && !dir.exists(GOLDEN_ROOT)) stop("GOLDEN_ROOT not found: ", GOLDEN_ROOT, call. = FALSE)
uses_stats <- any(vapply(PAIRS, function(p) identical(p$old_root, "stats") && nzchar(RUN_FOLDERS[[p$run]]), logical(1)))
if (uses_stats && !nzchar(GOLDEN_STATS_FOLDER))
  stop("Please set GOLDEN_STATS_FOLDER (the Tools/GoldenData folder) at the top of this script.", call. = FALSE)
if (nzchar(GOLDEN_STATS_FOLDER) && !dir.exists(GOLDEN_STATS_FOLDER)) stop("GOLDEN_STATS_FOLDER not found: ", GOLDEN_STATS_FOLDER, call. = FALSE)
for (rn in names(RUN_FOLDERS)) if (nzchar(RUN_FOLDERS[[rn]]) && !dir.exists(RUN_FOLDERS[[rn]]))
  stop("NEW_RUNS$", rn, " folder not found: ", RUN_FOLDERS[[rn]], call. = FALSE)
if (!nzchar(REPORT_FOLDER)) {
  rf <- unlist(RUN_FOLDERS); rf <- rf[nzchar(rf)]
  REPORT_FOLDER <- if (nzchar(NEW_FOLDER)) NEW_FOLDER else if (length(rf)) rf[[1]] else getwd()
}

# Read what a NEW run used, from the run_parameters.json inside its folder
read_run_settings <- function(folder, module) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Please install the 'jsonlite' package: install.packages(\"jsonlite\")", call. = FALSE)
  }
  f <- file.path(folder, "run_parameters.json")
  if (!file.exists(f)) return(NULL)
  rec <- tryCatch(jsonlite::read_json(f, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(rec) || is.null(rec[[module]])) return(list(found = FALSE))
  sec <- rec[[module]]
  flat <- function(v) {
    if (is.null(v)) return(NA_character_)
    v <- unlist(v)
    if (length(v) == 0) return("")
    if (is.logical(v)) return(paste(tolower(as.character(v)), collapse = ","))
    paste(as.character(v), collapse = ",")
  }
  get <- function(key) {
    if (key == "FOLDER_NAME") return(basename(folder))
    if (key %in% names(sec) && !is.list(sec[[key]])) return(flat(sec[[key]]))
    flat(sec[["Settings"]][[key]])
  }
  list(found = TRUE, get = get)
}
# Check one run folder against RUN_SPEC. Returns list(ok, lines).
check_run <- function(rn) {
  folder <- RUN_FOLDERS[[rn]]; spec <- RUN_SPEC[[rn]]
  out <- c(paste0("RUN ", rn, " folder: ", folder))
  st_cache <- list()
  st_for <- function(module) {
    if (is.null(st_cache[[module]])) st_cache[[module]] <<- read_run_settings(folder, module)
    st_cache[[module]]
  }
  st <- st_for(spec$module)
  if (is.null(st)) return(list(ok = FALSE, lines = c(out, "  !! No run_parameters.json in this folder: it is not an output folder of the right script.")))
  if (!isTRUE(st$found)) return(list(ok = FALSE, lines = c(out, paste0("  !! run_parameters.json has no '", spec$module, "' section: this folder was not made by that script."))))
  ok <- TRUE
  out <- c(out, paste0("  module ", spec$module, " found"))
  for (e in spec$expect) {
    est <- if (is.null(e$module)) st else st_for(e$module)
    if (!isTRUE(est$found)) {
      ok <- FALSE
      out <- c(out, paste0("  !! ", e$key, ": the run record has no '", e$module, "' section (the input folder was not made by that script)"))
      next
    }
    v <- est$get(e$key)
    good <- switch(e$op,
                   eq    = identical(unname(v), unname(e$value)),
                   has   = !is.na(v) && grepl(e$value, v, fixed = TRUE),
                   hasnt = !is.na(v) && !grepl(e$value, v, fixed = TRUE),
                   atleast = { vn <- suppressWarnings(as.numeric(v)); !is.na(vn) && vn >= as.numeric(e$value) })
    if (!good) ok <- FALSE
    shown <- if (is.na(v)) "(not recorded)" else paste0("'", v, "'")
    out <- c(out, paste0("  ", if (good) "OK  " else "!!  ", if (!is.null(e$module)) paste0("[", e$module, "] ") else "", e$key, " = ", shown, "   expected ",
                         switch(e$op, eq = "exactly ", has = "to contain ", hasnt = "NOT to contain ", atleast = "at least "), "'", e$value, "'"))
  }
  list(ok = ok, lines = out)
}

full_path <- function(f, folder) {
  if (grepl("^([A-Za-z]:[/\\\\]|/|~)", f)) f else file.path(folder, f)
}

report <- character(0)
say <- function(...) {
  line <- paste0(...)
  report <<- c(report, line)
  cat(line, "\n", sep = "")
}

read_as_text <- function(path) {
  d <- read.csv(path, colClasses = "character", check.names = FALSE,
                na.strings = if (NAN_AS_NA) c("NA", "", "NaN") else c("NA", ""),
                stringsAsFactors = FALSE, comment.char = "")
  # drop the row-number column that write.csv adds (empty or "X" heading, values 1..n)
  if (ncol(d) > 1 && colnames(d)[1] %in% c("", "X") &&
      identical(suppressWarnings(as.integer(d[[1]])), seq_len(nrow(d)))) {
    d <- d[, -1, drop = FALSE]
  }
  d
}

# TRUE where two text vectors hold the same value (numbers compared as numbers)
same_value <- function(a, b, tol) {
  a_num <- suppressWarnings(as.numeric(a))
  b_num <- suppressWarnings(as.numeric(b))
  both_na   <- is.na(a) & is.na(b)
  both_num  <- !is.na(a_num) & !is.na(b_num)
  num_equal <- both_num & ((a_num == b_num) | (abs(a_num - b_num) <= tol))
  num_equal[is.na(num_equal)] <- FALSE
  txt_equal <- !is.na(a) & !is.na(b) & (a == b)
  both_na | num_equal | txt_equal
}

key_string <- function(d, key) do.call(paste, c(lapply(key, function(k) d[[k]]), sep = " | "))



#### One comparison ############################################################
# Returns list(lines, stats). `perturb = TRUE` is the negative control: one value
# in the NEW table is changed before comparing.
compare_pair <- function(p, perturb = FALSE) {
  lines <- character(0)
  add <- function(...) lines <<- c(lines, paste0(...))
  stats <- list(rows_old = NA, rows_new = NA, rows_matched = NA, values_compared = 0,
                max_abs_diff = 0, n_unexpected = 0, n_expected = 0, structural = 0,
                md5_old = NA, md5_new = NA, md5_expected = NA, md5_ok = NA)

  golden <- nzchar(p$old_rel)
  tol <- if (is.na(p$tolerance)) NUMERIC_TOLERANCE else p$tolerance
  stats$tol <- tol
  root <- if (identical(p$old_root, "stats")) GOLDEN_STATS_FOLDER else GOLDEN_ROOT
  root_name <- if (identical(p$old_root, "stats")) "GOLDEN_STATS_FOLDER" else "GOLDEN_ROOT"
  old_dir <- if (nzchar(p$old_folder)) p$old_folder else OLD_FOLDER
  new_dir <- if (golden) RUN_FOLDERS[[p$run]] else if (nzchar(p$new_folder)) p$new_folder else NEW_FOLDER
  old_path <- if (golden) file.path(root, p$old_rel) else full_path(p$old, old_dir)
  new_path <- full_path(p$new, new_dir)
  if (golden && file.exists(old_path)) {
    stats$md5_old <- unname(tools::md5sum(old_path))
    stats$md5_expected <- p$old_md5
    stats$md5_ok <- identical(stats$md5_old, p$old_md5)
    if (!stats$md5_ok) {
      add("  !! WRONG OLD FILE: the file found at the expected location is not the published one.")
      add("     expected MD5 ", p$old_md5)
      add("     found    MD5 ", stats$md5_old)
      add("     (expected location: ", root_name, "/", p$old_rel, ")")
      add("  -- not compared")
      stats$structural <- 1; stats$not_compared <- TRUE
      return(list(lines = lines, stats = stats))
    }
    add("  OLD file MD5 matches the expected one (", substr(stats$md5_old, 1, 8), "): this is the published file")
  }
  if (!file.exists(old_path) || !file.exists(new_path)) {
    which_missing <- c(if (!file.exists(old_path)) old_path, if (!file.exists(new_path)) new_path)
    for (mp in which_missing) {
      add("  !! File not found: ", mp)
      if (golden && identical(mp, old_path)) add("     Expected inside ", root_name, " at: ", p$old_rel)
      d <- dirname(mp)
      if (dir.exists(d)) {
        near <- list.files(d, pattern = "\\.csv$")
        add("     The folder holds ", length(near), " csv file(s)", if (length(near)) ":" else ".")
        for (nf in head(near, 12)) add("       ", nf)
      } else add("     (the folder itself does not exist)")
    }
    add("  -- not compared")
    stats$structural <- 1
    stats$not_compared <- TRUE
    return(list(lines = lines, stats = stats))
  }
  stats$md5_old <- unname(tools::md5sum(old_path))
  stats$md5_new <- unname(tools::md5sum(new_path))
  if (is.null(stats$md5_ok)) { stats$md5_ok <- NA; stats$md5_expected <- NA }
  old <- read_as_text(old_path)
  new <- read_as_text(new_path)

  # rename OLD columns
  for (nm in names(p$rename_old)) {
    from <- p$rename_old[[nm]]
    if (from != nm && from %in% colnames(old)) colnames(old)[colnames(old) == from] <- nm
  }

  # columns
  only_old <- setdiff(colnames(old), colnames(new))
  only_new <- setdiff(colnames(new), colnames(old))
  declared_new <- if (isTRUE(p$extra_new_ok)) only_new else intersect(only_new, p$new_columns)
  only_new <- setdiff(only_new, declared_new)
  declared_old <- intersect(only_old, p$dropped_columns)
  only_old <- setdiff(only_old, declared_old)
  if (length(declared_old) > 0) add("  ~  (declared) columns only in OLD (dropped on purpose): ", paste(declared_old, collapse = ", "))
  if (length(only_old) > 0) { add("  !! Columns only in OLD: ", paste(only_old, collapse = ", ")); stats$structural <- stats$structural + 1 }
  if (length(only_new) > 0) { add("  !! Columns only in NEW: ", paste(only_new, collapse = ", ")); stats$structural <- stats$structural + 1 }
  if (length(declared_new) > 0) add("  ~  (declared) columns only in NEW: ", paste(declared_new, collapse = ", "))
  if (length(p$rename_old) > 0) {
    add("  ~  OLD columns renamed: ", paste0(p$rename_old, " -> ", names(p$rename_old), collapse = "; "))
  }
  shared <- setdiff(intersect(colnames(old), colnames(new)), p$ignore)
  if (length(p$ignore) > 0) add("  ~  Not compared (ignored): ", paste(p$ignore, collapse = ", "))

  stats$rows_old <- nrow(old); stats$rows_new <- nrow(new)
  add("  Rows: OLD = ", nrow(old), ", NEW = ", nrow(new))
  if (nrow(old) != nrow(new) && !(nzchar(p$extra_new_rows_by) && nrow(new) > nrow(old))) stats$structural <- stats$structural + 1

  if (perturb) {
    # negative control: change the last numeric column (not a key) by 0.001 in one row
    num_cols <- shared[vapply(shared, function(cn) !all(is.na(suppressWarnings(as.numeric(new[[cn]])))), logical(1))]
    num_cols <- setdiff(num_cols, p$key)
    if (length(num_cols) > 0 && nrow(new) > 0) {
      cn <- num_cols[length(num_cols)]
      new[[cn]][1] <- as.character(as.numeric(new[[cn]][1]) + 0.001)
    }
  }

  key <- p$key
  use_key <- !(length(key) == 1 && is.na(key)) && all(key %in% colnames(old)) && all(key %in% colnames(new))
  if (!(length(key) == 1 && is.na(key)) && !use_key) {
    add("  !! Key column(s) not found in both files: ", paste(key, collapse = ", "), "; comparing row by row instead")
    stats$structural <- stats$structural + 1
  }
  if (use_key) {
    ko <- key_string(old, key); kn <- key_string(new, key)
    if (anyDuplicated(ko) || anyDuplicated(kn)) {
      add("  !! Key column(s) '", paste(key, collapse = "+"), "' have duplicates; comparing row by row instead")
      stats$structural <- stats$structural + 1
      use_key <- FALSE
    }
  }
  if (use_key) {
    in_both  <- intersect(ko, kn)
    miss_new <- setdiff(ko, kn); miss_old <- setdiff(kn, ko)
    if (length(miss_new) > 0) { add("  !! ", length(miss_new), " row(s) only in OLD, e.g. ", paste(head(miss_new, 3), collapse = " ; ")); stats$structural <- stats$structural + 1 }
    if (length(miss_old) > 0 && nzchar(p$extra_new_rows_by) && p$extra_new_rows_by %in% colnames(old) && p$extra_new_rows_by %in% colnames(new)) {
      ecol <- p$extra_new_rows_by
      allowed <- !(new[[ecol]][match(miss_old, kn)] %in% old[[ecol]])
      if (any(allowed)) {
        add("  ~  (declared) ", sum(allowed), " row(s) only in NEW, for ", ecol, " value(s) that are not in OLD: ",
            paste(head(sort(unique(new[[ecol]][match(miss_old[allowed], kn)])), 4), collapse = ", "),
            if (length(unique(new[[ecol]][match(miss_old[allowed], kn)])) > 4) ", ..." else "",
            "  (", length(unique(new[[ecol]][match(miss_old[allowed], kn)])), " ", ecol, " value(s) in all)")
      }
      miss_old <- miss_old[!allowed]
    }
    if (length(miss_old) > 0) { add("  !! ", length(miss_old), " row(s) only in NEW, e.g. ", paste(head(miss_old, 3), collapse = " ; ")); stats$structural <- stats$structural + 1 }
    same_order <- identical(ko[ko %in% in_both], kn[kn %in% in_both])
    add("  Rows matched on ", paste(key, collapse = " + "), ": ", length(in_both),
        if (same_order) " (same row order)" else " (row order differs; not a problem)")
    old_a <- old[match(in_both, ko), , drop = FALSE]
    new_a <- new[match(in_both, kn), , drop = FALSE]
    row_label <- in_both
  } else {
    n <- min(nrow(old), nrow(new))
    old_a <- old[seq_len(n), , drop = FALSE]
    new_a <- new[seq_len(n), , drop = FALSE]
    row_label <- paste("row", seq_len(n))
    add("  Rows compared by position: ", n)
  }
  stats$rows_matched <- nrow(old_a)

  # values
  n_diff_cols <- 0
  for (cn in shared) {
    eq <- same_value(old_a[[cn]], new_a[[cn]], tol)
    stats$values_compared <- stats$values_compared + length(eq)
    on <- suppressWarnings(as.numeric(old_a[[cn]])); nn <- suppressWarnings(as.numeric(new_a[[cn]]))
    d <- abs(on - nn); d <- d[is.finite(d)]
    if (length(d) > 0) stats$max_abs_diff <- max(stats$max_abs_diff, d)
    if (all(eq)) next
    n_diff_cols <- n_diff_cols + 1
    expected <- cn %in% p$expected
    if (expected) stats$n_expected <- stats$n_expected + sum(!eq) else stats$n_unexpected <- stats$n_unexpected + sum(!eq)
    add("  ", if (expected) "~  (expected) " else "!! ", "Column '", cn, "': ",
        sum(!eq), " of ", length(eq), " value(s) differ")
    for (i in head(which(!eq), MAX_EXAMPLES)) {
      add("       ", substr(row_label[i], 1, 60), "   OLD: ",
          ifelse(is.na(old_a[[cn]][i]), "NA", old_a[[cn]][i]),
          "   NEW: ", ifelse(is.na(new_a[[cn]][i]), "NA", new_a[[cn]][i]))
    }
  }
  if (n_diff_cols == 0 && stats$structural == 0) {
    add("  OK: all ", stats$values_compared, " values identical (tolerance ", tol,
        "; largest difference ", format(stats$max_abs_diff, digits = 3), ")")
  }
  list(lines = lines, stats = stats)
}



#### Run all comparisons #######################################################
say("COMPARISON REPORT  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "   (Compare_GoldenFiles.R v3.5)")
if (nzchar(GOLDEN_ROOT)) say("GOLDEN_ROOT: ", GOLDEN_ROOT)
if (nzchar(OLD_FOLDER)) say("OLD: ", OLD_FOLDER)
if (nzchar(NEW_FOLDER)) say("NEW: ", NEW_FOLDER)
say("Tolerance: ", NUMERIC_TOLERANCE, " (a pair may set its own, shown in its result)   NaN treated as missing: ", NAN_AS_NA)
if (nzchar(GOLDEN_STATS_FOLDER)) say("GOLDEN_STATS_FOLDER: ", GOLDEN_STATS_FOLDER)
say(strrep("-", 72))

records <- list()
skipped <- character(0)
PAIRS <- Filter(function(p) {
  if (nzchar(p$old_rel)) {
    if (!nzchar(RUN_FOLDERS[[p$run]])) { skipped <<- c(skipped, p$label); return(FALSE) }
    return(TRUE)
  }
  has_old <- nzchar(p$old_folder) || (!p$require_folders && nzchar(OLD_FOLDER)) || grepl("^([A-Za-z]:[/\\\\]|/|~)", p$old)
  has_new <- nzchar(p$new_folder) || (!p$require_folders && nzchar(NEW_FOLDER)) || grepl("^([A-Za-z]:[/\\\\]|/|~)", p$new)
  if (!(has_old && has_new)) skipped <<- c(skipped, p$label)
  has_old && has_new
}, PAIRS)
run_ok <- list()
for (rn in names(RUN_FOLDERS)) {
  if (!nzchar(RUN_FOLDERS[[rn]])) { say("RUN ", rn, ": blank in NEW_RUNS -> not checked"); next }
  chk <- check_run(rn); run_ok[[rn]] <- chk$ok
  for (l in chk$lines) say(l)
  if (!chk$ok) say("  !! RUN ", rn, " SETTINGS DO NOT MATCH: wrong folder (or a different setting was used). Its pairs are not compared.")
}
if (length(skipped) > 0) say("Skipped (folder not given for that run): ", paste(skipped, collapse = "; "))
for (p in PAIRS) {
  say("")
  say("PAIR: ", p$label)
  say("  OLD file: ", p$old)
  say("  NEW file: ", p$new)
  if (nzchar(p$old_rel) && isFALSE(run_ok[[p$run]])) {
    res <- list(lines = "  !! not compared: the NEW run folder did not pass the settings check above",
                stats = list(rows_old = NA, rows_new = NA, rows_matched = NA, values_compared = 0, max_abs_diff = 0,
                             n_unexpected = 0, n_expected = 0, structural = 1, md5_old = NA, md5_new = NA,
                             md5_expected = NA, md5_ok = NA))
  } else res <- compare_pair(p)
  for (l in res$lines) say(l)
  s <- res$stats
  status <- if (s$structural == 0 && s$n_unexpected == 0) "PASS" else "FAIL"
  say("  => ", status)
  records[[length(records) + 1]] <- data.frame(
    Pair = p$label, Old_File = p$old, New_File = p$new,
    Run = if (nzchar(p$run)) p$run else NA, Old_Rel_Path = if (nzchar(p$old_rel)) p$old_rel else NA,
    Old_MD5_Expected = s$md5_expected, Old_MD5_Is_Published_File = s$md5_ok,
    Run_Settings_OK = if (nzchar(p$run)) isTRUE(run_ok[[p$run]]) else NA,
    New_Run_Folder = if (nzchar(p$run)) RUN_FOLDERS[[p$run]] else NA,
    Old_MD5 = s$md5_old, New_MD5 = s$md5_new,
    Rows_Old = s$rows_old, Rows_New = s$rows_new, Rows_Matched = s$rows_matched,
    Values_Compared = s$values_compared, Max_Abs_Difference = s$max_abs_diff,
    Unexpected_Differences = s$n_unexpected, Expected_Differences = s$n_expected,
    Structural_Problems = s$structural, Status = status,
    Tolerance = if (is.null(s$tol)) NUMERIC_TOLERANCE else s$tol, stringsAsFactors = FALSE)
}

neg_status <- "not run"
compared_ok <- which(vapply(records, function(r) r$Structural_Problems == 0 && r$Values_Compared > 0, logical(1)))
if (NEGATIVE_CONTROL && length(compared_ok) == 0) {
  say("")
  say("NEGATIVE CONTROL not run: no pair could be compared (see the missing files above).")
}
if (NEGATIVE_CONTROL && length(compared_ok) > 0) {
  say("")
  say("NEGATIVE CONTROL (one value in a copy of the NEW table of pair '", PAIRS[[compared_ok[1]]]$label, "' is changed by 0.001)")
  ctl <- compare_pair(PAIRS[[compared_ok[1]]], perturb = TRUE)
  detected <- ctl$stats$n_unexpected + ctl$stats$n_expected > 0
  neg_status <- if (detected) "detected" else "NOT DETECTED"
  say("  Change ", if (detected) "was detected: the checks can fail, so the PASS results mean something."
      else "was NOT detected: do not trust these results (e.g. every column is ignored or declared expected).")
}

overall_pass <- all(vapply(records, function(r) r$Status == "PASS", logical(1))) && !(neg_status %in% c("NOT DETECTED", "not run"))
say("")
say(strrep("-", 72))
say(if (overall_pass) "RESULT: PASS -- no unexpected differences; negative control behaved as it should."
    else "RESULT: FAIL -- see lines starting with '!!' (or the negative-control message).")

stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
report_path <- file.path(REPORT_FOLDER, paste0("comparison_report_", stamp, ".txt"))
record_path <- file.path(REPORT_FOLDER, paste0("golden_record_", stamp, ".csv"))
writeLines(report, report_path)
rec_df <- do.call(rbind, records)
rec_df$Negative_Control <- neg_status
rec_df$Checked_At <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
write.csv(rec_df, record_path, row.names = FALSE)
cat("\nSaved:\n  ", report_path, "\n  ", record_path, "\n", sep = "")
