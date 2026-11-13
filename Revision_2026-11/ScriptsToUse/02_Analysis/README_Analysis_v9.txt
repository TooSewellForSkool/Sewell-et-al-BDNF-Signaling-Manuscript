================================================================================
README: Analysis layer          [01_PercAboveThreshold.R (v3), 02_AssignThresholdCategories.R, 03_Statistics.R (v2), 04_TimecourseStatistics.R (v1)]
================================================================================
Status: 01 v2 and 02 passed the golden check on 2026-10-23 (13 runs, 67 of 67 pairs identical; comparison_report_2026-10-23_14-53-47).
        01 v3 (2026-10-25) = v2 + FILL_EMPTY_CLUSTERS (default FALSE); golden run E below.
        03_Statistics v1 (2026-10-25): published method only (Welch t-tests against NoTreat 0 hr, no adjustment). Golden runs P1-P4
        (Tables S3, S4A-C) all passed (comparison_report_2026-10-25_22-19-17; 78 of 78 pairs; 03 v2 rechecked 23:20, 78 of 78; P3 at tolerance 5e-4 because two published cells
        carry 3 significant figures).
        03_Statistics v2 (2026-10-25) = v1 + design summary with litters + assumption diagnostics; the p-values are unchanged.
        Needs Litter_Map.csv (next to this README). Check it against your records.

WHAT THIS LAYER DOES
--------------------
Turns the cleaned, clustered data (DataCleanup + LeidenClustering) into the TABLES
behind the figures. It makes no plots; plots belong to the Visualization layer and
read these tables only.

SCRIPT 1: 01_PercAboveThreshold.R
---------------------------------
Replaces the table-making code of 10.1, 12.1 and the PseudoBulk/Clusters mode of 12.1.1.
For every sample and cluster it calculates, for every marker:
  1. a threshold = the 95th percentile of the baseline samples (TimePoint 0, e.g.
     NoTreat 0 hr). The percentile is taken for each baseline file, then averaged.
  2. the percent of cells above that threshold (strictly greater than)
  3. fold change and log2 fold change against the control condition(s)
  ANALYSIS_BY = "PseudoBulk": one threshold from all cells (Cluster = 0)
  ANALYSIS_BY = "Clusters"  : one threshold per cluster

BEFORE YOU START
  * Input: a folder with expression_matrix_analysis.csv, expression_matrix_other.csv,
    filenums.csv, metadata.csv, panel.csv and cluster_R<N>_assigns.csv (or
    ..._assigns_merged.csv). Folders made by the DataCleanup scripts are checked
    against the run record; the published (older) folders have none and are
    accepted with a note.
  * Packages: svDialogs, rstudioapi, data.table, jsonlite.
  * metadata.csv: the published layout (Condition = "BDNF_15min", Treatment = "BDNF",
    columns Expt and Reordered) and the new layout (01_Get_Metadata v8: Treatment =
    "BDNF_15min", Condition = "BDNF", Expt_ID, no Reordered) are both understood. The
    mapping used is written to the run record. Without a Reordered column the plotting
    order is derived (baseline first, then treatment, time, replicate).

HOW TO RUN
  1. Open 01_PercAboveThreshold.R, edit USER SETTINGS (each setting is explained
     beside it), and click Source.
  2. A new folder <date>_<time>_PercAboveThresh_<mode>_<samples> appears inside the
     input folder.

WHAT YOU GET (in the new folder)
  Expt<set>_PercAboveThresh_95thPerc_RAWvalues.csv         percent above threshold per sample (and cluster)
  Expt<set>_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv same, averaged over replicates (N_Replicates given)
  ..._FC.csv, ..._FC_RepsAvg.csv                           fold change vs control
  ..._log2FC.csv, ..._log2FC_RepsAvg.csv                   log2 fold change
  Thresholds.csv                                           threshold per marker and cluster, with n baseline files/cells
  run_parameters.json                                      settings, files, cells kept, markers, warnings, checksums

  Columns of every table: File_Num, File_Name, Sample_Condition, Reordered, TimePoint,
  Replicate, Cluster, [N_Cells in RAWvalues], then markers.

SCRIPT 2: 02_AssignThresholdCategories.R   (replaces the data half of 13.1_ExtractByMarkerExpression_HiLo_v5)
---------------------------------------------------------------------------------------------------------------
Sorts every cell into a category for one or two markers by cutting the marker at the values
you choose (e.g. TrkB_i at 0.15, 1.4, 3.7 -> levels 1-4; two markers -> "14" = p75NTR level 1,
TrkB_i level 4). It writes ThresholdCategory.csv (one row per cell), the cut values, the percent of
cells per category and the cells per sample per category. Nothing from the input folder is copied
(13.1 copied everything into every output folder). Plots go to the Visualization layer.
To analyse within categories: run 02, then in 01 set THRESHOLD_CATEGORY_FOLDER to the folder 02 made
(and, if you want only some categories, CATEGORIES_TO_ANALYSE, e.g. c("11","41","14","44")).
With categories, 01 does everything inside each category: threshold from that category's baseline
cells, percent of that category's cells, fold change against the same category's baseline.
01 v2 with THRESHOLD_CATEGORY_FOLDER = "" gives the same tables as v1.

SCRIPT 3: 03_Statistics.R   (replaces the statistics code inside the LineGraphAnalysis_... plotting scripts)
--------------------------------------------------------------------------------------------------------
Compares every group at every time point with a reference group (NoTreat at 0 hr), marker by marker, and writes
the p-values as tables that the plotting scripts read. v1 reproduces the published method exactly:
  t.test(group, reference, paired = FALSE) in R's default form = WELCH's t-test (the legends say "Student's"),
  unpaired, no multiple-testing adjustment. Samples are treated as independent (replicate numbers are not litters).
Input: a folder made by 01_PercAboveThreshold.R (it needs that script's run record). Choose STAT_TABLE
("RAWvalues", "FC" or "log2FC"; the *_RepsAvg tables are never tested). To compare clusters or categories instead of
treatments set COMPARE_BY; the script stops if the selection would mix clusters, categories or treatments inside
one compared group. A new folder <date>_<time>_Statistics_... appears INSIDE the 01 folder with:
  Statistics_long.csv          one row per marker x group x time point: n of group and reference, means, test, p, stars
  Statistics_pvalues_wide.csv  the layout of the Supplementary Tables (rows = time point x group, columns = markers)
  Statistics_design.csv        samples (and replicate numbers) behind every group and time point
  Statistics_litters.csv       (v2) the litter of every sample used, from Litter_Map.csv
  Statistics_diagnostics.csv   (v2) Shapiro-Wilk on residuals, Brown-Forsythe, variance ratio, skewness, mean-SD trend per marker and
                               family (one treatment over time + the reference); also on log2 scale when STAT_TABLE = "FC".
                               Reported only: nothing is skipped or chosen from them. With n = 2-4 per group they have little power.
  run_parameters.json          settings, input table checksum, warnings, output checksums
Read the change log at the top of the script (A1-A6, C1-C3) before using it for anything new.

SCRIPT 4: 04_TimecourseStatistics.R   (NEW for the revision; not part of the published analysis)
--------------------------------------------------------------------------------------------------
Answers the reviewer's comment (use ANOVA across time points; adjust for multiple testing). Input: a folder made by
03_Statistics.R v2 (it re-reads the same table and selection that 03 recorded). For every marker and every "family"
(one treatment, one cluster/category: its time points and the NoTreat 0 hr reference) it writes:
  Timecourse_long.csv         Omnibus_p, p_unadjusted (pairwise Welch, as published), p_primary (adjusted within
                              treatment x marker: Dunnett, or Holm on Welch), q_BH_all_markers (alternative family:
                              BH over all markers x time points of the treatment), stars on p_primary
  Timecourse_families.csv     the suggestion (classic ANOVA + Dunnett vs Welch ANOVA + Welch/Holm) and why, the checks
                              behind it, the omnibus p and its BH-adjusted version across markers, Suggested_Scale
  Timecourse_pvalues_wide.csv p_primary in the layout the plotting scripts read
The suggestion follows the design first (independent groups, one factor, many-to-one comparison) and then the spread
checks; normality is reported but never changes the test. OMNIBUS and POSTHOC override it. Samples are treated as
independent; litters are listed (Statistics_litters.csv) but not modelled (a mixed model is a possible later upgrade).
Dunnett uses random numbers (package multcomp): SEED is set and recorded.

GOLDEN CHECKS AGAINST THE PUBLISHED RESULTS
  Tools/Compare_GoldenFiles.R (v3.2) finds every published file inside GOLDEN_ROOT by its manifest
  path and checks its MD5, and checks every new run folder's settings. You only fill in GOLDEN_ROOT
  and the folders in NEW_RUNS (leave a run blank to skip it). Send me comparison_report_*.txt and
  golden_record_*.csv. Settings not listed below stay at the script's default:
  INPUT_FOLDER = the published input folder, CLUSTER_ROUND 1, SAMPLE_SET c(1,2,3,4), THRESHOLD_PROB 0.95,
  THRESHOLD_TIMEPOINT 0, CTRL_PER_TIMEPOINT FALSE, CLUSTERS_TO_REMOVE c().

  DONE 2026-10-23 with 01 v1:
   Run A  01: ANALYSIS_BY "PseudoBulk", CONDITIONS_TO_EXCLUDE c()            vs 10.1 45-sample tables (7)
   Run B  01: ANALYSIS_BY "PseudoBulk", CONDITIONS_TO_EXCLUDE c("CompMed")   vs Fig 3C tables (7)
   Run C  01: ANALYSIS_BY "Clusters", CLUSTERS_MERGED TRUE                   vs 12.1 tables (7)
  TO DO (01 v2 and 02). Please also re-run A, B, C with v2 (same settings): they must still pass.
   Run D95  01 v2: as A, plus LOG2FC_INF_AS_NA TRUE                          vs 12.1.1 Fig 3D tables (7)
   Run D75  01 v2: as D95, THRESHOLD_PROB 0.75                               vs 12.1.1 SuppFig 3 (75th) tables (7)
   Run D50  01 v2: as D95, THRESHOLD_PROB 0.5                                vs 12.1.1 SuppFig 3 (50th) tables (7)
   Run T1   02: SELECTED_MARKERS c(x = "TrkB_i_174Yb"), CUT_VALUES list(x = c(0.15, 1.4, 3.7))
            (CUTS_JUSTIFIED TRUE)                                             vs published ThresholdCategory.csv (Fig 5 TrkB_i folder)
   Run T2   02: c(x = "p75NTR_161Dy"), list(x = c(0.45, 1.45, 3.95))         vs published ThresholdCategory.csv (Fig 5 p75NTR folder)
   Run T3   02: c(x = "p75NTR_161Dy", y = "TrkB_i_174Yb"),
            list(x = c(0.45, 1.45, 3.95), y = c(0.15, 1.4, 3.7))             vs published ThresholdCategory.csv (Fig 6A folder)
   Run T4   02: c(x = "p75NTR_161Dy", y = "TrkB_i_174Yb"),
            list(x = c(1.45), y = c(0.15, 1.4, 3.7))                          vs published ThresholdCategory.csv (SuppFig 12 folder)
   Run K1   01 v2: THRESHOLD_CATEGORY_FOLDER = the T1 folder, CATEGORIES_TO_ANALYSE c(),
            LOG2FC_INF_AS_NA TRUE                                             vs Fig 5 TrkB_i tables, categories 1-4 (7)
   Run K2   01 v2: THRESHOLD_CATEGORY_FOLDER = the T3 folder, CATEGORIES_TO_ANALYSE c("11","41","14","44"),
            LOG2FC_INF_AS_NA TRUE                                             vs Fig 6A tables (7)
   Run K3   01 v2: THRESHOLD_CATEGORY_FOLDER = the T2 folder, CATEGORIES_TO_ANALYSE c(),
            LOG2FC_INF_AS_NA TRUE                                             vs p75NTR categories 1-4 tables (7)
  Run 02 before the K runs (they need its folder). All runs use the same input folder.
  (All of the above were run on 2026-10-23: 67/67 pairs identical.)

  03_Statistics v1, Table S3 (run P1): run on 2026-10-25. All 220 published p-values matched (largest difference
   4.5e-8, i.e. the rounding of the published table). The first report said FAIL only because the new table also tests
   the 20 markers that Table S3 does not show; Compare_GoldenFiles v3.3 declares those extra markers (every published
   marker must still match). Re-run the comparison with v3.3 and the SAME P1 folder; no new run is needed.
   Run P1   03: INPUT_FOLDER = the Run A folder, STAT_TABLE "RAWvalues", TEST "welch", COMPARE_BY "Treatment",
            REFERENCE_TREATMENT "NoTreat", REFERENCE_TIMEPOINT 0, REFERENCE_SCOPE "pooled",
            TREATMENTS c("NoTreat", "ContDepr", "BDNF"), CLUSTERS c(), CATEGORIES c(), MARKERS c()
            -> vs Tools/GoldenData/PublishedStats_TableS3_long.csv (tolerance 1e-7, published values are rounded)

  TO DO (01 v3, run E): the 20-cluster tables of the Fig 3 folder (round 1, NOT merged; 45 samples x 20 clusters = 900 rows)
   Run E    01 v3: ANALYSIS_BY "Clusters", CLUSTER_ROUND 1, CLUSTERS_MERGED FALSE, CONDITIONS_TO_EXCLUDE c(),
            THRESHOLD_PROB 0.95, LOG2FC_INF_AS_NA TRUE, FILL_EMPTY_CLUSTERS TRUE,
            THRESHOLD_CATEGORY_FOLDER ""                                      vs 12.1.1 Fig 3 cluster tables (7 pairs)
   (12.1.1 added a row for every sample x cluster, with NA where a sample has no cells in a cluster - clusters 19 and 20 are
    very small; FILL_EMPTY_CLUSTERS TRUE does the same. Run C is the 4 merged clusters of Supplementary Fig 6.)

  TO DO (03_Statistics, Supplementary Table 4; all three use the Run E folder as INPUT_FOLDER):
   STAT_TABLE "FC", TEST "welch", COMPARE_BY "Cluster", REFERENCE_TREATMENT "NoTreat", REFERENCE_TIMEPOINT 0,
   REFERENCE_SCOPE "pooled", CATEGORIES c(), MARKERS c()  (the published p = 1 at time point 0 for every cluster shows
   the tests were run on the fold-change table)
   Run P2   TREATMENTS c("NoTreat", "BDNF"),    CLUSTERS c(1, 2, 4, 7, 13)    -> Table S4A (Fig 3C)
   Run P3   TREATMENTS c("NoTreat", "CompMed"), CLUSTERS c(1, 3, 15, 17)      -> Table S4B (Fig 3D)
   Run P4   TREATMENTS c("NoTreat", "BDNF"),    CLUSTERS c(9, 17, 18)         -> Table S4C (Fig 3E)
   (the treatment for each panel is taken from the published folder names; tell me if a panel used another)
  In Compare_GoldenFiles.R (v3.3) set GOLDEN_STATS_FOLDER, and E, P2, P3, P4 (the folder 03 made = the Statistics_... folder
  INSIDE the Run E folder). Runs A-K may stay as they are. Not yet covered: the Fig 6 / Supp Fig 12 p-value tables.

IMPORTANT NOTES
  * Read the change log at the top of the script (items A1-A8) before relying on it.
  * The threshold is calculated from baseline files even if they are listed in
    FILES_TO_EXCLUDE (as published); see change log A3.
  * Samples and results are joined by File_Num, not by row position (change log A2).
  * No random numbers are used.
  * Threshold categories are in 01 v2 + 02 (see above). Not reproduced on purpose: the mixed mode of
    12.1 (threshold from category cells, percentages from all cells); no published figure used it.
  * 12.1.1 turns -Inf log2 fold changes into NA before averaging replicates; 10.1/12.1 tables keep -Inf.
    LOG2FC_INF_AS_NA chooses which (change log A10). No table checked so far contained -Inf.

COMMON PROBLEMS
  "No baseline cells found"       No analysed sample has TimePoint = THRESHOLD_TIMEPOINT.
  "...cluster file belongs to a different dataset"   Row counts differ from filenums.csv.
  "CTRL_PER_TIMEPOINT is TRUE, so NORM_CONDITIONS must list..."   List the controls.
  Warning "no 'Reordered' column"  Plotting order was derived; add Reordered to metadata.csv to fix it.


HOW TO RUN THE v2 / v1 CHECKS (2026-10-25)
  Golden check of 03 v2: in Tools/Compare_GoldenFiles.R set Run P1 (and, if you re-run them, P2-P4) to the NEW 03_Statistics v2
  output folder (the <date>_Statistics_... folder inside the Run A folder). The 22:52 comparison still pointed at the v1 folder.
  04_TimecourseStatistics v1 has no published results to match. It is checked against an independent Python calculation
  (Inventory/Reference04_python_timecourse.py): send me Timecourse_long.csv and Timecourse_families.csv.
  Settings for that check: 03 with the P1 settings (Run A folder, STAT_TABLE "RAWvalues", RUN_DIAGNOSTICS TRUE), then 04 with INPUT_FOLDER = that 03 folder, all else default.


MEAN EXPRESSION MODE (01_PercAboveThreshold.R version 4, 2026-10-27)
--------------------------------------------------------------------
MEASURE <- "mean_expression" makes the same tables from the MEAN of each marker per sample and cluster (the quantity of the
published Fig 5F, SelectPlots_v9) instead of the percent of cells above a threshold. Fold change is against the baseline mean of the
same cluster. MEAN_SCALE "as_stored" (default) averages the values as stored (arcsinh scale), as Fig 5F did. Output names:
Expt<set>_MeanExpr_FC.csv etc. (no Thresholds.csv). The default MEASURE "percent_above" gives exactly the version 3 tables.
Check: Inventory/Fig5F_points_extracted_from_PDFs.csv (replicate values recovered from the published Fig 5F PDFs).

SCRIPT 05: 05_ConditionStatistics.R   (NEW for the revision, 2026-10-26; replaces the t-tests drawn by 16_Barplot_v7.r)
---------------------------------------------------------------------------------------------------------------------
Compares CONDITIONS at one time point (ContDepr, BDNFTest, BDNF+K252aTest) for every marker in every cluster/category.
Input:  a folder holding the fold-change table (*_FC.csv, not RepsAvg) from 01_PercAboveThreshold.R (or the merged-category table).
Per marker x cluster/category (version 2): Welch's ANOVA (omnibus), then ALL THREE pairs with Welch t-tests + Holm, the SAME method for every marker
(decision 2026-10-26; Tukey HSD and the classic ANOVA are reported in p_Tukey / Omnibus_p_classic for the supplementary tables; the spread-based suggestion stays in Conditions_families.csv),
then ALL THREE pairs: ContDepr-BDNFTest, ContDepr-BDNF+K252aTest, BDNFTest-BDNF+K252aTest.
Columns: p_unadjusted (Welch, = the published p), p_Holm, p_primary (= p_Holm by default), p_Tukey, q_BH_all_markers
(BH over all markers x pairs of the cluster/category; family size is written in Family_Size_BH), Omnibus_p, q_omnibus_BH_markers.
Output folder <time>_Conditions_...: Conditions_long.csv, Conditions_families.csv (diagnostics + suggestion), Conditions_groupstats.csv
(n, mean, SD), Conditions_input_values.csv (the samples), run_parameters.json. No random numbers.
Check: Inventory/Reference05_python_conditions.py (same rules in Python) and its outputs Reference05v2_pb.csv / Reference05v2_merged3cat.csv / Reference05v2_clusters20.csv. Golden: the published bar-graph p-values (Welch, unadjusted) reproduce in all 41 published panels (2026-10-26).
Note: with n = 3 per group the variance ratio is > 3 for ~85-90% of families, which is why one method (Welch) is used everywhere.

TIP: KEEP THE PROJECT FOLDER PATH SHORT
---------------------------------------
Each step writes its output folder INSIDE the folder of the step before it, so the full path of a late figure can
become very long (about 500 characters in the Fig 6 runs). Very long paths can fail: Windows stops at 260
characters, and R's pdf() device cuts paths over 511 characters (the figure scripts now work around the latter by
writing the PDF to a short temporary path, but other programs, zip tools and Illustrator may still struggle).
Work in a short project folder, e.g. C:\cytof or ~/cytof, and avoid long folder names above it. If a run says a file
"was not found after saving", the path is the first thing to check.
