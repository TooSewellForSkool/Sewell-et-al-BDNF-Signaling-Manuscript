# Fig6D: how it was made

Run record(s): 2026-11-03_16-00-35_Linegraphs_ByCluster_Cat23_p_primary_Fig6D   [file label: 6D]

## Run folder `2026-11-03_16-00-35_Linegraphs_ByCluster_Cat23_p_primary_Fig6D   [file label: 6D]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `03_Statistics` | 2 |
| 3 | `04_TimecourseStatistics` | 1.2 |
| 4 | `01_PlotLinegraphs` | 1.12 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: percent_above
- `MEAN_SCALE`: null
- `ANALYSIS_BY`: Clusters
- `CLUSTER_ROUND`: 1
- `CLUSTERS_MERGED`: false
- `Cluster_File`: cluster_R1_assigns.csv
- `SAMPLE_SET_COLUMN`: Expt
- `SAMPLE_SET`: [1, 2, 3, 4]
- `CONDITIONS_TO_EXCLUDE`: []
- `CLUSTERS_TO_REMOVE`: []
- `FILES_TO_EXCLUDE`: []
- `THRESHOLD_PROB`: 0.95
- `THRESHOLD_TIMEPOINT`: 0
- `CTRL_PER_TIMEPOINT`: false
- `NORM_CONDITIONS_Used`: ["NoTreat_0hr"]
- `THRESHOLD_CATEGORY_FOLDER`: BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/2026-10-29_15-50-33_Thr...
- `CATEGORIES_TO_ANALYSE`: []
- `LOG2FC_INF_AS_NA`: true
- `FILL_EMPTY_CLUSTERS`: false
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: 5 * sinh(x) on every Plot = 1 channel
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv` 6fc86bf0; `Expt1234_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` 81bc4f0e; `Expt1234_PercAboveThresh_95thPerc_FC.csv` 476736c6; `Expt1234_PercAboveThresh_95thPerc_FC_RepsAvg.csv` 4fb6c980; `Expt1234_PercAboveThresh_95thPerc_log2FC.csv` 6c57d014; `Expt1234_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` fc39fc15; `Thresholds.csv` 437c2c63

### 03_Statistics (v2)

- `STAT_TABLE`: FC
- `TEST`: welch
- `COMPARE_BY`: Cluster
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `REFERENCE_SCOPE`: pooled
- `RUN_DIAGNOSTICS`: true
- `LITTER_MAP`: 
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: [1, 2, 4, 18]
- `CATEGORIES`: ["23"]
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 476736c6911cdd26bf82955cf79afd0e
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 4748
- `Table_Rows_Used`: 76
- `Reference_Rows`: 16
- `Comparisons`: 528
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 402de454; `Statistics_pvalues_wide.csv` 523a8279; `Statistics_design.csv` 8580f4f8; `Statistics_litters.csv` 02533e0c; `Statistics_diagnostics.csv` 21c17aa4

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: ["1", "2", "4", "18"]
- `CATEGORIES`: ["23"]
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 476736c6911cdd26bf82955cf79afd0e
- `Statistics_Folder`: 2026-10-29_17-07-02_Statistics_FC_welch_byCluster_Cl1-2-4-18_Cat23_Fig6D
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 88
- `Comparisons`: 440
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 89
- `Significant_Primary`: 37
- `Significant_BH_All_Markers`: 16
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 0e560f2b; `Timecourse_families.csv` 35adfdcc; `Timecourse_pvalues_wide.csv` a2d16b2f

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Cluster
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: ["1", "2", "4", "18"]
- `CATEGORY`: 23
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 5
- `PANEL_SIZE_IN`: 5
- `CAPTION_ON_FIGURE`: false
- `BACKGROUND`: white
- `TRANSPARENT_ON`: white
- `Palette_File`: Palette_Colors.csv
- `POINT_SHAPES`: auto
- `Statistics_Folder`: 2026-10-29_17-07-22_Timecourse_FC_welch-welch_holm_Fig6D
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 0e560f2b1911e8e60962dba3ff54f493
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 476736c6911cdd26bf82955cf79afd0e
- `Panels`: 5
- `Starred_Group_TimePoints`: 23
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid-Cat23.pdf` a12fbf57; `Linegraph_Clusters_ByCluster-Grid-Cat23.png` cf88c9b4; `Plotted_data.csv` fb2eaa16; `Plotted_points.csv` 8c6d5bf3; `Plotted_axes.csv` b258739a; `Plotted_style.csv` 8553f282; `Figure_legend.txt` ddea793a
