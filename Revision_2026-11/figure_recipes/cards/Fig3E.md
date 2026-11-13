# Fig3E: how it was made

Run record(s): 2026-11-02_08-48-33_Linegraphs_ByCluster_p_primary_Fig3E   [file label: 3E]

## Run folder `2026-11-02_08-48-33_Linegraphs_ByCluster_p_primary_Fig3E   [file label: 3E]`

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
- `THRESHOLD_CATEGORY_FOLDER`: 
- `CATEGORIES_TO_ANALYSE`: []
- `LOG2FC_INF_AS_NA`: true
- `FILL_EMPTY_CLUSTERS`: true
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: 5 * sinh(x) on every Plot = 1 channel
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv` e5e5ee9d; `Expt1234_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` 9ca0d4d7; `Expt1234_PercAboveThresh_95thPerc_FC.csv` a5330fbb; `Expt1234_PercAboveThresh_95thPerc_FC_RepsAvg.csv` 33b974e6; `Expt1234_PercAboveThresh_95thPerc_log2FC.csv` 2732dcb2; `Expt1234_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` f05c5288; `Thresholds.csv` 01e36724

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
- `CLUSTERS`: [1, 2, 4, 7, 9, 13, 17, 18]
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 900
- `Table_Rows_Used`: 152
- `Reference_Rows`: 32
- `Comparisons`: 1056
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 7b8d9050; `Statistics_pvalues_wide.csv` a7b94f99; `Statistics_design.csv` 1318645d; `Statistics_litters.csv` 02533e0c; `Statistics_diagnostics.csv` f4137396

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: ["1", "2", "4", "7", "9", "13", "17", "18"]
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Statistics_Folder`: 2026-10-29_16-59-57_Statistics_FC_welch_byCluster_Cl1-2-4-7-9-13-17-18_Fig3_BDNF
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 176
- `Comparisons`: 880
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 164
- `Significant_Primary`: 90
- `Significant_BH_All_Markers`: 44
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 28e218f9; `Timecourse_families.csv` 07ab6d16; `Timecourse_pvalues_wide.csv` b224f919

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
- `CLUSTERS`: ["9", "17", "18"]
- `CATEGORY`: none
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 5
- `PANEL_SIZE_IN`: 5
- `CAPTION_ON_FIGURE`: false
- `BACKGROUND`: white
- `TRANSPARENT_ON`: white
- `Palette_File`: Palette_Colors.csv
- `POINT_SHAPES`: auto
- `Statistics_Folder`: 2026-10-29_17-00-07_Timecourse_FC_welch-welch_holm_Fig3_BDNF
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 28e218f95514b6eadade93ab88384c77
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Panels`: 5
- `Starred_Group_TimePoints`: 28
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` 09fa968b; `Linegraph_Clusters_ByCluster-Grid.png` 293aaab9; `Plotted_data.csv` 0170cd8d; `Plotted_points.csv` 4d94298c; `Plotted_axes.csv` 5f172f4b; `Plotted_style.csv` a49181cd; `Figure_legend.txt` f1e16f1c
