# Fig3D: how it was made

Run record(s): 2026-11-02_08-49-57_Linegraphs_ByCluster_p_primary_Fig3D   [file label: 3D]

## Run folder `2026-11-02_08-49-57_Linegraphs_ByCluster_p_primary_Fig3D   [file label: 3D]`

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
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: [1, 3, 15, 17]
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 900
- `Table_Rows_Used`: 68
- `Reference_Rows`: 16
- `Comparisons`: 528
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 7c8ef6e0; `Statistics_pvalues_wide.csv` 1067ad7f; `Statistics_design.csv` 9655eadc; `Statistics_litters.csv` b8c8b57c; `Statistics_diagnostics.csv` 71e3786a

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: ["1", "3", "15", "17"]
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Statistics_Folder`: 2026-10-29_17-00-50_Statistics_FC_welch_byCluster_Cl1-3-15-17_Fig3_CompMed
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 88
- `Comparisons`: 440
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 89
- `Significant_Primary`: 38
- `Significant_BH_All_Markers`: 31
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 7c0111cd; `Timecourse_families.csv` f9f19790; `Timecourse_pvalues_wide.csv` 5e67e1d6

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Cluster
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: ["1", "3", "15", "17"]
- `CATEGORY`: none
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 1
- `PANEL_SIZE_IN`: 5
- `CAPTION_ON_FIGURE`: false
- `BACKGROUND`: white
- `TRANSPARENT_ON`: white
- `Palette_File`: Palette_Colors.csv
- `POINT_SHAPES`: auto
- `Statistics_Folder`: 2026-10-29_17-00-55_Timecourse_FC_welch-welch_holm_Fig3_CompMed
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 7c0111cd1638704bc4d78a0d1c41bf54
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Panels`: 4
- `Starred_Group_TimePoints`: 16
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` 7ad16ea6; `Linegraph_Clusters_ByCluster-Grid.png` 18e19bca; `Plotted_data.csv` 7b8b14e5; `Plotted_points.csv` 9f69f0d5; `Plotted_axes.csv` 2cc6fc49; `Plotted_style.csv` 5831fe00; `Figure_legend.txt` 01b2f067
