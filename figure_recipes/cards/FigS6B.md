# FigS6B: how it was made

Run record(s): 2026-11-03_01-02-43_Linegraphs_ByCluster_p_primary_FigS6   [file label: S6B]

## Run folder `2026-11-03_01-02-43_Linegraphs_ByCluster_p_primary_FigS6   [file label: S6B]`

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
- `CLUSTERS_MERGED`: true
- `Cluster_File`: cluster_R1_assigns_merged.csv
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
- `FILL_EMPTY_CLUSTERS`: false
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: 5 * sinh(x) on every Plot = 1 channel
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv` 122094de; `Expt1234_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` 5ec1fa1a; `Expt1234_PercAboveThresh_95thPerc_FC.csv` 1c688465; `Expt1234_PercAboveThresh_95thPerc_FC_RepsAvg.csv` 20c2341b; `Expt1234_PercAboveThresh_95thPerc_log2FC.csv` 2938ce3c; `Expt1234_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` b5a7898f; `Thresholds.csv` 6e024968

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
- `CLUSTERS`: [2, 3]
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 1c688465d089b243e17a2733118591df
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 180
- `Table_Rows_Used`: 38
- `Reference_Rows`: 8
- `Comparisons`: 264
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` d7d7861a; `Statistics_pvalues_wide.csv` 8692773d; `Statistics_design.csv` d6657e14; `Statistics_litters.csv` 02533e0c; `Statistics_diagnostics.csv` 0328905f

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: ["2", "3"]
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 1c688465d089b243e17a2733118591df
- `Statistics_Folder`: 2026-11-03_00-41-43_Statistics_FC_welch_byCluster_Cl2-3_FigS6
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 44
- `Comparisons`: 220
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 61
- `Significant_Primary`: 33
- `Significant_BH_All_Markers`: 17
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 50bc38c3; `Timecourse_families.csv` aa36f6c5; `Timecourse_pvalues_wide.csv` 78cad05c

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
- `CLUSTERS`: ["2", "3"]
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
- `Statistics_Folder`: 2026-11-03_00-43-19_Timecourse_FC_welch-welch_holm_FigS6
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 50bc38c3ff6f43259b1ee04fcc72389a
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 1c688465d089b243e17a2733118591df
- `Panels`: 5
- `Starred_Group_TimePoints`: 24
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` 67708b3e; `Linegraph_Clusters_ByCluster-Grid.png` acad6723; `Plotted_data.csv` 310f0ea8; `Plotted_points.csv` 8b5ceec5; `Plotted_axes.csv` 7d142406; `Plotted_style.csv` c61bd4b8; `Figure_legend.txt` d1a98e8c
