# Fig5F: how it was made

Run record(s): 2026-11-01_14-13-44_Linegraphs_ByCluster_p_primary_Fig5F_CompMed-Top   [file label: 5F_CompMed-Top], 2026-11-01_14-14-04_Linegraphs_ByCluster_p_primary_Fig5F_CompMed-Bot   [file label: 5F_CompMed-Bot], 2026-11-01_14-15-24_Linegraphs_ByCluster_p_primary_Fig5F_BDNF-Top   [file label: 5F_BDNF-Top], 2026-11-01_14-15-49_Linegraphs_ByCluster_p_primary_Fig5F_BDNF-Bot   [file label: 5F_BDNF-Bot], 2026-11-01_14-17-22_Linegraphs_ByCluster_p_primary_Fig5F_ContDepr-Top   [file label: 5F_ContDepr-Top], 2026-11-01_14-17-44_Linegraphs_ByCluster_p_primary_Fig5F_ContDepr-Bot   [file label: 5F_ContDepr-Bot]

## Run folder `2026-11-01_14-13-44_Linegraphs_ByCluster_p_primary_Fig5F_CompMed-Top   [file label: 5F_CompMed-Top]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `03_Statistics` | 2 |
| 3 | `04_TimecourseStatistics` | 1.2 |
| 4 | `01_PlotLinegraphs` | 1.12 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: mean_expression
- `MEAN_SCALE`: as_stored
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
- `FILL_EMPTY_CLUSTERS`: false
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: none (mean of the values as stored, as Fig 5F)
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_MeanExpr_RAWvalues.csv` 1975beca; `Expt1234_MeanExpr_RAWvalues_RepsAvg.csv` 481b4d18; `Expt1234_MeanExpr_FC.csv` 950a5b7f; `Expt1234_MeanExpr_FC_RepsAvg.csv` c67a47b1; `Expt1234_MeanExpr_log2FC.csv` 6d75945c; `Expt1234_MeanExpr_log2FC_RepsAvg.csv` de4616f0

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
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 891
- `Table_Rows_Used`: 338
- `Reference_Rows`: 79
- `Comparisons`: 120
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` fb904e53; `Statistics_pvalues_wide.csv` de5c1fb7; `Statistics_design.csv` d03fb3d4; `Statistics_litters.csv` b8c8b57c; `Statistics_diagnostics.csv` bf9c7c99

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Statistics_Folder`: 2026-11-01_14-08-30_Statistics_FC_welch_byCluster_Fig5F_CompMed
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 20
- `Comparisons`: 100
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 12
- `Significant_Primary`: 3
- `Significant_BH_All_Markers`: 3
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 99ce288a; `Timecourse_families.csv` a8a140d4; `Timecourse_pvalues_wide.csv` c8bbf22d

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Cluster
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: [0.5, 1.5]
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: ["1", "15"]
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
- `Statistics_Folder`: 2026-11-01_14-08-42_Timecourse_FC_welch-welch_holm_Fig5F_CompMed
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 99ce288ae6b8a302a199fe5d9bc2044d
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Panels`: 1
- `Starred_Group_TimePoints`: 0
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` 0fd613d9; `Linegraph_Clusters_ByCluster-Grid.png` 2f34af1e; `Plotted_data.csv` f5d36f48; `Plotted_points.csv` aec22de1; `Plotted_axes.csv` bd16831c; `Plotted_style.csv` 876eec3f; `Figure_legend.txt` d74ab880

## Run folder `2026-11-01_14-14-04_Linegraphs_ByCluster_p_primary_Fig5F_CompMed-Bot   [file label: 5F_CompMed-Bot]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `03_Statistics` | 2 |
| 3 | `04_TimecourseStatistics` | 1.2 |
| 4 | `01_PlotLinegraphs` | 1.12 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: mean_expression
- `MEAN_SCALE`: as_stored
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
- `FILL_EMPTY_CLUSTERS`: false
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: none (mean of the values as stored, as Fig 5F)
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_MeanExpr_RAWvalues.csv` 1975beca; `Expt1234_MeanExpr_RAWvalues_RepsAvg.csv` 481b4d18; `Expt1234_MeanExpr_FC.csv` 950a5b7f; `Expt1234_MeanExpr_FC_RepsAvg.csv` c67a47b1; `Expt1234_MeanExpr_log2FC.csv` 6d75945c; `Expt1234_MeanExpr_log2FC_RepsAvg.csv` de4616f0

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
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 891
- `Table_Rows_Used`: 338
- `Reference_Rows`: 79
- `Comparisons`: 120
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` fb904e53; `Statistics_pvalues_wide.csv` de5c1fb7; `Statistics_design.csv` d03fb3d4; `Statistics_litters.csv` b8c8b57c; `Statistics_diagnostics.csv` bf9c7c99

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Statistics_Folder`: 2026-11-01_14-08-30_Statistics_FC_welch_byCluster_Fig5F_CompMed
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 20
- `Comparisons`: 100
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 12
- `Significant_Primary`: 3
- `Significant_BH_All_Markers`: 3
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 99ce288a; `Timecourse_families.csv` a8a140d4; `Timecourse_pvalues_wide.csv` c8bbf22d

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Cluster
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: [0.5, 1.5]
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: ["2", "13", "4"]
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
- `Statistics_Folder`: 2026-11-01_14-08-42_Timecourse_FC_welch-welch_holm_Fig5F_CompMed
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 99ce288ae6b8a302a199fe5d9bc2044d
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Panels`: 1
- `Starred_Group_TimePoints`: 1
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` 0798ba24; `Linegraph_Clusters_ByCluster-Grid.png` dd90e1ec; `Plotted_data.csv` 02a080bc; `Plotted_points.csv` 166b82d0; `Plotted_axes.csv` bd16831c; `Plotted_style.csv` eae98fa1; `Figure_legend.txt` f12ecaee

## Run folder `2026-11-01_14-15-24_Linegraphs_ByCluster_p_primary_Fig5F_BDNF-Top   [file label: 5F_BDNF-Top]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `03_Statistics` | 2 |
| 3 | `04_TimecourseStatistics` | 1.2 |
| 4 | `01_PlotLinegraphs` | 1.12 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: mean_expression
- `MEAN_SCALE`: as_stored
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
- `FILL_EMPTY_CLUSTERS`: false
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: none (mean of the values as stored, as Fig 5F)
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_MeanExpr_RAWvalues.csv` 1975beca; `Expt1234_MeanExpr_RAWvalues_RepsAvg.csv` 481b4d18; `Expt1234_MeanExpr_FC.csv` 950a5b7f; `Expt1234_MeanExpr_FC_RepsAvg.csv` c67a47b1; `Expt1234_MeanExpr_log2FC.csv` 6d75945c; `Expt1234_MeanExpr_log2FC_RepsAvg.csv` de4616f0

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
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 891
- `Table_Rows_Used`: 375
- `Reference_Rows`: 79
- `Comparisons`: 120
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 1caadaef; `Statistics_pvalues_wide.csv` a9907541; `Statistics_design.csv` 9dfe6041; `Statistics_litters.csv` 02533e0c; `Statistics_diagnostics.csv` 4fd96f5b

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Statistics_Folder`: 2026-11-01_14-04-59_Statistics_FC_welch_byCluster_Fig5F_BDNF
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 20
- `Comparisons`: 100
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 39
- `Significant_Primary`: 24
- `Significant_BH_All_Markers`: 31
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` fd01d514; `Timecourse_families.csv` 27f4a712; `Timecourse_pvalues_wide.csv` 106c3dde

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Cluster
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: [0.5, 1.5]
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: ["1", "15"]
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
- `Statistics_Folder`: 2026-11-01_14-05-38_Timecourse_FC_welch-welch_holm_Fig5F_BDNF
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: fd01d514b08ce351ef86c4eb825576e2
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Panels`: 1
- `Starred_Group_TimePoints`: 0
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` 814f6ce2; `Linegraph_Clusters_ByCluster-Grid.png` 09687ed7; `Plotted_data.csv` a334a20a; `Plotted_points.csv` 6cb6e8eb; `Plotted_axes.csv` bd16831c; `Plotted_style.csv` 876eec3f; `Figure_legend.txt` 3310ef46

## Run folder `2026-11-01_14-15-49_Linegraphs_ByCluster_p_primary_Fig5F_BDNF-Bot   [file label: 5F_BDNF-Bot]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `03_Statistics` | 2 |
| 3 | `04_TimecourseStatistics` | 1.2 |
| 4 | `01_PlotLinegraphs` | 1.12 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: mean_expression
- `MEAN_SCALE`: as_stored
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
- `FILL_EMPTY_CLUSTERS`: false
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: none (mean of the values as stored, as Fig 5F)
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_MeanExpr_RAWvalues.csv` 1975beca; `Expt1234_MeanExpr_RAWvalues_RepsAvg.csv` 481b4d18; `Expt1234_MeanExpr_FC.csv` 950a5b7f; `Expt1234_MeanExpr_FC_RepsAvg.csv` c67a47b1; `Expt1234_MeanExpr_log2FC.csv` 6d75945c; `Expt1234_MeanExpr_log2FC_RepsAvg.csv` de4616f0

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
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 891
- `Table_Rows_Used`: 375
- `Reference_Rows`: 79
- `Comparisons`: 120
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 1caadaef; `Statistics_pvalues_wide.csv` a9907541; `Statistics_design.csv` 9dfe6041; `Statistics_litters.csv` 02533e0c; `Statistics_diagnostics.csv` 4fd96f5b

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Statistics_Folder`: 2026-11-01_14-04-59_Statistics_FC_welch_byCluster_Fig5F_BDNF
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 20
- `Comparisons`: 100
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 39
- `Significant_Primary`: 24
- `Significant_BH_All_Markers`: 31
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` fd01d514; `Timecourse_families.csv` 27f4a712; `Timecourse_pvalues_wide.csv` 106c3dde

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Cluster
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: [0.5, 1.5]
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "BDNF"]
- `CLUSTERS`: ["2", "13", "4"]
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
- `Statistics_Folder`: 2026-11-01_14-05-38_Timecourse_FC_welch-welch_holm_Fig5F_BDNF
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: fd01d514b08ce351ef86c4eb825576e2
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Panels`: 1
- `Starred_Group_TimePoints`: 8
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` ee86c3d4; `Linegraph_Clusters_ByCluster-Grid.png` d396707d; `Plotted_data.csv` b3728684; `Plotted_points.csv` cb9870e4; `Plotted_axes.csv` bd16831c; `Plotted_style.csv` eae98fa1; `Figure_legend.txt` 44a96314

## Run folder `2026-11-01_14-17-22_Linegraphs_ByCluster_p_primary_Fig5F_ContDepr-Top   [file label: 5F_ContDepr-Top]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `03_Statistics` | 2 |
| 3 | `04_TimecourseStatistics` | 1.2 |
| 4 | `01_PlotLinegraphs` | 1.12 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: mean_expression
- `MEAN_SCALE`: as_stored
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
- `FILL_EMPTY_CLUSTERS`: false
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: none (mean of the values as stored, as Fig 5F)
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_MeanExpr_RAWvalues.csv` 1975beca; `Expt1234_MeanExpr_RAWvalues_RepsAvg.csv` 481b4d18; `Expt1234_MeanExpr_FC.csv` 950a5b7f; `Expt1234_MeanExpr_FC_RepsAvg.csv` c67a47b1; `Expt1234_MeanExpr_log2FC.csv` 6d75945c; `Expt1234_MeanExpr_log2FC_RepsAvg.csv` de4616f0

### 03_Statistics (v2)

- `STAT_TABLE`: FC
- `TEST`: welch
- `COMPARE_BY`: Cluster
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `REFERENCE_SCOPE`: pooled
- `RUN_DIAGNOSTICS`: true
- `LITTER_MAP`: 
- `TREATMENTS`: ["NoTreat", "ContDepr"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 891
- `Table_Rows_Used`: 336
- `Reference_Rows`: 79
- `Comparisons`: 120
- `Comparisons_NA`: 1
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 0f5a0f74; `Statistics_pvalues_wide.csv` ef0470bc; `Statistics_design.csv` d2933c41; `Statistics_litters.csv` 7188602a; `Statistics_diagnostics.csv` aed53dd0

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "ContDepr"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Statistics_Folder`: 2026-11-01_14-00-06_Statistics_FC_welch_byCluster_Fig5F_ContDepr
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 20
- `Comparisons`: 100
- `Comparisons_NA`: 1
- `Significant_Unadjusted`: 8
- `Significant_Primary`: 2
- `Significant_BH_All_Markers`: 2
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` ae8892de; `Timecourse_families.csv` 7182172d; `Timecourse_pvalues_wide.csv` c1822e6a

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Cluster
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: [0.5, 1.5]
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "ContDepr"]
- `CLUSTERS`: ["1", "15"]
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
- `Statistics_Folder`: 2026-11-01_14-00-46_Timecourse_FC_welch-welch_holm_Fig5F_ContDepr
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: ae8892deda76f3b871013e9ec9491a3f
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Panels`: 1
- `Starred_Group_TimePoints`: 0
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` e6c33bdb; `Linegraph_Clusters_ByCluster-Grid.png` b8394bde; `Plotted_data.csv` d1a0119f; `Plotted_points.csv` babef909; `Plotted_axes.csv` bd16831c; `Plotted_style.csv` 876eec3f; `Figure_legend.txt` 59165ac4

## Run folder `2026-11-01_14-17-44_Linegraphs_ByCluster_p_primary_Fig5F_ContDepr-Bot   [file label: 5F_ContDepr-Bot]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `03_Statistics` | 2 |
| 3 | `04_TimecourseStatistics` | 1.2 |
| 4 | `01_PlotLinegraphs` | 1.12 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: mean_expression
- `MEAN_SCALE`: as_stored
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
- `FILL_EMPTY_CLUSTERS`: false
- `Quantile_Type`: R default (type 7)
- `Back_Transform`: none (mean of the values as stored, as Fig 5F)
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_MeanExpr_RAWvalues.csv` 1975beca; `Expt1234_MeanExpr_RAWvalues_RepsAvg.csv` 481b4d18; `Expt1234_MeanExpr_FC.csv` 950a5b7f; `Expt1234_MeanExpr_FC_RepsAvg.csv` c67a47b1; `Expt1234_MeanExpr_log2FC.csv` 6d75945c; `Expt1234_MeanExpr_log2FC_RepsAvg.csv` de4616f0

### 03_Statistics (v2)

- `STAT_TABLE`: FC
- `TEST`: welch
- `COMPARE_BY`: Cluster
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `REFERENCE_SCOPE`: pooled
- `RUN_DIAGNOSTICS`: true
- `LITTER_MAP`: 
- `TREATMENTS`: ["NoTreat", "ContDepr"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 891
- `Table_Rows_Used`: 336
- `Reference_Rows`: 79
- `Comparisons`: 120
- `Comparisons_NA`: 1
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 0f5a0f74; `Statistics_pvalues_wide.csv` ef0470bc; `Statistics_design.csv` d2933c41; `Statistics_litters.csv` 7188602a; `Statistics_diagnostics.csv` aed53dd0

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "ContDepr"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Statistics_Folder`: 2026-11-01_14-00-06_Statistics_FC_welch_byCluster_Fig5F_ContDepr
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 20
- `Comparisons`: 100
- `Comparisons_NA`: 1
- `Significant_Unadjusted`: 8
- `Significant_Primary`: 2
- `Significant_BH_All_Markers`: 2
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` ae8892de; `Timecourse_families.csv` 7182172d; `Timecourse_pvalues_wide.csv` c1822e6a

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Cluster
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: [0.5, 1.5]
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "ContDepr"]
- `CLUSTERS`: ["2", "13", "4"]
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
- `Statistics_Folder`: 2026-11-01_14-00-46_Timecourse_FC_welch-welch_holm_Fig5F_ContDepr
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: ae8892deda76f3b871013e9ec9491a3f
- `Input_Table`: Expt1234_MeanExpr_FC.csv
- `Input_Table_MD5`: 950a5b7f3f281a2b4ed908cd98a0e785
- `Panels`: 1
- `Starred_Group_TimePoints`: 1
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByCluster-Grid.pdf` 0882d314; `Linegraph_Clusters_ByCluster-Grid.png` 6dc26e27; `Plotted_data.csv` 70dbec61; `Plotted_points.csv` 19fa78ac; `Plotted_axes.csv` bd16831c; `Plotted_style.csv` eae98fa1; `Figure_legend.txt` e0c8cf18
