# FigS8: how it was made

Run record(s): 2026-11-03_12-23-05_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl1], 2026-11-03_12-24-31_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl2], 2026-11-03_12-25-13_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl4], 2026-11-03_12-26-23_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl7], 2026-11-03_12-27-05_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl13]

## Run folder `2026-11-03_12-23-05_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl1]`

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
- `COMPARE_BY`: Treatment
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `REFERENCE_SCOPE`: pooled
- `RUN_DIAGNOSTICS`: true
- `LITTER_MAP`: 
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: [1]
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 900
- `Table_Rows_Used`: 45
- `Reference_Rows`: 4
- `Comparisons`: 352
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` c3b92c3c; `Statistics_pvalues_wide.csv` d866733c; `Statistics_design.csv` 6bde00c0; `Statistics_litters.csv` 5dfe0b4b; `Statistics_diagnostics.csv` 5854e0af

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["1"]
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Statistics_Folder`: 2026-11-03_12-22-20_Statistics_FC_welch_byTreatment_Cl1_FigS8
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 66
- `Comparisons`: 330
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 66
- `Significant_Primary`: 21
- `Significant_BH_All_Markers`: 18
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` c677389e; `Timecourse_families.csv` 8c52dc2b; `Timecourse_pvalues_wide.csv` b698e65f

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Treatment
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["1"]
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
- `Statistics_Folder`: 2026-11-03_12-22-38_Timecourse_FC_welch-welch_holm_FigS8
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: c677389ece015df209a06e4c6845f4d0
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Panels`: 6
- `Starred_Group_TimePoints`: 9
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByTreatment-Grid.pdf` e24a9792; `Linegraph_Clusters_ByTreatment-Grid.png` 82335098; `Plotted_data.csv` da421612; `Plotted_points.csv` 69f8046a; `Plotted_axes.csv` e067fe79; `Plotted_style.csv` e97c56ab; `Figure_legend.txt` 030fdea2

## Run folder `2026-11-03_12-24-31_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl2]`

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
- `COMPARE_BY`: Treatment
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `REFERENCE_SCOPE`: pooled
- `RUN_DIAGNOSTICS`: true
- `LITTER_MAP`: 
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: [2]
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 900
- `Table_Rows_Used`: 45
- `Reference_Rows`: 4
- `Comparisons`: 352
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` e264071b; `Statistics_pvalues_wide.csv` 51b2b8db; `Statistics_design.csv` 6bde00c0; `Statistics_litters.csv` 5dfe0b4b; `Statistics_diagnostics.csv` 12cc148f

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["2"]
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Statistics_Folder`: 2026-11-03_12-23-40_Statistics_FC_welch_byTreatment_Cl2_FigS8
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 66
- `Comparisons`: 330
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 82
- `Significant_Primary`: 51
- `Significant_BH_All_Markers`: 38
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` fba6aa97; `Timecourse_families.csv` acd9ed56; `Timecourse_pvalues_wide.csv` 953d1f65

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Treatment
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["2"]
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
- `Statistics_Folder`: 2026-11-03_12-23-49_Timecourse_FC_welch-welch_holm_FigS8
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: fba6aa9714b04498c0b4aedf557bd2b0
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Panels`: 6
- `Starred_Group_TimePoints`: 22
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByTreatment-Grid.pdf` af431af3; `Linegraph_Clusters_ByTreatment-Grid.png` 33ce5aea; `Plotted_data.csv` 245e53e7; `Plotted_points.csv` f5832158; `Plotted_axes.csv` e067fe79; `Plotted_style.csv` e97c56ab; `Figure_legend.txt` 02353abf

## Run folder `2026-11-03_12-25-13_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl4]`

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
- `COMPARE_BY`: Treatment
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `REFERENCE_SCOPE`: pooled
- `RUN_DIAGNOSTICS`: true
- `LITTER_MAP`: 
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: [4]
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 900
- `Table_Rows_Used`: 45
- `Reference_Rows`: 4
- `Comparisons`: 352
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 003afc1a; `Statistics_pvalues_wide.csv` 5e744f30; `Statistics_design.csv` 6bde00c0; `Statistics_litters.csv` 5dfe0b4b; `Statistics_diagnostics.csv` 97a2aaab

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["4"]
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Statistics_Folder`: 2026-11-03_12-24-56_Statistics_FC_welch_byTreatment_Cl4_FigS8
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 66
- `Comparisons`: 330
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 76
- `Significant_Primary`: 33
- `Significant_BH_All_Markers`: 22
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 6124252a; `Timecourse_families.csv` 3464d84d; `Timecourse_pvalues_wide.csv` 37b2f4fd

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Treatment
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["4"]
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
- `Statistics_Folder`: 2026-11-03_12-25-02_Timecourse_FC_welch-welch_holm_FigS8
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 6124252afe4d6111959a3d2cfecff4c6
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Panels`: 6
- `Starred_Group_TimePoints`: 20
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByTreatment-Grid.pdf` e4442f0c; `Linegraph_Clusters_ByTreatment-Grid.png` 1345c66b; `Plotted_data.csv` 4c45d863; `Plotted_points.csv` 95be560c; `Plotted_axes.csv` e067fe79; `Plotted_style.csv` e97c56ab; `Figure_legend.txt` 94fd758c

## Run folder `2026-11-03_12-26-23_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl7]`

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
- `COMPARE_BY`: Treatment
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `REFERENCE_SCOPE`: pooled
- `RUN_DIAGNOSTICS`: true
- `LITTER_MAP`: 
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: [7]
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 900
- `Table_Rows_Used`: 45
- `Reference_Rows`: 4
- `Comparisons`: 352
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` fcbd1604; `Statistics_pvalues_wide.csv` 7364918a; `Statistics_design.csv` 6bde00c0; `Statistics_litters.csv` 5dfe0b4b; `Statistics_diagnostics.csv` 44e0362c

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["7"]
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Statistics_Folder`: 2026-11-03_12-25-40_Statistics_FC_welch_byTreatment_Cl7_FigS8
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 66
- `Comparisons`: 330
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 72
- `Significant_Primary`: 37
- `Significant_BH_All_Markers`: 17
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 6dfa24fb; `Timecourse_families.csv` 08a816d6; `Timecourse_pvalues_wide.csv` 5995a815

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Treatment
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["7"]
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
- `Statistics_Folder`: 2026-11-03_12-25-45_Timecourse_FC_welch-welch_holm_FigS8
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 6dfa24fb483fa3746db5a5cbec2708b2
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Panels`: 6
- `Starred_Group_TimePoints`: 17
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByTreatment-Grid.pdf` aa408074; `Linegraph_Clusters_ByTreatment-Grid.png` ff242c78; `Plotted_data.csv` 9a4e2a2e; `Plotted_points.csv` ad0c2bbc; `Plotted_axes.csv` e067fe79; `Plotted_style.csv` e97c56ab; `Figure_legend.txt` 5a9b1be3

## Run folder `2026-11-03_12-27-05_Linegraphs_ByTreatment_p_primary_FigS8   [file label: S8_Cl13]`

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
- `COMPARE_BY`: Treatment
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `REFERENCE_SCOPE`: pooled
- `RUN_DIAGNOSTICS`: true
- `LITTER_MAP`: 
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: [13]
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 900
- `Table_Rows_Used`: 45
- `Reference_Rows`: 4
- `Comparisons`: 352
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` f5636650; `Statistics_pvalues_wide.csv` 7c0cc89d; `Statistics_design.csv` 6bde00c0; `Statistics_litters.csv` 5dfe0b4b; `Statistics_diagnostics.csv` f0d0b3f9

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: FC
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["13"]
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Statistics_Folder`: 2026-11-03_12-26-45_Statistics_FC_welch_byTreatment_Cl13_FigS8
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 66
- `Comparisons`: 330
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 62
- `Significant_Primary`: 35
- `Significant_BH_All_Markers`: 18
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` f1e9ba06; `Timecourse_families.csv` bb7bf81b; `Timecourse_pvalues_wide.csv` 33d19953

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Treatment
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "ContDepr", "CompMed", "BDNF"]
- `CLUSTERS`: ["13"]
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
- `Statistics_Folder`: 2026-11-03_12-26-51_Timecourse_FC_welch-welch_holm_FigS8
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: f1e9ba0694fab7d835146a68aca64ff6
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: a5330fbb7aef84b5bd94c173daec2330
- `Panels`: 6
- `Starred_Group_TimePoints`: 16
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_Clusters_ByTreatment-Grid.pdf` 30f6c76f; `Linegraph_Clusters_ByTreatment-Grid.png` 0d1861f5; `Plotted_data.csv` dcc41960; `Plotted_points.csv` e0069364; `Plotted_axes.csv` e067fe79; `Plotted_style.csv` e97c56ab; `Figure_legend.txt` d297eb77
