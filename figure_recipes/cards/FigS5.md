# FigS5: how it was made

> **Note (S5):** the statistics reported in the paper come from the run labelled `S5(onlySignalingMarkers)` (2026-11-10, signaling markers only, 01 > 03 > 04). The figure file comes from the run labelled `S5` (2026-11-01, all markers); it was not redrawn because the stars did not change.


Run record(s): 2026-11-01_13-36-45_Linegraphs_ByTreatment_p_primary_FigS5   [file label: S5], 2026-11-10_09-22-34_Timecourse_RAWvalues_welch-welch_holm   [file label: S5(onlySignalingMarkers)]

## Run folder `2026-11-01_13-36-45_Linegraphs_ByTreatment_p_primary_FigS5   [file label: S5]`

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
- `ANALYSIS_BY`: PseudoBulk
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
- `Back_Transform`: 5 * sinh(x) on every Plot = 1 channel
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [0]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv` 1709bdac; `Expt1234_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` 15177dc2; `Expt1234_PercAboveThresh_95thPerc_FC.csv` 5326f8b8; `Expt1234_PercAboveThresh_95thPerc_FC_RepsAvg.csv` 5838e7e9; `Expt1234_PercAboveThresh_95thPerc_log2FC.csv` 73e8ffaf; `Expt1234_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` 41859761; `Thresholds.csv` 263ba6bd

### 03_Statistics (v2)

- `STAT_TABLE`: RAWvalues
- `TEST`: welch
- `COMPARE_BY`: Treatment
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
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv
- `Input_Table_MD5`: 1709bdac0e35d04dfea5db62ac0f6f22
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 45
- `Table_Rows_Used`: 17
- `Reference_Rows`: 4
- `Comparisons`: 240
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` b3f35030; `Statistics_pvalues_wide.csv` 723e94f6; `Statistics_design.csv` fbed732f; `Statistics_litters.csv` b8c8b57c; `Statistics_diagnostics.csv` aa9e87f8

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: RAWvalues
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv
- `Input_Table_MD5`: 1709bdac0e35d04dfea5db62ac0f6f22
- `Statistics_Folder`: 2026-10-29_16-58-29_Statistics_RAWvalues_welch_byTreatment_FigS5
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 40
- `Comparisons`: 200
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 77
- `Significant_Primary`: 36
- `Significant_BH_All_Markers`: 31
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 2bb46d02; `Timecourse_families.csv` a9916af5; `Timecourse_pvalues_wide.csv` 82b6c7e7

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Treatment
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "CompMed", "BDNF"]
- `CLUSTERS`: []
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
- `Statistics_Folder`: 2026-10-29_16-58-33_Timecourse_RAWvalues_welch-welch_holm_FigS5
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 2bb46d027a0664e2dcccc73ba174991d
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv
- `Input_Table_MD5`: 1709bdac0e35d04dfea5db62ac0f6f22
- `Panels`: 20
- `Starred_Group_TimePoints`: 19
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_PseudoBulk_ByTreatment-Grid.pdf` 735c901a; `Linegraph_PseudoBulk_ByTreatment-Grid.png` d6977684; `Plotted_data.csv` 7833943e; `Plotted_points.csv` 46ec1627; `Plotted_axes.csv` 97ef1771; `Plotted_style.csv` 72334d32; `Figure_legend.txt` dc188940

## Run folder `2026-11-10_09-22-34_Timecourse_RAWvalues_welch-welch_holm   [file label: S5(onlySignalingMarkers)]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `03_Statistics` | 2 |
| 3 | `04_TimecourseStatistics` | 1.2 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: percent_above
- `MEAN_SCALE`: null
- `ANALYSIS_BY`: PseudoBulk
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
- `Back_Transform`: 5 * sinh(x) on every Plot = 1 channel
- `Cells_In_Folder`: 680466
- `Cells_Analysed`: 575689
- `Clusters_Analysed`: [0]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv` 1709bdac; `Expt1234_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` 15177dc2; `Expt1234_PercAboveThresh_95thPerc_FC.csv` 5326f8b8; `Expt1234_PercAboveThresh_95thPerc_FC_RepsAvg.csv` 5838e7e9; `Expt1234_PercAboveThresh_95thPerc_log2FC.csv` 73e8ffaf; `Expt1234_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` 41859761; `Thresholds.csv` 263ba6bd

### 03_Statistics (v2)

- `STAT_TABLE`: RAWvalues
- `TEST`: welch
- `COMPARE_BY`: Treatment
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
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv
- `Input_Table_MD5`: 1709bdac0e35d04dfea5db62ac0f6f22
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 45
- `Table_Rows_Used`: 17
- `Reference_Rows`: 4
- `Comparisons`: 132
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` 54477014; `Statistics_pvalues_wide.csv` 48516a31; `Statistics_design.csv` fbed732f; `Statistics_litters.csv` b8c8b57c; `Statistics_diagnostics.csv` 17919d01

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: RAWvalues
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "CompMed"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv
- `Input_Table_MD5`: 1709bdac0e35d04dfea5db62ac0f6f22
- `Statistics_Folder`: 2026-11-10_09-22-28_Statistics_RAWvalues_welch_byTreatment_FigS5_onlySignalingMarkers
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 22
- `Comparisons`: 110
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 46
- `Significant_Primary`: 28
- `Significant_BH_All_Markers`: 25
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 67e93f23; `Timecourse_families.csv` a88bae08; `Timecourse_pvalues_wide.csv` 98e6ab45
