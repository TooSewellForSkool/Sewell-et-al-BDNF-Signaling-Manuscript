# FigS6C: how it was made

Run record(s): 2026-10-29_22-39-31_Bargraphs_p_primary_FigS6C_Neuronal   [file label: S6C_Neuronal], 2026-10-29_22-39-49_Bargraphs_p_primary_FigS6C_Non-neuronal   [file label: S6C_Non-neuronal]

## Run folder `2026-10-29_22-39-31_Bargraphs_p_primary_FigS6C_Neuronal   [file label: S6C_Neuronal]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `05_ConditionStatistics` | 2 |
| 3 | `03_PlotBargraphs` | 2 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: percent_above
- `MEAN_SCALE`: null
- `ANALYSIS_BY`: Clusters
- `CLUSTER_ROUND`: 1
- `CLUSTERS_MERGED`: true
- `Cluster_File`: cluster_R1_assigns_merged.csv
- `SAMPLE_SET_COLUMN`: Expt
- `SAMPLE_SET`: [5]
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
- `Cells_Analysed`: 104777
- `Clusters_Analysed`: [1, 2, 3, 4]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt5_PercAboveThresh_95thPerc_RAWvalues.csv` c25db68d; `Expt5_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` 25536845; `Expt5_PercAboveThresh_95thPerc_FC.csv` 19384a0e; `Expt5_PercAboveThresh_95thPerc_FC_RepsAvg.csv` c53f1f27; `Expt5_PercAboveThresh_95thPerc_log2FC.csv` 4a33b077; `Expt5_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` 7370c735; `Thresholds.csv` accbac44

### 05_ConditionStatistics (v2)

- `CONDITIONS`: ["ContDepr", "BDNFTest", "BDNF+K252aTest"]
- `CLUSTERS`: []
- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `Pairs`: ["ContDepr vs BDNFTest", "ContDepr vs BDNF+K252aTest", "BDNFTest vs BDNF+K252aTest"]
- `Primary_Family`: the 3 pairs of one marker in one cluster/category (default: Welch t-tests + Holm for every marker; Tukey in p_Tukey)
- `Alternative_Family`: Benjamini-Hochberg over all markers x pairs of one cluster/category, on p_unadjusted
- `Input_Table`: Expt5_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 19384a0e67d2afa4aa75b37bcfb0ad3e
- `Families_Analysed`: 88
- `Comparisons`: 264
- `Comparisons_NA`: 0
- `Significant_Tukey`: 81
- `Significant_Unadjusted`: 80
- `Significant_Primary`: 60
- `Significant_BH_All_Markers`: 42
- `Random_Seed`: not used (no random numbers in this script)
- Output files (MD5): `Conditions_long.csv` fdc7fde6; `Conditions_families.csv` 140ae946; `Conditions_groupstats.csv` 51466ef0; `Conditions_input_values.csv` 7c621b38

### 03_PlotBargraphs (v2)

- `P_COLUMN`: p_primary
- `CLUSTERS`: ["2"]
- `MIN_N_FOR_STARS`: 3
- `SHOW_NS`: true
- `STEP_INCREASE`: 0.15
- `COLOR_OPTION`: C
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 6
- `PANEL_SIZE_IN`: 5
- `Statistics_Folder`: 2026-10-29_17-16-53_Conditions_welch-welch_holm_FigS6C
- `Panels`: 5
- `Figures`: 5
- `Starred_Comparisons`: 12
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `BarPlot_Clusters_2_pERK-Grid.pdf` ce0b8597; `BarPlot_Clusters_2_pERK-Grid.png` 478707f1; `BarPlot_Clusters_2_pAKt-Grid.pdf` 8a8ac37d; `BarPlot_Clusters_2_pAKt-Grid.png` f2fccba1; `BarPlot_Clusters_2_pCreb-Grid.pdf` c275301f; `BarPlot_Clusters_2_pCreb-Grid.png` 866336e5; `BarPlot_Clusters_2_pS6-Grid.pdf` f41e6acb; `BarPlot_Clusters_2_pS6-Grid.png` 9989a85b

## Run folder `2026-10-29_22-39-49_Bargraphs_p_primary_FigS6C_Non-neuronal   [file label: S6C_Non-neuronal]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `05_ConditionStatistics` | 2 |
| 3 | `03_PlotBargraphs` | 2 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: percent_above
- `MEAN_SCALE`: null
- `ANALYSIS_BY`: Clusters
- `CLUSTER_ROUND`: 1
- `CLUSTERS_MERGED`: true
- `Cluster_File`: cluster_R1_assigns_merged.csv
- `SAMPLE_SET_COLUMN`: Expt
- `SAMPLE_SET`: [5]
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
- `Cells_Analysed`: 104777
- `Clusters_Analysed`: [1, 2, 3, 4]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt5_PercAboveThresh_95thPerc_RAWvalues.csv` c25db68d; `Expt5_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` 25536845; `Expt5_PercAboveThresh_95thPerc_FC.csv` 19384a0e; `Expt5_PercAboveThresh_95thPerc_FC_RepsAvg.csv` c53f1f27; `Expt5_PercAboveThresh_95thPerc_log2FC.csv` 4a33b077; `Expt5_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` 7370c735; `Thresholds.csv` accbac44

### 05_ConditionStatistics (v2)

- `CONDITIONS`: ["ContDepr", "BDNFTest", "BDNF+K252aTest"]
- `CLUSTERS`: []
- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `Pairs`: ["ContDepr vs BDNFTest", "ContDepr vs BDNF+K252aTest", "BDNFTest vs BDNF+K252aTest"]
- `Primary_Family`: the 3 pairs of one marker in one cluster/category (default: Welch t-tests + Holm for every marker; Tukey in p_Tukey)
- `Alternative_Family`: Benjamini-Hochberg over all markers x pairs of one cluster/category, on p_unadjusted
- `Input_Table`: Expt5_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 19384a0e67d2afa4aa75b37bcfb0ad3e
- `Families_Analysed`: 88
- `Comparisons`: 264
- `Comparisons_NA`: 0
- `Significant_Tukey`: 81
- `Significant_Unadjusted`: 80
- `Significant_Primary`: 60
- `Significant_BH_All_Markers`: 42
- `Random_Seed`: not used (no random numbers in this script)
- Output files (MD5): `Conditions_long.csv` fdc7fde6; `Conditions_families.csv` 140ae946; `Conditions_groupstats.csv` 51466ef0; `Conditions_input_values.csv` 7c621b38

### 03_PlotBargraphs (v2)

- `P_COLUMN`: p_primary
- `CLUSTERS`: ["3"]
- `MIN_N_FOR_STARS`: 3
- `SHOW_NS`: true
- `STEP_INCREASE`: 0.15
- `COLOR_OPTION`: C
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 6
- `PANEL_SIZE_IN`: 5
- `Statistics_Folder`: 2026-10-29_17-16-53_Conditions_welch-welch_holm_FigS6C
- `Panels`: 5
- `Figures`: 5
- `Starred_Comparisons`: 9
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `BarPlot_Clusters_3_pERK-Grid.pdf` 70dd3c9e; `BarPlot_Clusters_3_pERK-Grid.png` f1a6325e; `BarPlot_Clusters_3_pAKt-Grid.pdf` c3e12c42; `BarPlot_Clusters_3_pAKt-Grid.png` 5591e6d9; `BarPlot_Clusters_3_pCreb-Grid.pdf` dc043e17; `BarPlot_Clusters_3_pCreb-Grid.png` 014a8923; `BarPlot_Clusters_3_pS6-Grid.pdf` 5847bf08; `BarPlot_Clusters_3_pS6-Grid.png` 90e5e40a
