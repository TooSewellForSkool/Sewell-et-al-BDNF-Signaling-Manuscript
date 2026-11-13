# FigS7B: how it was made

Run record(s): 2026-10-29_22-42-00_Bargraphs_p_primary_FigS7B_pAkt   [file label: S7B_pAkt], 2026-10-29_22-43-51_Bargraphs_p_primary_FigS7B_pCreb   [file label: S7B_pCreb], 2026-10-29_22-44-09_Bargraphs_p_primary_FigS7B_pGSK3b   [file label: S7B_pGSK3b]

## Run folder `2026-10-29_22-42-00_Bargraphs_p_primary_FigS7B_pAkt   [file label: S7B_pAkt]`

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
- `CLUSTERS_MERGED`: false
- `Cluster_File`: cluster_R1_assigns.csv
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
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt5_PercAboveThresh_95thPerc_RAWvalues.csv` 9b6646cb; `Expt5_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` ff61286f; `Expt5_PercAboveThresh_95thPerc_FC.csv` 32cdcdd4; `Expt5_PercAboveThresh_95thPerc_FC_RepsAvg.csv` c0512b02; `Expt5_PercAboveThresh_95thPerc_log2FC.csv` 039dd74c; `Expt5_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` fca845f0; `Thresholds.csv` bb34d5a8

### 05_ConditionStatistics (v2)

- `CONDITIONS`: ["ContDepr", "BDNFTest", "BDNF+K252aTest"]
- `CLUSTERS`: ["1", "2", "4", "3", "17", "18"]
- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `Pairs`: ["ContDepr vs BDNFTest", "ContDepr vs BDNF+K252aTest", "BDNFTest vs BDNF+K252aTest"]
- `Primary_Family`: the 3 pairs of one marker in one cluster/category (default: Welch t-tests + Holm for every marker; Tukey in p_Tukey)
- `Alternative_Family`: Benjamini-Hochberg over all markers x pairs of one cluster/category, on p_unadjusted
- `Input_Table`: Expt5_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 32cdcdd41baa9fa6eab0562ba3131b2a
- `Families_Analysed`: 132
- `Comparisons`: 396
- `Comparisons_NA`: 0
- `Significant_Tukey`: 105
- `Significant_Unadjusted`: 101
- `Significant_Primary`: 64
- `Significant_BH_All_Markers`: 31
- `Random_Seed`: not used (no random numbers in this script)
- Output files (MD5): `Conditions_long.csv` d6f7cce1; `Conditions_families.csv` 0350073b; `Conditions_groupstats.csv` 961b4d72; `Conditions_input_values.csv` 716ad0b5

### 03_PlotBargraphs (v2)

- `P_COLUMN`: p_primary
- `CLUSTERS`: ["1", "2", "4", "3", "17", "18"]
- `MIN_N_FOR_STARS`: 3
- `SHOW_NS`: true
- `STEP_INCREASE`: 0.15
- `COLOR_OPTION`: C
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 6
- `PANEL_SIZE_IN`: 5
- `Statistics_Folder`: 2026-10-29_17-16-19_Conditions_welch-welch_holm_Fig3F
- `Panels`: 6
- `Figures`: 1
- `Starred_Comparisons`: 10
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `BarPlot_Clusters_1-2-4-3-17-18_pAKt-Grid.pdf` 484893aa; `BarPlot_Clusters_1-2-4-3-17-18_pAKt-Grid.png` a567fdfb; `Plotted_data.csv` ff83b903; `Figure_legend.txt` 9ecc3338

## Run folder `2026-10-29_22-43-51_Bargraphs_p_primary_FigS7B_pCreb   [file label: S7B_pCreb]`

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
- `CLUSTERS_MERGED`: false
- `Cluster_File`: cluster_R1_assigns.csv
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
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt5_PercAboveThresh_95thPerc_RAWvalues.csv` 9b6646cb; `Expt5_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` ff61286f; `Expt5_PercAboveThresh_95thPerc_FC.csv` 32cdcdd4; `Expt5_PercAboveThresh_95thPerc_FC_RepsAvg.csv` c0512b02; `Expt5_PercAboveThresh_95thPerc_log2FC.csv` 039dd74c; `Expt5_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` fca845f0; `Thresholds.csv` bb34d5a8

### 05_ConditionStatistics (v2)

- `CONDITIONS`: ["ContDepr", "BDNFTest", "BDNF+K252aTest"]
- `CLUSTERS`: ["1", "2", "4", "3", "17", "18"]
- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `Pairs`: ["ContDepr vs BDNFTest", "ContDepr vs BDNF+K252aTest", "BDNFTest vs BDNF+K252aTest"]
- `Primary_Family`: the 3 pairs of one marker in one cluster/category (default: Welch t-tests + Holm for every marker; Tukey in p_Tukey)
- `Alternative_Family`: Benjamini-Hochberg over all markers x pairs of one cluster/category, on p_unadjusted
- `Input_Table`: Expt5_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 32cdcdd41baa9fa6eab0562ba3131b2a
- `Families_Analysed`: 132
- `Comparisons`: 396
- `Comparisons_NA`: 0
- `Significant_Tukey`: 105
- `Significant_Unadjusted`: 101
- `Significant_Primary`: 64
- `Significant_BH_All_Markers`: 31
- `Random_Seed`: not used (no random numbers in this script)
- Output files (MD5): `Conditions_long.csv` d6f7cce1; `Conditions_families.csv` 0350073b; `Conditions_groupstats.csv` 961b4d72; `Conditions_input_values.csv` 716ad0b5

### 03_PlotBargraphs (v2)

- `P_COLUMN`: p_primary
- `CLUSTERS`: ["1", "2", "4", "3", "17", "18"]
- `MIN_N_FOR_STARS`: 3
- `SHOW_NS`: true
- `STEP_INCREASE`: 0.15
- `COLOR_OPTION`: C
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 6
- `PANEL_SIZE_IN`: 5
- `Statistics_Folder`: 2026-10-29_17-16-19_Conditions_welch-welch_holm_Fig3F
- `Panels`: 6
- `Figures`: 1
- `Starred_Comparisons`: 6
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `BarPlot_Clusters_1-2-4-3-17-18_pCreb-Grid.pdf` b74fb47e; `BarPlot_Clusters_1-2-4-3-17-18_pCreb-Grid.png` 75202e2d; `Plotted_data.csv` a33b7dd8; `Figure_legend.txt` 9ecc3338

## Run folder `2026-10-29_22-44-09_Bargraphs_p_primary_FigS7B_pGSK3b   [file label: S7B_pGSK3b]`

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
- `CLUSTERS_MERGED`: false
- `Cluster_File`: cluster_R1_assigns.csv
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
- `Clusters_Analysed`: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt5_PercAboveThresh_95thPerc_RAWvalues.csv` 9b6646cb; `Expt5_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` ff61286f; `Expt5_PercAboveThresh_95thPerc_FC.csv` 32cdcdd4; `Expt5_PercAboveThresh_95thPerc_FC_RepsAvg.csv` c0512b02; `Expt5_PercAboveThresh_95thPerc_log2FC.csv` 039dd74c; `Expt5_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` fca845f0; `Thresholds.csv` bb34d5a8

### 05_ConditionStatistics (v2)

- `CONDITIONS`: ["ContDepr", "BDNFTest", "BDNF+K252aTest"]
- `CLUSTERS`: ["1", "2", "4", "3", "17", "18"]
- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `Pairs`: ["ContDepr vs BDNFTest", "ContDepr vs BDNF+K252aTest", "BDNFTest vs BDNF+K252aTest"]
- `Primary_Family`: the 3 pairs of one marker in one cluster/category (default: Welch t-tests + Holm for every marker; Tukey in p_Tukey)
- `Alternative_Family`: Benjamini-Hochberg over all markers x pairs of one cluster/category, on p_unadjusted
- `Input_Table`: Expt5_PercAboveThresh_95thPerc_FC.csv
- `Input_Table_MD5`: 32cdcdd41baa9fa6eab0562ba3131b2a
- `Families_Analysed`: 132
- `Comparisons`: 396
- `Comparisons_NA`: 0
- `Significant_Tukey`: 105
- `Significant_Unadjusted`: 101
- `Significant_Primary`: 64
- `Significant_BH_All_Markers`: 31
- `Random_Seed`: not used (no random numbers in this script)
- Output files (MD5): `Conditions_long.csv` d6f7cce1; `Conditions_families.csv` 0350073b; `Conditions_groupstats.csv` 961b4d72; `Conditions_input_values.csv` 716ad0b5

### 03_PlotBargraphs (v2)

- `P_COLUMN`: p_primary
- `CLUSTERS`: ["1", "2", "4", "3", "17", "18"]
- `MIN_N_FOR_STARS`: 3
- `SHOW_NS`: true
- `STEP_INCREASE`: 0.15
- `COLOR_OPTION`: C
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 6
- `PANEL_SIZE_IN`: 5
- `Statistics_Folder`: 2026-10-29_17-16-19_Conditions_welch-welch_holm_Fig3F
- `Panels`: 6
- `Figures`: 1
- `Starred_Comparisons`: 4
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `BarPlot_Clusters_1-2-4-3-17-18_pGSK3b-Grid.pdf` a1e1247e; `BarPlot_Clusters_1-2-4-3-17-18_pGSK3b-Grid.png` bd3be362; `Plotted_data.csv` 9a51d3e3; `Figure_legend.txt` 9ecc3338
