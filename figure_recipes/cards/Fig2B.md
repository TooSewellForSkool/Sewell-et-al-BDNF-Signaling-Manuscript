# Fig2B: how it was made

Run record(s): 2026-10-29_22-34-26_Bargraphs_p_primary_Fig2B   [file label: 2B]

## Run folder `2026-10-29_22-34-26_Bargraphs_p_primary_Fig2B   [file label: 2B]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `01_PercAboveThreshold` | 4 |
| 2 | `05_ConditionStatistics` | 2 |
| 3 | `03_PlotBargraphs` | 2 |

### 01_PercAboveThreshold (v4)

- `MEASURE`: percent_above
- `MEAN_SCALE`: null
- `ANALYSIS_BY`: PseudoBulk
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
- `Clusters_Analysed`: [0]
- `Markers_Derived_Rule`: published rule (A4)
- `Random_Seed`: not used in this script
- Output files (MD5): `Expt5_PercAboveThresh_95thPerc_RAWvalues.csv` 8503a13f; `Expt5_PercAboveThresh_95thPerc_RAWvalues_RepsAvg.csv` 38507f1d; `Expt5_PercAboveThresh_95thPerc_FC.csv` c0a08f35; `Expt5_PercAboveThresh_95thPerc_FC_RepsAvg.csv` b71b80be; `Expt5_PercAboveThresh_95thPerc_log2FC.csv` 8f81efc7; `Expt5_PercAboveThresh_95thPerc_log2FC_RepsAvg.csv` 7f1b6b13; `Thresholds.csv` 1fd616e8

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
- `Input_Table_MD5`: c0a08f35ba28253dd65a87da27516617
- `Families_Analysed`: 22
- `Comparisons`: 66
- `Comparisons_NA`: 0
- `Significant_Tukey`: 30
- `Significant_Unadjusted`: 29
- `Significant_Primary`: 22
- `Significant_BH_All_Markers`: 19
- `Random_Seed`: not used (no random numbers in this script)
- Output files (MD5): `Conditions_long.csv` 3f5ddd94; `Conditions_families.csv` 00ff2562; `Conditions_groupstats.csv` c1e1bfb3; `Conditions_input_values.csv` 3cd25812

### 03_PlotBargraphs (v2)

- `P_COLUMN`: p_primary
- `CLUSTERS`: ["0"]
- `MIN_N_FOR_STARS`: 3
- `SHOW_NS`: true
- `STEP_INCREASE`: 0.15
- `COLOR_OPTION`: C
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 6
- `PANEL_SIZE_IN`: 5
- `Statistics_Folder`: 2026-10-29_17-15-26_Conditions_welch-welch_holm_Fig2B
- `Panels`: 5
- `Figures`: 5
- `Starred_Comparisons`: 11
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `BarPlot_Clusters_0_pERK-Grid.pdf` 7bbc38bd; `BarPlot_Clusters_0_pERK-Grid.png` 9910ac18; `BarPlot_Clusters_0_pAKt-Grid.pdf` 68d7ad7c; `BarPlot_Clusters_0_pAKt-Grid.png` 485df5d2; `BarPlot_Clusters_0_pCreb-Grid.pdf` 6184068b; `BarPlot_Clusters_0_pCreb-Grid.png` 7a6c1ca8; `BarPlot_Clusters_0_pS6-Grid.pdf` e7a00a4d; `BarPlot_Clusters_0_pS6-Grid.png` 039b8f2e
