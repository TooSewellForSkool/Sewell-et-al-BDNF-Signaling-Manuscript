# FigS4B: how it was made

Run record(s): 2026-11-01_13-33-27_Linegraphs_ByTreatment_p_primary_FigS4B   [file label: S4B]

## Run folder `2026-11-01_13-33-27_Linegraphs_ByTreatment_p_primary_FigS4B   [file label: S4B]`

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
- `TREATMENTS`: ["NoTreat", "ContDepr", "BDNF"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Adjustment`: none (as published)
- `Test_Detail`: stats::t.test(x, y, paired = FALSE)  [Welch; the published method]
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv
- `Input_Table_MD5`: 1709bdac0e35d04dfea5db62ac0f6f22
- `Upstream_Module`: 01_PercAboveThreshold
- `Upstream_Script_Version`: 4
- `Table_Rows_Total`: 45
- `Table_Rows_Used`: 32
- `Reference_Rows`: 4
- `Comparisons`: 440
- `Comparisons_NA`: 0
- `Litter_Map_MD5`: d96da78cda3d6c0966a5714271c65621
- `Random_Seed`: not used in this script
- Output files (MD5): `Statistics_long.csv` d1c64df4; `Statistics_pvalues_wide.csv` 93bf28b4; `Statistics_design.csv` 5d1b8a35; `Statistics_litters.csv` 4248c20a; `Statistics_diagnostics.csv` bb17dd62

### 04_TimecourseStatistics (v1.2)

- `OMNIBUS`: welch
- `POSTHOC`: welch_holm
- `VAR_RATIO_LIMIT`: 3
- `SEED`: 1
- `STAT_TABLE`: RAWvalues
- `REFERENCE_TREATMENT`: NoTreat
- `REFERENCE_TIMEPOINT`: 0
- `TREATMENTS`: ["NoTreat", "ContDepr", "BDNF"]
- `CLUSTERS`: []
- `CATEGORIES`: []
- `Primary_Family`: one treatment x one marker (x one cluster/category): time points against the reference
- `Alternative_Family`: Benjamini-Hochberg over all markers x time points of one treatment (x cluster/category), on p_unadjusted
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv
- `Input_Table_MD5`: 1709bdac0e35d04dfea5db62ac0f6f22
- `Statistics_Folder`: 2026-10-29_16-58-00_Statistics_RAWvalues_welch_byTreatment_FigS4B
- `Statistics_Script_Version`: 2
- `Families_Analysed`: 80
- `Comparisons`: 400
- `Comparisons_NA`: 0
- `Significant_Unadjusted`: 77
- `Significant_Primary`: 35
- `Significant_BH_All_Markers`: 12
- `Random_Seed`: 1
- Output files (MD5): `Timecourse_long.csv` 2a7eb719; `Timecourse_families.csv` 7dee08d1; `Timecourse_pvalues_wide.csv` 3046f237

### 01_PlotLinegraphs (v1.12)

- `P_COLUMN`: p_primary
- `PLOT_BY`: Treatment
- `SD_STYLE`: ribbon
- `ERRORBAR_MIN_N`: 3
- `MIN_N_FOR_STARS`: 3
- `Y_LIMITS`: []
- `STAR_STEP`: auto
- `NO_STARS_TREATMENTS`: []
- `TREATMENTS`: ["NoTreat", "ContDepr", "BDNF"]
- `CLUSTERS`: []
- `CATEGORY`: none
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_COLUMNS`: 3
- `PANEL_SIZE_IN`: 5
- `CAPTION_ON_FIGURE`: false
- `BACKGROUND`: white
- `TRANSPARENT_ON`: white
- `Palette_File`: Palette_Colors.csv
- `POINT_SHAPES`: auto
- `Statistics_Folder`: 2026-10-29_16-58-08_Timecourse_RAWvalues_welch-welch_holm_FigS4B
- `Statistics_Source`: 04_TimecourseStatistics
- `Statistics_File_MD5`: 2a7eb7198e8114e56e57694877186da1
- `Input_Table`: Expt1234_PercAboveThresh_95thPerc_RAWvalues.csv
- `Input_Table_MD5`: 1709bdac0e35d04dfea5db62ac0f6f22
- `Panels`: 15
- `Starred_Group_TimePoints`: 5
- `N_Mismatch_With_Statistics`: 0
- `PNG_Device`: ragg::agg_png
- `Random_Seed`: not used in this script
- Output files (MD5): `Linegraph_PseudoBulk_ByTreatment-Grid.pdf` 9b21dac0; `Linegraph_PseudoBulk_ByTreatment-Grid.png` 79ba9c77; `Plotted_data.csv` f7d196f0; `Plotted_points.csv` 31557e5f; `Plotted_axes.csv` 26d2555f; `Plotted_style.csv` 0391c2ae; `Figure_legend.txt` 5f797838
