# Fig5A: how it was made

Run record(s): 2026-10-29_22-48-05_CategoryHistograms_Fig5A   [file label: 5A]

## Run folder `2026-10-29_22-48-05_CategoryHistograms_Fig5A   [file label: 5A]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `02_AssignThresholdCategories` | 2 |
| 2 | `04_PlotCategoryHistograms` | 1.1 |

### 02_AssignThresholdCategories (v2)

- `SELECTED_MARKERS`: {"x": "TrkB_i_174Yb"}
- `CUTS_BY_QUANTILE`: false
- `CUT_QUANTILES`: []
- `Cut_Values_Used`: {"x": [0.15, 1.4, 3.7]}
- `CUTS_JUSTIFIED`: true
- `CUT_SHIFT`: 0
- `Rule`: level = 1 + number of cut values <= the cell's value (a cell on a cut goes to the higher level)
- `Cells`: 680466
- `Random_Seed`: not used in this script
- Output files (MD5): `ThresholdCategory.csv` 802760a4; `Category_Definitions.csv` 6981f9c9; `Category_Percentages.csv` 4a386a4d; `Histogram_Data.csv` 88f456bf; `Category_Cells_By_File.csv` f4e425ee

### 04_PlotCategoryHistograms (v1.1)

- `SCALE_NOTE`: arcsinh(x / 5)
- `LEVEL_NAMES`: []
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_WIDTH_IN`: 6
- `PANEL_HEIGHT_IN`: 4
- `Category_Folder`: 2026-10-29_15-48-58_ThreshCat_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7
- `Cut_Shift_In_Category_Run`: 0
- `Cells`: 680466
- `Random_Seed`: not used in this script
- Output files (MD5): `CategoryHistograms.pdf` bb261938; `CategoryHistograms.png` 52bd1602; `Plotted_data.csv` 1485f97b; `Plotted_levels.csv` d085c3a7; `Figure_legend.txt` 045990c8
