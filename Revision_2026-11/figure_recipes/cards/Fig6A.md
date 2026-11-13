# Fig6A: how it was made

Run record(s): 2026-10-29_22-48-59_CategoryHistograms_Fig6A   [file label: 6A]

## Run folder `2026-10-29_22-48-59_CategoryHistograms_Fig6A   [file label: 6A]`

Chain (each step reads the previous step's output; settings are what to type into the USER SETTINGS block of that script):

| Step | Script | Version |
|---|---|---|
| 1 | `02_AssignThresholdCategories` | 2 |
| 2 | `04_PlotCategoryHistograms` | 1.1 |

### 02_AssignThresholdCategories (v2)

- `SELECTED_MARKERS`: {"x": "p75NTR_161Dy", "y": "TrkB_i_174Yb"}
- `CUTS_BY_QUANTILE`: false
- `CUT_QUANTILES`: []
- `Cut_Values_Used`: {"x": [0.45, 1.45, 3.95], "y": [0.15, 1.4, 3.7]}
- `CUTS_JUSTIFIED`: true
- `CUT_SHIFT`: 0
- `Rule`: level = 1 + number of cut values <= the cell's value (a cell on a cut goes to the higher level)
- `Cells`: 680466
- `Random_Seed`: not used in this script
- Output files (MD5): `ThresholdCategory.csv` b5a8c290; `Category_Definitions.csv` 8045d913; `Category_Percentages.csv` 4c0b19c2; `Histogram_Data.csv` dc3ea628; `Category_Cells_By_File.csv` 4bd2cd52

### 04_PlotCategoryHistograms (v1.1)

- `SCALE_NOTE`: arcsinh(x / 5)
- `LEVEL_NAMES`: []
- `OUTPUT_FORMATS`: ["pdf", "png"]
- `PNG_DPI`: 300
- `PANEL_WIDTH_IN`: 6
- `PANEL_HEIGHT_IN`: 4
- `Category_Folder`: 2026-10-29_15-50-12_ThreshCat_p75NTR_161Dy-Gr0.45-Gr1.45-Gr3.95_TrkB_i_174Yb-Gr0.15-Gr1.4-Gr3.7
- `Cut_Shift_In_Category_Run`: 0
- `Cells`: 680466
- `Random_Seed`: not used in this script
- Output files (MD5): `CategoryHistograms.pdf` 9aad0033; `CategoryHistograms.png` 344fb368; `Plotted_data.csv` 65371b49; `Plotted_levels.csv` e273c2b2; `Figure_legend.txt` 5dbfc46b
