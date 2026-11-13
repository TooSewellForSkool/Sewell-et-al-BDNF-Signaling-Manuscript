# How to reproduce a figure panel

1. Install R 4.5.2 and the packages named in the card's `Package_Versions` (data.table, jsonlite, ggplot2, ragg, gridExtra, viridis).
2. Download the repo and the data tables (links go in `DATA.md`, to be written once the files are uploaded). Keep the folder names as released.
3. Open `figure_recipes/cards/<panel>.md` (e.g. `Fig2C.md`). It lists the scripts to run, **in order**, with their versions.
4. For each script: open it in RStudio, set the USER SETTINGS block to the values listed on the card, set the input folder to the previous step's output folder, and run it. The script writes a new dated folder and a `run_parameters.json`.
5. Check your result: the card lists the MD5 of each output file. Compare with `Tools/Make_FileManifest.R`. Table outputs should match exactly. PDF/PNG files can differ by a few bytes between machines; compare `Plotted_data.csv` instead.
6. Cell-level steps (01_PercAboveThreshold onward) use random seeds only where the card lists one.

Panels without a card (schematics, UMAPs, legacy-unchanged panels) are listed in `figure_recipes/Legacy_Figure_Scripts.md` and `Panels_still_to_confirm.md`.
