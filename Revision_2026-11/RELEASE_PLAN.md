# What goes in the release

## In the git repository (small, under 25 MB per file)
- `ScriptsToUse/` (1.9 MB), `LICENSE` (GPL-3.0), `README.md` (from PIPELINE_OVERVIEW.md), `HOW_TO_REPRODUCE.md`, `DATA.md`
- `figure_recipes/` (Figure_Recipes.csv + cards), `run_records/` (one `run_parameters.json` per final panel, flat; folder name = panel)
- Small stage tables (`Expt..._RAWvalues.csv`, `_FC.csv`, `Statistics_*.csv`, `Plotted_data.csv`): a few MB; they let a reader remake every statistic and plot **without** the cell-level data.
- `Tools/GoldenData/`, `Palette_Colors.csv`, `MD5SUMS.txt`

## As release assets (large, up to 2 GiB each, NOT in the repository)
- Cell-level tables (`expression_matrix_*.csv`, 0.18-0.54 GB each; batch-corrected ones 0.19 and 0.53 GB), compressed with zip/gz. Plus `cluster_R*_assigns.csv`, `metadata.csv`, `filenums.csv` (1-2 MB; fine in the repo).

## Not included
Raw FCS files (state where they can be requested or whether they are deposited), personal local paths, exploratory runs with no panel label (79 records).
