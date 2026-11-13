# Data for the revision (release asset)

The cell-level input tables are too large for the repository. They are attached to the GitHub release **`v2.0-revision`** as one zip file.

- Download: https://github.com/TooSewellForSkool/Sewell-et-al-BDNF-Signaling-Manuscript/releases/tag/v2.0-revision (under "Assets")
- File: `Sewell_BDNFsignaling_InputData_v2.0-revision.zip`  (701,764,678 bytes, about 670 MB)
- MD5 of the zip: `7f5678123af80302fe53523a18d46dcc`
- Unzipping gives one folder, `InputData/`, holding the files listed below.

680,466 cells in total. One row per cell in every per-cell file, in the same row order.

## What is in the zip

| File | What it is |
|---|---|
| `expression_matrix_analysis (R1).csv` (MD5 `0b2e27d76c1e3f2f2f839b59f40dc24b`) | Round 1 (cell identity) set. Clustering channels (17): Tuj1, CXC3CR1, Vimentin, OligoO4, s100b, HB9, CGRP, NeuN, CHAT, nNOS, Sox2, p75NTR, DCX, MAP2, Islet1, BLBP, TrkB_i |
| `expression_matrix_other (R1).csv` | The remaining 55 channels for round 1 |
| `expression_matrix_analysis (R2).csv` (MD5 `b8e6e5ba3123d895c54da47d819a5bbe`) | Round 2 (signaling) set. Clustering channels (19): pPLK1, H3K9ac, pGSK3b, p-p38, pPLCG-2, pSTAT5, p-p53, pSTAT3, pSrc, pATM, p-cJUN, yH2AX, pNFKB, pAKt, pTrkB, pERK, CC3, pS6, pCreb |
| `expression_matrix_other (R2).csv` | The remaining 53 channels for round 2 |
| `metadata.csv` | Sample information per source file (experiment, treatment, time point, replicate) |
| `filenums.csv` | Source file number of each cell (links a cell to a row of `metadata.csv`) |
| `panel.csv` | Channel (marker) list |
| `cluster_R1_assigns.csv` | Round 1 (cell identity) Leiden cluster number of each cell: 20 clusters |
| `cluster_R1_groups.csv`, `cluster_R2_groups.csv` | Blank templates for naming clusters (cluster numbers only; the label columns are empty). The cell-type groups used in the paper are in `cluster_R1_merge-categories.csv` |
| `cluster_R1_merge.csv` | Rules for merging round 1 clusters (Supplementary Fig. 6), as read by the scripts |
| `cluster_R1_merge-categories.csv` | Readable version of the same rules (reference only; the scripts do not read it) |
| `cluster_R1_assigns_merged.csv` | Round 1 assignments after merging |
| `cluster_R2_assigns.csv`, `cluster_R2_merge.csv`, `cluster_R2_merge-categories.csv`, `cluster_R2_assigns_merged.csv` | The same five files for round 2 (signaling clusters) |
| `umap_xy_R1.csv`, `umap_xy_R2.csv`, `umap_xyz_3D_R1.csv` | UMAP coordinates used for the UMAP figures |
| `Palette_Colors.csv` | One colour per cluster, used across the paper |

## Round 1 and round 2 file names
The two rounds each have their own pair of expression matrices, both covering all 72 channels but split differently into "analysis" and "other". The scripts look for the plain names `expression_matrix_analysis.csv` and `expression_matrix_other.csv`. Before running a script, either:
- copy the pair you need into your input folder and remove the `(R1)` / `(R2)` part of the names, or
- keep the names as released and set `EXPR_ANALYSIS_NAME` and `EXPR_OTHER_NAME` near the top of `02_Analysis/01_PercAboveThreshold` to the full file names.

In each pair, the `analysis` file holds the channels used for that round of clustering and the `other` file holds all remaining channels. Combined, either pair gives the same 72-channel data for every cell (confirmed by the author), so the percent-above-threshold results do not depend on which pair is used. Use the R1 pair with the round 1 cluster files and the R2 pair with the round 2 cluster files for consistency with the clustering.

## How the files are used
The first analysis step, `02_Analysis/01_PercAboveThreshold`, reads the expression matrices, `metadata.csv`, `filenums.csv` and a cluster file. From there, follow the card for the panel you want in `figure_recipes/cards/`.

Unzip into one folder and point the script's input folder at it. Check the zip with the MD5 above before use.
