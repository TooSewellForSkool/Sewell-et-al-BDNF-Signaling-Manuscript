================================================================================
README: Visualization layer          [01_PlotLinegraphs.R (v1)]      (02_CompareLinegraphs.R v1)
================================================================================
Status: 01_PlotLinegraphs v1 written 2026-10-25; not yet run. Golden check pending (see below).

WHAT THIS LAYER DOES
--------------------
Draws figures from tables made in the Analysis layer. It calculates no statistics: the stars are READ from a
Statistics folder (03_Statistics.R or 04_TimecourseStatistics.R). Which p-value draws the stars is a setting.

SCRIPT 1: 01_PlotLinegraphs.R
-----------------------------
Line graphs (one panel per marker): small dots = samples, markers/line = mean, band = +/- SD, stars = p-value below
0.05 / 0.01 / 0.001. The published look is kept (default grey ggplot theme, colours, 5-inch panels, 5 per row).
Input:  STATS_FOLDER = a folder made by 03_Statistics.R (v2+) or 04_TimecourseStatistics.R. The table behind it is found
        automatically (and its checksum checked). 03 must have used STAT_TABLE "FC" or "log2FC".
Choose: P_COLUMN
          03 folder:  "p_value"                     (the published test: Welch t-test, unadjusted)
          04 folder:  "p_unadjusted" | "p_primary" | "q_BH_all_markers"
        MARKERS, TREATMENTS (plotting order), CLUSTERS (c() = pseudobulk), PLOT_BY ("auto" = published rule).
        OUTPUT_FORMATS c("pdf","png") (pdf stays editable in Illustrator), PNG_DPI (300 print, 150 slides).
Output: a new folder <date>_<time>_Linegraphs_By<...>_<p column> INSIDE the Statistics folder with
          Linegraph_<...>-Grid.pdf / .png      the figure
          Figure_legend.txt                    legend text: test, p-value column, adjustment family, what n counts, litters
          Plotted_data.csv                     one row per marker x group x time point as drawn: n, mean, SD, p, stars
          Plotted_points.csv                   every plotted sample value
          Plotted_axes.csv, Plotted_style.csv  y limits and colours (02_CompareLinegraphs.R redraws the panels from these)
          run_parameters.json                  settings, checksums, the legend text, warnings
        CAPTION_ON_FIGURE = TRUE prints a short version of the legend at the bottom of the figure.

Changes from the published script: read the change log (A1-A7) at the top of the script. Two to know:
  * grid rows used round() and silently dropped panels for some marker counts; now ceiling() (A2).
  * the grid of p-value tables is replaced by Plotted_data.csv (A4).
The stars on the published figures were moved by hand in Illustrator; the script reproduces the published
PLACEMENT RULE, not the hand edits.

SETTINGS TO REPRODUCE THE PUBLISHED PSEUDOBULK FIGURE (Fig 2C / Supp Fig 4B style)
  1. 03_Statistics.R (v2): INPUT_FOLDER = the Run A folder, STAT_TABLE "FC", TEST "welch", COMPARE_BY "Treatment",
     TREATMENTS c("NoTreat","ContDepr","BDNF") (add "CompMed" for the 4-treatment version), RUN_DIAGNOSTICS as you like.
  2. 01_PlotLinegraphs.R: STATS_FOLDER = that folder, P_COLUMN "auto", MARKERS = the markers of the figure, CLUSTERS c().
GOLDEN CHECK (planned): Plotted_data.csv vs an independent Python calculation of the same figure
  (Inventory/Reference_PlottedData_python_pseudobulk_NoTreat-ContDepr-BDNF_FC_welch.csv: n, mean, SD, p, stars for all 40
  markers, treatments NoTreat/ContDepr/BDNF), then the published RepsAvg tables and Supplementary Tables S3/S4 for the stars.

SCRIPT 2: 02_CompareLinegraphs.R   (OPTIONAL; nothing depends on it)
-------------------------------------------------------------------
Puts two runs of 01_PlotLinegraphs.R side by side, marker by marker (A left = reference, B right = new version).
Both runs must come from the same table, selection and markers: the script stops if the lines (n, mean, SD) are not
identical, so only the stars differ. Needs the current 01 (it writes Plotted_axes.csv and Plotted_style.csv).
Settings: FOLDER_A, FOLDER_B (output folders of 01), LABEL_A/LABEL_B (default = statistics source + p-value column),
          MARKERS, ONLY_CHANGED_MARKERS (draw only markers where a star changed), PAIRS_PER_ROW, OUTPUT_FORMATS, PNG_DPI.
Output: a folder inside FOLDER_B with Linegraph_Comparison.pdf/.png, Comparison_changes.csv (every group/time point:
        p and stars in A and B, Change = same / gained / lost / stronger / weaker), Comparison_summary.txt, run_parameters.json.
Typical revision use: A = 01 run on the 03 folder (P_COLUMN "p_value", the published test); B = 01 run on the 04 folder
(P_COLUMN "p_primary", or "q_BH_all_markers").
