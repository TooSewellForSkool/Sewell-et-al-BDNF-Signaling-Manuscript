PALETTE  (Reviewer #3: colour legibility)  -- written 2026-10-28
=================================================================
WHERE THE COLOURS LIVE
  Helpers/Palette_Colors.csv   ONE line per colour: Set (Cluster / Treatment / Category / Tier / MergedGroup / Merge),
                               ID, Label, Color (hex), Note, CoPlotGroup.
  Change a colour THERE; every plotting script reads it. Never type colours into a plotting script.
  Helpers/01_Palette.R         the functions the scripts use (background, automatic contrast adjustment, line types).

HOW IT WORKS
  * CoPlotGroup = clusters that are drawn together in a line graph (1,7,2,13,4 | 12,8,6,11 | 15,3 | 9,17 |
    5,10,14,16,18 | 19 | 20). Colours inside a group are chosen to be told apart, also with colour blindness.
  * Line type inside a group (solid, long dash, dash, dots, dot-dash) is fixed per cluster, so a cluster looks the
    same in every figure and two similar colours can still be separated by their line.
  * BACKGROUND (in each plotting script): "white", "transparent", "published" (grey panel) or any colour.
    A colour whose contrast with the background is below the minimum is lightened (dark background) or darkened
    (light background) automatically. Minimum: 2.0 on light, 3.0 on dark (MIN_RATIO_* in 01_Palette.R).
    Deliberate exception: the light end of the blue gradient (cluster 1, about 2.1:1 on white).
  * "transparent" needs TRANSPARENT_ON = "white" or "black": the colour it will be placed on.

CHECK THE COLOURS
  Source Palette/ColourMap.R (choose BACKGROUND at the top). It writes ColourMap_<background>.png/.pdf (one line per
  colour on the background), ColourMap_CVD_<background>.png (colour-blindness simulation) and
  Palette_Check_<background>.csv (per group: contrast with background, contrast between neighbours, smallest colour
  difference under every vision type).

STATUS: clusters and treatments are wired into 01_PlotLinegraphs.R (v1.10) and 02_CompareLinegraphs.R (v1.3).
  Still to do: 03_PlotBargraphs.R, 04_PlotCategoryHistograms.R, UMAP scripts.
