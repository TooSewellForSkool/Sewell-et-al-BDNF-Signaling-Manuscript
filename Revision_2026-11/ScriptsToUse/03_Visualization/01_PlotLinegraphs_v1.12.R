################################################################################
# 01_PlotLinegraphs.R
################################################################################
# Layer:    Visualization (no statistics are calculated here; p-values are READ from the Analysis layer)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-25
# Edited:   2026-10-31
# Version:  1.12   (NEW; replaces the plotting half of LineGraphAnalysis_... v10-v11)
#
# Purpose:
# * Draw the line graphs of the paper: one panel per marker, mean +/- SD of the samples at each time point,
#   the individual samples as small dots, and stars where a time point differs from the reference.
# * The stars come from the statistics folder you choose (made by 03_Statistics.R or 04_TimecourseStatistics.R),
#   and you choose WHICH p-value draws them (P_COLUMN). The figure legend text that says which test, which
#   p-value and which adjustment was used is written to Figure_legend.txt (and can be printed on the figure).
# * Also writes what was plotted (Plotted_data.csv, Plotted_points.csv) so that two versions of a figure can be
#   compared by 02_CompareLinegraphs.R, and so the plotting can be checked against the published tables.
#
# Look: 5-inch panels, 5 columns, star size and placement are as published. Since v1.10 the COLOURS come from the shared
# palette (Helpers/Palette_Colors.csv) and the background can be chosen (BACKGROUND); see change log A8.
#
# ==============================================================================
# CHANGE LOG: LineGraphAnalysis_... (published) -> v1. v1.1 figure saved with limitsize = FALSE (no effect on the look; avoids refusal of figures taller than 50 in). v1.2 star layer no longer drawn in the legend (removes the large "a" key glyph). v1.3 BH legend states the family size (markers x time points = comparisons). v1.4 MIN_N_FOR_STARS (default 3): no star where a group has fewer than 3 samples; p-values kept; legend names the time points; Plotted_data.csv gains Suppressed_Low_N. v1.5 one expression category per figure (A7). v1.6 figure saved with a WHITE background (bg = "white"): before, the area outside the panels (and the title) was transparent, so viewers with a dark background showed it black; on white paper/slides it looks the same as before. v1.7 the PDF is written through a short temporary path and copied: R's pdf() device cuts paths longer than 511 characters, which saved the PDF as "...Grid-Cat2" (no ending) in deeply nested folders. v1.8 the legend names the measure that was plotted, read from the run record of 01_PercAboveThreshold (MEASURE): "percentage of cells above threshold" or "mean expression of the marker" (Fig 5F); before, it always said percentage above threshold. v1.9 new settings Y_LIMITS and STAR_STEP (default c() / NA = published behaviour, unchanged): for mean-expression fold changes the published axis (-1 to 3) and star spacing (3 units) put the stars outside the plot. v1.10 colours and line types from the shared palette (Helpers/01_Palette.R, Palette_Colors.csv); new settings BACKGROUND and TRANSPARENT_ON; default page white without the grey panel (A8). v1.11 accepts a Statistics folder made from the RAWvalues table (Fig 2C, S4B, S5: percent of all cells in the sample above threshold, not a fold change): y axis label, legend text and lower axis limit (0 instead of minus one SD) follow the table; nothing else changes (A9). v1.12 POINT_SHAPES: the mean points get a shape per cluster (like the line type), so lines stay apart in greyscale and for colour-blind readers; "none" gives the round dots of v1.11 (A10).
# ==============================================================================
# Convention:  [RESULT] = can change what is drawn / the numbers;  [BEHAVIOR] = changes how the script runs
#
# A. Needs your attention (nothing here is applied silently)
#  A1 [RESULT: none intended] Means, SDs, axis limits, colours, star positions, panel layout and titles follow the
#     published script (the version with the stats table, "...timepointComparison__v2__stats_table_.R").
#  A2 [RESULT, FIXED BUG] The published grid used  nrow = round(markers / 5).  With 1-2 markers over a multiple of
#     5 (e.g. 12 markers -> round(2.4) = 2 rows = 10 panels) the last panels were silently NOT drawn. Now
#     ceiling() is used. A 19-marker figure (4 rows) is unaffected.
#  A3 [RESULT, declared] The stars are no longer calculated here: they are read from the Statistics folder, and
#     the thresholds (* < 0.05, ** < 0.01, *** < 0.001) are applied to the p-value column you choose. The
#     published figures correspond to the 03_Statistics.R output ("p_value": Welch's t-test, unadjusted).
#  A4 [BEHAVIOR] The published script also exported a grid of p-value tables (PDF). That is replaced by
#     Plotted_data.csv (it has the p-value of every plotted time point).
#  A5 [RESULT, published behaviour kept] The mean and SD are calculated WITHOUT removing NA values, as before: a
#     group that contains an NA (e.g. an empty cluster in one sample) has no mean and is not drawn.
#  A6 [BEHAVIOR] Axis limits come from the whole table (all treatments and clusters), as in the published script.
#  A7 [RESULT, declared] v1.5: tables with expression categories (ThreshCat) are drawn ONE CATEGORY PER FIGURE: the category
#     is the single one in the Statistics folder (03_Statistics.R CATEGORIES). Rows, statistics and the AXIS LIMITS (A6) then
#     come from that category only (the legacy Fig 6 script is not available to check whether it used the whole table;
#     limits can differ from the published figure by this). The category is in the file and folder names.
#
#  A8 [RESULT: colours/line types only, declared] v1.10: colours and line types come from Helpers/01_Palette.R + Palette_Colors.csv
#     (Reviewer #3, colour legibility). Clusters: colour AND line type depend on the cluster (same in every figure), and
#     clusters that are plotted together get colours that differ also for colour-blind readers. Treatments keep their
#     published colours (NoTreat, ContDepr, CompMed, BDNF). The default page is now white WITHOUT the grey panel
#     (BACKGROUND <- "white"); "published" restores the grey panel, "transparent" gives a figure without a background.
#     Colours too faint for the chosen background are lightened/darkened automatically (listed in the console and in
#     run_parameters.json). Means, SDs, stars, axes and layout are not touched.
#
#  A9 [RESULT, declared] v1.11: STAT_TABLE "RAWvalues" is accepted. The y axis then reads "% of all cells above threshold (marker)"
#     (or "mean expression (marker)" when the table is a mean-expression table) and the legend says the values are not
#     normalised to a baseline. The lower axis limit is 0 (the values cannot be negative) instead of minus one SD; the upper limit,
#     star spacing and everything else follow the same rule as for fold changes. Y_LIMITS and STAR_STEP still override.
#     Reason: the published Fig 2C y axis was labelled "mean fold change" although the plotted values are percentages.
#
#  A10 [RESULT: how the mean points look only, declared] v1.12: new setting POINT_SHAPES (default "auto"). Clusters: the shape of the
#     mean points follows the cluster's position in its co-plotted group, like the line type (round dot, triangle, square, diamond,
#     cross; Helpers/01_Palette.R). Treatments: round dots as before (optional "Shape" column in Palette_Colors.csv). "none" = all round
#     dots (identical to v1.11); or give your own R point codes, one per group in plotting order, e.g. c(19, 17, 15). The small sample
#     dots stay round, and the 0 hr point shared by all lines stays the round NoTreat dot. Means, SDs, stars, axes and layout are not
#     touched; the shapes are written to Plotted_style.csv, run_parameters.json and the legend.
#
# B. Efficiency / readability / reproducibility
#  B1 One run record (run_parameters.json): settings, the statistics folder and its checksums, outputs, session.
#  B2 PDF (text and points stay editable in Illustrator; dots are real circles, not font symbols) and PNG in one go.
#  B3 No random numbers are used.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
STATS_FOLDER   <- ""            # a folder made by 03_Statistics.R (v2+) or 04_TimecourseStatistics.R ("" = choose in a pop-up)

# Which p-value draws the stars?
#   "auto"  = "p_value" for a 03_Statistics.R folder (published test) / "p_primary" for a 04 folder
#   03 folder: "p_value"
#   04 folder: "p_unadjusted", "p_primary" or "q_BH_all_markers"  (see the header of 04_TimecourseStatistics.R)
P_COLUMN       <- "p_primary"

# What is drawn (c() = everything the statistics folder has)
MARKERS        <-   c("pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La")
  # c("pERK_171Yb", "pCreb_176Lu", "pS6_175Lu", "p-p38_142Nd", "pGSK3b_139La")
  # c("pERK_171Yb", "pCreb_176Lu", "pAKt_168Er", "pPLK1_113In", "H3K9ac_115In",
  #                   "pSTAT3_153Eu", "pSTAT5_147Sm","pATM_159Tb","yH2AX_165Ho","p-p53_151Eu",
  #                   "pPLCG-2_144Nd", "pSrc_154Sm","p-cJUN_162Dy", "pNFKB_166Er","CC3_173Yb")
# c("pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La")
# c("pPLK1_113In", "yH2AX_165Ho", "TrkB_s_158Gd",
#   "pATM_159Tb", "pPLCG-2_144Nd", "pSrc_154Sm",
#   "H3K9ac_115In", "p-p38_142Nd", "CC3_173Yb",
#   "pSTAT3_153Eu", "pSTAT5_147Sm", "p-p53_151Eu",
#   "p-cJUN_162Dy", "pNFKB_166Er", "pTrkB_170Er")

TREATMENTS     <- c("NoTreat", "ContDepr", "BDNF")           # in plotting order, e.g. c("NoTreat", "ContDepr", "CompMed", "BDNF"); c() = the treatments of the Statistics folder
CLUSTERS       <- c()           # e.g. c(1, 2, 4, 18); c() = pseudobulk
PLOT_BY        <- "auto"        # "auto" (published rule), "Treatment" or "Cluster"
                                #   Treatment = one colour per treatment;  Cluster = one colour per cluster
                                #   (needs exactly one treatment besides the reference)

# How the spread and the stars are drawn
SD_STYLE        <- "ribbon"          # "ribbon" (published: shaded mean +/- SD band) or "errorbar" (SD bars at each time point)
ERRORBAR_MIN_N  <- 3                 # "errorbar" only: SD bars are drawn where a group has at least this many samples
                                     #   (an SD from 2 samples is not drawn; the mean is still shown)
MIN_N_FOR_STARS <- 3                 # no star is drawn for a group with fewer than this many samples (0 = draw all stars).
                                     #   Default 3 (decision 2026-10-26: no statistics are reported on n = 2).
                                     #   The p-values stay in Plotted_data.csv and the Statistics tables; only the star is withheld.
Y_LIMITS        <- c()               # c() = published rule (axis limits from the whole table; star rows spaced by 3 units).
                                     #   Fits percent-above-threshold values. For MEAN-EXPRESSION fold changes (values around 1) the
                                     #   published rule gives an axis of about -1 to 3 and the stars fall outside it. Then give the limits,
                                     #   e.g. c(0.4, 1.6): the stars are stacked from the top (one row per cluster/treatment).
STAR_STEP       <- NA                # distance between star rows when Y_LIMITS is set (NA = 5 % of the axis range)
NO_STARS_TREATMENTS <- c()           # treatments drawn WITHOUT stars even if the Statistics folder has them,
                                     #   e.g. c("CompMed") for a supplementary panel that shows no statistics
                                     #   (only for PLOT_BY "Treatment"; the legend says so)

# Colours and background (colours are defined in Helpers/Palette_Colors.csv)
POINT_SHAPES    <- "auto"            # shape of the MEAN points: "auto" (by cluster, like the line type; treatments stay round dots),
                                     #   "none" (all round dots, as before v1.12) or your own R point codes, one per group in plotting order,
                                     #   e.g. c(19, 17, 15) (19 = circle, 17 = triangle, 15 = square, 18 = diamond, 4 = cross)
BACKGROUND      <- "white"           # "white", "transparent", "published" (grey panel of the original figures) or any colour, e.g. "#1E1E1E"
TRANSPARENT_ON  <- "white"           # only for BACKGROUND "transparent": "white" or "black" = what the figure will be placed on
                                     #   (decides the contrast check and whether text and axes are drawn black or white)

# Output
OUTPUT_FORMATS  <- c("pdf", "png")   # "pdf" (editable in Illustrator), "png" (presentations)
PNG_DPI         <- 300               # 300 for print; 150 is enough for slides (and a smaller file)
PANEL_COLUMNS   <- 5                 # panels per row
PANEL_SIZE_IN   <- 5                 # width and height of one panel (inches)
CAPTION_ON_FIGURE <- FALSE           # TRUE: print a short version of the legend at the bottom of the figure
OUTPUT_LABEL    <- ""                # optional text added to the output folder name

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "01_PlotLinegraphs.R"
SCRIPT_VERSION           <- "1.12"
MODULE_BASE_NAME         <- "01_PlotLinegraphs"
REQUIRED_HELPERS_VERSION <- 7
REQUIRED_PALETTE_VERSION <- 3
ID_COLUMNS <- c("File_Num", "File_Name", "Sample_Condition", "Reordered", "TimePoint",
                "Replicate", "Cluster", "ThreshCat", "N_Cells", "N_Replicates")



#### Load shared helper functions (do not edit this block) #####################
local({
  start_dir <- tryCatch({
    p <- rstudioapi::getSourceEditorContext()$path
    if (is.character(p) && nzchar(p)) dirname(p) else getwd()
  }, error = function(e) getwd())
  search_dir   <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  helpers_file <- NA_character_
  for (i in 1:5) {
    candidates <- c(file.path(search_dir, "Helpers", "00_Helpers.R"),
                    file.path(search_dir, "00_Helpers.R"))
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates) > 0) { helpers_file <- candidates[1]; break }
    search_dir <- dirname(search_dir)
  }
  if (is.na(helpers_file)) {
    stop("Could not find 00_Helpers.R.\n",
         "It should be in a folder called 'Helpers' inside ScriptsToUse.\n",
         "Looked upward from: ", start_dir,
         "\n(Is this script saved inside your ScriptsToUse folder, and is the ",
         "'rstudioapi' package installed? RStudio needs it to find this script.)",
         call. = FALSE)
  }
  assign("HELPERS_DIR", dirname(helpers_file), envir = .GlobalEnv)
  source(helpers_file, local = FALSE)
})
check_helpers_version(REQUIRED_HELPERS_VERSION)
source(file.path(HELPERS_DIR, "01_Palette.R"), local = FALSE)
if (PALETTE_HELPERS_VERSION < REQUIRED_PALETTE_VERSION) stop("01_Palette.R is too old. Copy the newest one into: ", HELPERS_DIR, call. = FALSE)



#### Set up environment ########################################################
pkg_versions <- dc_load_packages(cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite",
                                          "ggplot2", "gridExtra", "viridis"))
has_ragg <- requireNamespace("ragg", quietly = TRUE)

cat("\n========================================================\n")
cat("[START] Line graphs with stars from the Statistics folder:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

warnings_log <- character(0)
note_warning <- function(txt) {
  warnings_log <<- c(warnings_log, txt)
  cat("[WARNING] ", txt, "\n", sep = "")
}
one_of <- function(value, allowed, name) {
  if (!is.character(value) || length(value) != 1 || !value %in% allowed) {
    stop(name, " must be one of: ", paste0("\"", allowed, "\"", collapse = ", "), ".", call. = FALSE)
  }
}
one_of(PLOT_BY, c("auto", "Treatment", "Cluster"), "PLOT_BY")
one_of(SD_STYLE, c("ribbon", "errorbar"), "SD_STYLE")
one_of(TRANSPARENT_ON, c("white", "black"), "TRANSPARENT_ON")
bg <- palette_background(BACKGROUND, TRANSPARENT_ON)
if (!is.numeric(MIN_N_FOR_STARS) || length(MIN_N_FOR_STARS) != 1 || MIN_N_FOR_STARS < 0) {
  stop("MIN_N_FOR_STARS must be one number, 0 or more (0 = draw all stars).", call. = FALSE)
}
if (!is.numeric(ERRORBAR_MIN_N) || length(ERRORBAR_MIN_N) != 1 || ERRORBAR_MIN_N < 2) {
  stop("ERRORBAR_MIN_N must be one number, 2 or more.", call. = FALSE)
}
one_of(P_COLUMN, c("auto", "p_value", "p_unadjusted", "p_primary", "q_BH_all_markers"), "P_COLUMN")
if (length(OUTPUT_FORMATS) == 0 || !all(OUTPUT_FORMATS %in% c("pdf", "png"))) {
  stop("OUTPUT_FORMATS must be \"pdf\", \"png\" or both.", call. = FALSE)
}
for (nm in c("PNG_DPI", "PANEL_COLUMNS", "PANEL_SIZE_IN")) {
  v <- get(nm)
  if (!is.numeric(v) || length(v) != 1 || v <= 0) stop(nm, " must be one positive number.", call. = FALSE)
}
if (!is.logical(CAPTION_ON_FIGURE) || length(CAPTION_ON_FIGURE) != 1 || is.na(CAPTION_ON_FIGURE)) {
  stop("CAPTION_ON_FIGURE must be TRUE or FALSE.", call. = FALSE)
}
if ("png" %in% OUTPUT_FORMATS && !has_ragg) {
  note_warning("Package 'ragg' is not installed: PNG files are made with R's standard png device (works, slightly less crisp text). install.packages(\"ragg\") improves it.")
}



#### Choose the Statistics folder ##############################################
is_stat_folder <- function(d) (file.exists(file.path(d, "Statistics_long.csv")) ||
                               file.exists(file.path(d, "Timecourse_long.csv"))) &&
                              file.exists(file.path(d, RUN_RECORD_NAME))
if (nzchar(STATS_FOLDER)) {
  if (!dir.exists(STATS_FOLDER)) stop("STATS_FOLDER not found:\n  ", STATS_FOLDER, call. = FALSE)
  stat_dir <- normalizePath(STATS_FOLDER, winslash = "/", mustWork = TRUE)
} else {
  Original_WD <- get_working_dir("project")
  dirs <- list.dirs(Original_WD, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[!grepl("/[_.][^/]*", sub(Original_WD, "", dirs, fixed = TRUE))]
  dirs <- sort(dirs[vapply(dirs, is_stat_folder, logical(1))], decreasing = TRUE)
  if (is_stat_folder(Original_WD)) dirs <- c(Original_WD, dirs)
  dirs <- unique(dirs)
  if (length(dirs) == 0) {
    stop("No Statistics folder (with Statistics_long.csv or Timecourse_long.csv) was found in:\n  ", Original_WD,
         "\nRun 03_Statistics.R (or 04_TimecourseStatistics.R) first, or set STATS_FOLDER.", call. = FALSE)
  }
  labels <- basename(dirs)
  pick <- ask_choice("Which Statistics folder should be plotted?", labels, preselect = labels[1])
  stat_dir <- normalizePath(dirs[match(pick, labels)], winslash = "/", mustWork = TRUE)
}
cat("[INFO] Statistics folder:\n       ", stat_dir, "\n", sep = "")

is_timecourse <- file.exists(file.path(stat_dir, "Timecourse_long.csv"))
stat_file <- file.path(stat_dir, if (is_timecourse) "Timecourse_long.csv" else "Statistics_long.csv")
record <- require_modules(stat_dir, "03_Statistics", how_to_fix = "Run 03_Statistics.R first.")
s3 <- record[["03_Statistics"]]
s4 <- record[["04_TimecourseStatistics"]]
if (is_timecourse && is.null(s4)) stop("Timecourse_long.csv is here, but the run record has no 04_TimecourseStatistics section.", call. = FALSE)
as_chr <- function(v) if (is.null(v) || length(v) == 0) character(0) else as.character(unlist(v))

# p-value column
if (P_COLUMN == "auto") P_COLUMN <- if (is_timecourse) "p_primary" else "p_value"
allowed_p <- if (is_timecourse) c("p_unadjusted", "p_primary", "q_BH_all_markers") else "p_value"
if (!P_COLUMN %in% allowed_p) {
  stop("P_COLUMN = \"", P_COLUMN, "\" is not available in this Statistics folder (made by ",
       if (is_timecourse) "04_TimecourseStatistics.R" else "03_Statistics.R", ").\nChoose one of: ",
       paste(allowed_p, collapse = ", "), call. = FALSE)
}
cat("[INFO] Stars are drawn from the column: ", P_COLUMN, "\n", sep = "")

# The table behind the Statistics folder (made by 01_PercAboveThreshold.R)
table_dir  <- if (is_timecourse) dirname(dirname(stat_dir)) else dirname(stat_dir)
table_path <- file.path(table_dir, s3$Input_Table)
if (!file.exists(table_path)) {
  stop("The table behind the Statistics folder was not found:\n  ", table_path,
       "\nKeep the Statistics folder inside the folder of the table (and the 04 folder inside the 03 folder).", call. = FALSE)
}
if (!identical(unname(as.character(tools::md5sum(table_path))), as.character(s3$Input_Table_MD5))) {
  stop("The table has changed since the statistics were calculated (checksum differs):\n  ", table_path,
       "\nRun 03_Statistics.R (and 04) again.", call. = FALSE)
}
STAT_TABLE          <- s3$Settings$STAT_TABLE
REFERENCE_TREATMENT <- s3$Settings$REFERENCE_TREATMENT
REFERENCE_TIMEPOINT <- as.numeric(s3$Settings$REFERENCE_TIMEPOINT)
if (!STAT_TABLE %in% c("FC", "log2FC", "RAWvalues")) {
  stop("Line graphs need STAT_TABLE \"FC\", \"log2FC\" or \"RAWvalues\" in 03_Statistics.R; this folder used \"",
       STAT_TABLE, "\".", call. = FALSE)
}
measure_txt <- if (identical(as.character(record[["01_PercAboveThreshold"]]$Settings$MEASURE), "mean_expression"))
  "mean expression of the marker (values as stored, arcsinh scale)" else "percentage of cells above threshold"
cat("[INFO] Table: ", basename(table_path), "  (", STAT_TABLE, ")\n", sep = "")



#### Read the table ############################################################
tb <- data.table::fread(table_path, na.strings = c("", "NA"), data.table = TRUE)
if (names(tb)[1] %in% c("", "V1") && ncol(tb) > 1 &&
    identical(as.integer(tb[[1]]), seq_len(nrow(tb)))) tb[, (1) := NULL]
for (cn in c("Sample_Condition", "TimePoint", "Cluster", "Replicate")) {
  if (!cn %in% names(tb)) stop("The table has no column '", cn, "'.", call. = FALSE)
}
PLOT_CATEGORY <- NULL                                   # set below when the table has expression categories (A7)
if ("ThreshCat" %in% names(tb) && any(!is.na(tb$ThreshCat) & tb$ThreshCat != "")) {
  stat_cats <- as.character(unlist(s3$Settings$CATEGORIES))
  if (length(stat_cats) != 1) {
    stop("This table has expression categories (ThreshCat). Run 03_Statistics.R with exactly ONE category in CATEGORIES ",
         "(one figure per category); this Statistics folder has ", length(stat_cats), ".", call. = FALSE)
  }
  PLOT_CATEGORY <- stat_cats
  tb[, ThreshCat := as.character(ThreshCat)]
  tb <- tb[ThreshCat == PLOT_CATEGORY]
  if (nrow(tb) == 0) stop("No rows of category ", PLOT_CATEGORY, " in the table.", call. = FALSE)
  cat("[INFO] Expression category ", PLOT_CATEGORY, ": ", nrow(tb), " table rows (axis limits come from this category only)\n", sep = "")
}
tb[, Treatment := sub("_[^_]+$", "", Sample_Condition)]
tb[, Cluster := as.character(Cluster)]
all_markers <- setdiff(names(tb), c(ID_COLUMNS, "Treatment"))

# Axis limits from the WHOLE table (published behaviour, A6)
y_max <- vapply(all_markers, function(m) max(tb[[m]][is.finite(tb[[m]])], na.rm = TRUE), numeric(1))
y_min <- vapply(all_markers, function(m) if (STAT_TABLE == "RAWvalues") 0 else ceiling(stats::sd(tb[[m]], na.rm = TRUE)) * -1, numeric(1))   # A9: raw values start at 0

# Statistics
st <- data.table::fread(stat_file, na.strings = c("", "NA"), data.table = TRUE)
if (!P_COLUMN %in% names(st)) stop("Column ", P_COLUMN, " is not in ", basename(stat_file), ".", call. = FALSE)
st[, Cluster := as.character(Cluster)]
st[, TimePoint := as.numeric(TimePoint)]
if (!is.null(PLOT_CATEGORY) && "ThreshCat" %in% names(st)) {
  st <- st[as.character(ThreshCat) == PLOT_CATEGORY]
  if (nrow(st) == 0) stop("The Statistics file has no rows for category ", PLOT_CATEGORY, ".", call. = FALSE)
}



#### What is plotted ###########################################################
stats_treatments <- as_chr(s3$Settings$TREATMENTS)        # the treatments the statistics were calculated for
treatments_plot <- if (length(TREATMENTS) > 0) TREATMENTS else if (length(stats_treatments) > 0) {
  unique(c(REFERENCE_TREATMENT, stats_treatments))        # default: only what has statistics
} else unique(c(REFERENCE_TREATMENT, setdiff(unique(tb$Treatment), REFERENCE_TREATMENT)))
if (length(NO_STARS_TREATMENTS) > 0 && length(setdiff(NO_STARS_TREATMENTS, treatments_plot)) > 0) {
  stop("NO_STARS_TREATMENTS not among the plotted treatments: ",
       paste(setdiff(NO_STARS_TREATMENTS, treatments_plot), collapse = ", "), call. = FALSE)
}
if (REFERENCE_TREATMENT %in% NO_STARS_TREATMENTS) {
  stop("NO_STARS_TREATMENTS cannot contain the reference treatment (", REFERENCE_TREATMENT, ").", call. = FALSE)
}
unstarred <- setdiff(treatments_plot, REFERENCE_TREATMENT)
if (length(stats_treatments) > 0) unstarred <- setdiff(unstarred, stats_treatments)
unstarred <- union(unstarred, NO_STARS_TREATMENTS)        # drawn without stars (no statistics, or by your choice)
if (length(setdiff(unstarred, NO_STARS_TREATMENTS)) > 0) {
  note_warning(paste0("Treatment(s) drawn without statistics (the Statistics folder covers ",
                      paste(stats_treatments, collapse = ", "), "): ",
                      paste(setdiff(unstarred, NO_STARS_TREATMENTS), collapse = ", "),
                      ". Add them to NO_STARS_TREATMENTS to confirm this is intended, or re-run 03_Statistics.R with them in TREATMENTS."))
}
bad_t <- setdiff(treatments_plot, unique(tb$Treatment))
if (length(bad_t) > 0) {
  stop("TREATMENTS not found in the table: ", paste(bad_t, collapse = ", "),
       "\nTreatments in the table: ", paste(unique(tb$Treatment), collapse = ", "), call. = FALSE)
}
if (!REFERENCE_TREATMENT %in% treatments_plot) {
  stop("TREATMENTS must include the reference treatment (", REFERENCE_TREATMENT, ").", call. = FALSE)
}
clusters_plot <- as.character(CLUSTERS)
if (length(clusters_plot) > 0 && length(setdiff(clusters_plot, unique(tb$Cluster))) > 0) {
  stop("CLUSTERS not found in the table: ", paste(setdiff(clusters_plot, unique(tb$Cluster)), collapse = ", "), call. = FALSE)
}
markers_plot <- if (length(MARKERS) > 0) MARKERS else intersect(all_markers, unique(st$Marker))
if (length(setdiff(markers_plot, all_markers)) > 0) {
  stop("MARKERS not found in the table: ", paste(setdiff(markers_plot, all_markers), collapse = ", "), call. = FALSE)
}
no_stats <- setdiff(markers_plot, unique(st$Marker))
if (length(no_stats) > 0) {
  note_warning(paste0(length(no_stats), " marker(s) are not in the Statistics folder and are drawn without stars: ",
                      paste(utils::head(no_stats, 5), collapse = ", "), if (length(no_stats) > 5) ", ..." else ""))
}

# Published rule for PlotBy
plot_by <- PLOT_BY
if (plot_by == "auto") {
  plot_by <- if (length(clusters_plot) == 0) "Treatment" else if (length(treatments_plot) > 3) "Treatment" else "Cluster"
}
non_ref <- setdiff(treatments_plot, REFERENCE_TREATMENT)
if (plot_by == "Cluster") {
  if (length(clusters_plot) == 0) stop("PLOT_BY = \"Cluster\" needs CLUSTERS.", call. = FALSE)
  if (length(non_ref) != 1) {
    stop("PLOT_BY = \"Cluster\" needs exactly one treatment besides ", REFERENCE_TREATMENT,
         " in TREATMENTS (found ", length(non_ref), ").", call. = FALSE)
  }
}
if (plot_by == "Treatment" && length(clusters_plot) > 1) {
  stop("PLOT_BY = \"Treatment\" can draw only one cluster (or pseudobulk) at a time; CLUSTERS has ",
       length(clusters_plot), ".", call. = FALSE)
}
is_pseudobulk <- length(clusters_plot) == 0
cat("[INFO] Plotting by ", plot_by, ": ", length(markers_plot), " markers; treatments ",
    paste(treatments_plot, collapse = ", "),
    if (is_pseudobulk) "; pseudobulk" else paste0("; clusters ", paste(clusters_plot, collapse = ",")), "\n", sep = "")

# The published statistics must match how the groups are drawn
if (!is_timecourse) {
  cmp_by <- unique(st$Compare_By)
  if (length(cmp_by) != 1 || cmp_by != plot_by) {
    stop("This 03_Statistics.R folder compared by '", paste(cmp_by, collapse = ","), "' but the figure is plotted by '",
         plot_by, "'. Use a Statistics folder made with COMPARE_BY = \"", plot_by, "\".", call. = FALSE)
  }
}

# Selected data
tb_sel <- tb[Treatment %in% treatments_plot]
if (!is_pseudobulk) tb_sel <- tb_sel[Cluster %in% clusters_plot]
if (nrow(tb_sel) == 0) stop("The selection holds 0 table rows.", call. = FALSE)
tp_levels <- as.character(sort(unique(as.numeric(tb_sel$TimePoint))))
group_levels <- if (plot_by == "Treatment") treatments_plot else clusters_plot

# Colours and line types (shared palette; see Helpers/01_Palette.R)
clusters_merged <- tolower(as.character(unlist(record[["01_PercAboveThreshold"]][["Settings"]][["CLUSTERS_MERGED"]]))[1]) %in% c("true", "t", "yes")
make_colors <- function() {
  n <- length(group_levels)
  if (plot_by == "Cluster") {
    res <- palette_clusters(clusters_plot, bg, merged = clusters_merged)
    if (length(res$missing) > 0) note_warning(paste0("No palette colour for cluster(s) ", paste(res$missing, collapse = ", "),
                                                      " (the palette has clusters 1-20): drawn in grey. Add them to Palette_Colors.csv."))
  } else {
    res <- palette_treatments(group_levels, bg)
    if (is.null(res)) {                           # a treatment that is not in the palette: the published rule
      if (n > 3) stop("More than 3 treatments: colours are defined for NoTreat, ContDepr, CompMed, BDNF and BDNF+K252aTest (Palette_Colors.csv).", call. = FALSE)
      cols <- c("lightblue4", viridis::viridis(n = n, option = "A")[-3])[seq_len(n)]
      res <- list(colors = stats::setNames(cols, group_levels), linetypes = stats::setNames(rep("solid", n), group_levels),
                  shapes = stats::setNames(rep(PALETTE_SHAPES[1], n), group_levels),
                  lifted = character(0), missing = character(0))
      note_warning("A treatment is not in Palette_Colors.csv: the published colour rule was used for this figure.")
    }
  }
  if (length(res$lifted) > 0) cat("[INFO] Colour(s) of ", paste(res$lifted, collapse = ", "), " adjusted to reach the minimum contrast with the background (", bg$name, ").\n", sep = "")
  res$colors <- stats::setNames(res$colors[as.character(group_levels)], group_levels)
  res$linetypes <- stats::setNames(res$linetypes[as.character(group_levels)], group_levels)
  res$shapes <- stats::setNames(res$shapes[as.character(group_levels)], group_levels)
  if (is.character(POINT_SHAPES) && length(POINT_SHAPES) == 1 && POINT_SHAPES == "none") {
    res$shapes[] <- PALETTE_SHAPES[1]
  } else if (!(is.character(POINT_SHAPES) && length(POINT_SHAPES) == 1 && POINT_SHAPES == "auto")) {
    res$shapes <- palette_shapes_custom(POINT_SHAPES, group_levels)
  }
  res
}
palette_res <- make_colors()
group_colors <- palette_res$colors
group_linetypes <- palette_res$linetypes
group_shapes <- palette_res$shapes
shapes_vary <- length(unique(group_shapes)) > 1

# Stars
star_of <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
get_stats_for <- function(m) {
  s <- st[Marker == m]
  if (nrow(s) == 0) return(NULL)
  if (!is_timecourse) {                          # 03: Group = treatment or cluster
    s <- s[, .(Group = as.character(Group), TimePoint, N_Stats = N_Group, p = get(P_COLUMN))]
  } else if (plot_by == "Treatment") {           # 04: one row per treatment x (cluster) x time point
    if (!is_pseudobulk) s <- s[Cluster == clusters_plot[1]]
    s <- s[, .(Group = Treatment, TimePoint, N_Stats = N_Group, p = get(P_COLUMN))]
  } else {                                       # 04, by cluster: the one treatment besides the reference
    s <- s[Treatment == non_ref & Cluster %in% clusters_plot]
    s <- s[, .(Group = Cluster, TimePoint, N_Stats = N_Group, p = get(P_COLUMN))]
  }
  s <- s[Group %in% group_levels]
  if (plot_by == "Treatment" && length(NO_STARS_TREATMENTS) > 0) s <- s[!Group %in% NO_STARS_TREATMENTS]
  s[, Stars := star_of(p)]
  s[, Suppressed_Low_N := !is.na(p) & !is.na(N_Stats) & N_Stats < MIN_N_FOR_STARS & MIN_N_FOR_STARS > 0]
  s[Suppressed_Low_N == TRUE, Stars := ""]
  s
}



#### Build the panels ##########################################################
y_label <- if (STAT_TABLE == "FC") "Fold Change (%s / 0hr)" else if (STAT_TABLE == "RAWvalues") {
  if (identical(as.character(record[["01_PercAboveThreshold"]]$Settings$MEASURE), "mean_expression")) "Mean expression (%s)"
  else "%% of all cells above threshold (%s)"
} else "log2FC(%s / 0hr)"
panel_title_suffix <- if (is_pseudobulk) "Pseudobulk" else paste0("Cluster ", paste(clusters_plot, collapse = ","))
plotted_data <- list(); plotted_points <- list(); axes_list <- list(); panels <- list(); n_mismatch <- 0L; n_nostat <- 0L

for (m in markers_plot) {
  pts <- tb_sel[, .(Sample_Condition, Treatment, Cluster, TimePoint, Replicate, Value = get(m))]
  pts[, Group := if (plot_by == "Treatment") Treatment else Cluster]
  summ <- pts[, .(N = sum(!is.na(Value)), N_Rows = .N, Mean = mean(Value), SD = stats::sd(Value)),
              by = .(Treatment, Cluster, TimePoint, Group)]
  summ[, Is_Appended_0hr := FALSE]
  with0    <- unique(summ[abs(TimePoint) < 1e-9, Treatment])
  without0 <- setdiff(unique(summ$Treatment), with0)
  if (length(without0) > 0) {
    zero_rows <- summ[abs(TimePoint) < 1e-9]
    appended <- data.table::rbindlist(lapply(without0, function(trt) {
      z <- data.table::copy(zero_rows); z[, Treatment := trt]
      z[, Group := if (plot_by == "Treatment") trt else Group]
      z[, Is_Appended_0hr := TRUE]; z
    }))
    summ <- data.table::rbindlist(list(summ, appended), use.names = TRUE)
  }
  summ[, Group := factor(Group, levels = group_levels)]
  summ[, TimePoint_f := factor(as.character(TimePoint), levels = tp_levels)]
  pts[, Group := factor(Group, levels = group_levels)]
  pts[, TimePoint_f := factor(as.character(TimePoint), levels = tp_levels)]

  # stars
  s <- get_stats_for(m)
  y_upper <- y_max[[m]] + abs(y_min[[m]])
  y_lower <- y_min[[m]]
  y_top <- y_upper * 1.12; star_base <- y_upper * 1.06; star_step <- 3        # published rule
  if (length(Y_LIMITS) == 2) {
    y_lower <- Y_LIMITS[1]; y_top <- Y_LIMITS[2]; y_upper <- y_top
    star_step <- if (is.na(STAR_STEP)) 0.05 * (y_top - y_lower) else STAR_STEP
    star_base <- y_top - 0.5 * star_step
  }
  axes_list[[m]] <- data.table::data.table(Marker = m, y_lower = y_lower, y_upper = y_upper, Y_Label = sprintf(y_label, m))
  stars_df <- NULL
  if (!is.null(s)) {
    s[, TimePoint_f := factor(as.character(TimePoint), levels = tp_levels)]
    s[, Group := factor(Group, levels = group_levels)]
    s[, y := star_base - (star_step * (as.integer(Group) - 1))]
    stars_df <- as.data.frame(s[Stars != "" & !is.na(TimePoint_f)])
  }

  # record what was plotted (one row per group x time point, as drawn: appended 0 hr copies are flagged)
  pd <- data.table::copy(summ)
  pd[, Marker := m]
  pd[, TimePoint := as.numeric(as.character(TimePoint_f))]
  if (!is.null(s)) {
    pd <- merge(pd, s[, .(Group, TimePoint = as.numeric(as.character(TimePoint_f)), N_Stats, p, Stars, Suppressed_Low_N)],
                by = c("Group", "TimePoint"), all.x = TRUE, sort = FALSE)
    n_mismatch <- n_mismatch + sum(!is.na(pd$N_Stats) & pd$N_Stats != pd$N & !pd$Is_Appended_0hr)
  } else {
    pd[, `:=`(N_Stats = NA_integer_, p = NA_real_, Stars = NA_character_, Suppressed_Low_N = NA)]
  }
  pd[, `:=`(P_Source = basename(stat_dir), P_Column = P_COLUMN)]
  plotted_data[[m]] <- pd[, .(Marker, PlotBy = plot_by, Group = as.character(Group), Treatment, Cluster, TimePoint,
                              N, N_Rows, Mean, SD, Is_Appended_0hr, N_Stats, P_Source, P_Column, p, Stars, Suppressed_Low_N)]
  pp <- pts[, .(Marker = m, Group = as.character(Group), Treatment, Cluster, TimePoint, Replicate, Sample_Condition, Value)]
  plotted_points[[m]] <- pp

  summ_df <- as.data.frame(summ); pts_df <- as.data.frame(pts)
  summ_with0 <- summ_df[summ_df$Treatment %in% with0, ]
  summ_without0 <- summ_df[summ_df$Treatment %in% without0, ]
  p <- ggplot2::ggplot(pts_df, ggplot2::aes(x = TimePoint_f, y = Value, color = Group)) +
    ggplot2::geom_point(size = 1, alpha = 0.5, show.legend = !shapes_vary) +
    ggplot2::geom_point(data = if (shapes_vary) summ_without0[!summ_without0$Is_Appended_0hr, ] else summ_without0,
                        ggplot2::aes(y = Mean, shape = Group), size = 3) +
    ggplot2::geom_line(data = summ_without0, ggplot2::aes(y = Mean, group = Group, linetype = Group), linewidth = 1) +
    ggplot2::geom_point(data = summ_with0, ggplot2::aes(y = Mean, shape = Group), size = 3) +
    ggplot2::labs(title = paste0(m, ": ", panel_title_suffix), color = plot_by, fill = plot_by, linetype = plot_by, shape = plot_by,
                  x = "Time (hr)", y = sprintf(y_label, m)) +
    ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = group_linetypes, drop = FALSE) +
    ggplot2::scale_shape_manual(values = group_shapes, drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(y_lower, y_top)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    palette_theme(bg)
  if (SD_STYLE == "ribbon") {
    # inserted as the 2nd layer (after the sample dots), the order of the published script
    p$layers <- append(p$layers,
                       list(ggplot2::geom_ribbon(data = summ_df,
                                                 ggplot2::aes(x = TimePoint_f, y = Mean, ymin = Mean - SD, ymax = Mean + SD,
                                                              fill = Group, group = Group),
                                                 alpha = 0.2, color = NA, inherit.aes = FALSE)),
                       after = 1)
  } else {
    bars <- summ_df[!is.na(summ_df$SD) & summ_df$N >= ERRORBAR_MIN_N & !summ_df$Is_Appended_0hr, ]
    if (nrow(bars) > 0) {
      p <- p + ggplot2::geom_errorbar(data = bars,
                                      ggplot2::aes(x = TimePoint_f, ymin = Mean - SD, ymax = Mean + SD, color = Group),
                                      width = 0.2, inherit.aes = FALSE)
    }
  }
  if (!is.null(stars_df) && nrow(stars_df) > 0) {
    p <- p + ggplot2::geom_text(data = stars_df,
                                ggplot2::aes(x = TimePoint_f, y = y, label = Stars, color = Group),
                                inherit.aes = FALSE, size = 10, show.legend = FALSE)
  } else if (is.null(s)) {
    n_nostat <- n_nostat + 1L
  }
  panels[[m]] <- p
}
plotted_data   <- data.table::rbindlist(plotted_data)
plotted_points <- data.table::rbindlist(plotted_points)
if (n_mismatch > 0) {
  note_warning(paste0(n_mismatch, " plotted group/time point(s) have a different n than the statistics used (column N vs N_Stats in ",
                      "Plotted_data.csv). Check that the Statistics folder was made from the same selection."))
}
n_stars <- sum(plotted_data$Stars != "" & !plotted_data$Is_Appended_0hr, na.rm = TRUE)
cat("[INFO] ", length(panels), " panels; ", n_stars, " starred group/time points\n", sep = "")



#### Figure legend text ########################################################
fmt_range <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) "NA" else if (min(x) == max(x)) as.character(min(x)) else paste0(min(x), "-", max(x)) }
n_by_group <- plotted_data[Is_Appended_0hr == FALSE, .(n = fmt_range(N)), by = Group]
n_txt <- paste0("Samples per time point: ", fmt_range(plotted_data[Is_Appended_0hr == FALSE, N]),
                " (", paste0(n_by_group$Group, ": ", n_by_group$n, collapse = "; "), ").")
ref_txt <- paste0(REFERENCE_TREATMENT, " at ", REFERENCE_TIMEPOINT, " hr")
design_file <- file.path(if (is_timecourse) dirname(stat_dir) else stat_dir, "Statistics_design.csv")
litter_txt <- ""
if (file.exists(design_file)) {
  dz <- data.table::fread(design_file)
  if (all(c("N_Samples", "N_Litters") %in% names(dz))) {
    litter_txt <- paste0(" n counts independent culture wells (samples); the ", fmt_range(dz$N_Samples), " samples per group come from ",
                         fmt_range(dz$N_Litters), " litters (Statistics_litters.csv); samples are treated as independent.")
  }
}
test_txt <- if (!is_timecourse) {
  paste0("Stars: Welch's two-sample t-test (unpaired, unequal variances; R's t.test default) of each time point against ",
         ref_txt, ", separately for each marker; p-values are not adjusted for multiple testing",
         if (plot_by == "Cluster" && identical(s3$Settings$REFERENCE_SCOPE, "pooled"))
           "; the reference values are pooled over all plotted clusters" else "", ". ")
} else {
  posthoc_used <- paste(unique(st$PostHoc), collapse = "/"); omni_used <- paste(setdiff(unique(st$Omnibus_Test), NA), collapse = "/")
  switch(P_COLUMN,
    p_unadjusted = paste0("Stars: pairwise Welch t-test of each time point against ", ref_txt,
                          " (the cluster's own reference samples), separately for each marker, unadjusted for multiple testing. "),
    p_primary = paste0("Stars: post-hoc comparison of each time point with ", ref_txt, " (", posthoc_used,
                       "; welch_holm = pairwise Welch t-tests with Holm adjustment, dunnett = Dunnett's test), adjusted within each treatment and marker",
                       "; omnibus test across time points: ", omni_used, ". "),
    q_BH_all_markers = {
      fam_by <- intersect(c("Treatment", "Cluster", "ThreshCat"), names(st))
      fam <- st[!is.na(p_unadjusted), .(n_cmp = .N, n_mark = data.table::uniqueN(Marker), n_tp = data.table::uniqueN(TimePoint)), by = fam_by]
      rng <- function(x) if (min(x) == max(x)) as.character(min(x)) else paste0(min(x), "-", max(x))
      paste0("Stars: pairwise Welch t-tests of each time point against ", ref_txt,
             ", adjusted by Benjamini-Hochberg false discovery rate over all markers and time points of the treatment (q-values; ",
             rng(fam$n_mark), " markers x ", rng(fam$n_tp), " time points = ", rng(fam$n_cmp),
             " comparisons per family, so q depends on which markers were in the statistics run). ")
    })
}
legend_text <- paste0(
  "Points: individual samples; markers/line: mean; ",
  if (SD_STYLE == "ribbon") "shaded band: mean +/- SD. " else
    paste0("error bars: mean +/- SD (shown where a group has at least ", ERRORBAR_MIN_N, " samples). "), "Y axis: ",
  if (STAT_TABLE == "RAWvalues") paste0("the ", measure_txt, ", raw values (not normalised to a baseline); for percentages, of all cells in the sample. ")
  else if (STAT_TABLE == "FC") paste0("fold change of the ", measure_txt, " relative to the baseline (NoTreat 0 hr) mean. ")
  else paste0("log2 fold change of the ", measure_txt, " relative to the baseline (NoTreat 0 hr) mean. "),
  test_txt, "* ", if (P_COLUMN == "q_BH_all_markers") "q" else "p", " < 0.05, ** < 0.01, *** < 0.001. ",
  n_txt, litter_txt,
  if (any(plotted_data$Suppressed_Low_N %in% TRUE)) {
    sp <- unique(plotted_data[Suppressed_Low_N %in% TRUE, .(Group, TimePoint)])
    sp <- sp[, .(tp = paste(sort(unique(TimePoint)), collapse = ", ")), by = Group]
    paste0(" Time points with fewer than ", MIN_N_FOR_STARS, " samples in the group (",
           paste0(sp$Group, ": ", sp$tp, " hr", collapse = "; "), ") are shown without statistics.")
  } else "",
  if (shapes_vary) paste0(" Marker shapes and line types identify the ", tolower(plot_by), "s (", paste0(names(group_shapes), ": ",
                          { nm <- c("19" = "circle", "16" = "circle", "17" = "triangle", "15" = "square", "18" = "diamond", "4" = "cross")[as.character(group_shapes)]
                            ifelse(is.na(nm), paste0("shape ", group_shapes), nm) }, collapse = ", "), ").") else "",
  if (any(plotted_data$Is_Appended_0hr)) paste0(" Lines start at the shared ", REFERENCE_TREATMENT, " ", REFERENCE_TIMEPOINT,
                                                " hr baseline (the other treatments have no separate ", REFERENCE_TIMEPOINT, " hr sample).") else "",
  if (length(unstarred) > 0) paste0(" Shown without statistics: ", paste(unstarred, collapse = ", "), ".") else "")
caption_short <- paste0(sub(" Samples per time point.*$", "", legend_text))



#### Save ######################################################################
stage_label <- paste0("Linegraphs_By", plot_by, if (!is.null(PLOT_CATEGORY)) paste0("_Cat", PLOT_CATEGORY) else "", "_", P_COLUMN, if (nzchar(OUTPUT_LABEL)) paste0("_", OUTPUT_LABEL) else "")
Stage_WD <- make_stage_dir(stat_dir, stage_label)

n_slots <- PANEL_COLUMNS * ceiling(length(panels) / PANEL_COLUMNS)           # A2: ceiling, not round
panel_list <- unname(panels)
if (length(panel_list) < n_slots) {
  panel_list <- c(panel_list, replicate(n_slots - length(panel_list), ggplot2::ggplot() + ggplot2::theme_void(), simplify = FALSE))
}
n_rows <- n_slots / PANEL_COLUMNS
grid_panels <- gridExtra::arrangeGrob(grobs = panel_list, nrow = n_rows, ncol = PANEL_COLUMNS)
title_text <- if (is_pseudobulk) "PseudoBulk Plots" else paste0("Cluster", paste(clusters_plot, collapse = ","), " Plots")
title_grob <- grid::textGrob(title_text, gp = grid::gpar(fontsize = 18, fontface = "bold", col = bg$ink))
if (CAPTION_ON_FIGURE) {
  cap_grob <- grid::textGrob(paste(strwrap(caption_short, width = 38 * PANEL_COLUMNS), collapse = "\n"),
                             x = 0.01, hjust = 0, gp = grid::gpar(fontsize = 10, col = bg$ink))
  full_plot <- gridExtra::arrangeGrob(title_grob, grid_panels, cap_grob, ncol = 1, heights = c(0.05, 0.89, 0.06))
} else {
  full_plot <- gridExtra::arrangeGrob(title_grob, grid_panels, ncol = 1, heights = c(0.05, 0.95))
}
fig_width  <- PANEL_SIZE_IN * PANEL_COLUMNS
fig_height <- PANEL_SIZE_IN * n_rows

base_name <- paste0("Linegraph_", if (is_pseudobulk) "PseudoBulk_" else "Clusters_", "By", plot_by, "-Grid",
                    if (!is.null(PLOT_CATEGORY)) paste0("-Cat", PLOT_CATEGORY) else "")
out_files <- character(0)
for (fmt in OUTPUT_FORMATS) {
  f <- file.path(Stage_WD, paste0(base_name, ".", fmt))
  ok <- tryCatch({
    if (fmt == "pdf") {
      # R's pdf() device cuts file paths longer than 511 characters (the file is then saved WITHOUT its ending), and
      # the nested Analysis folders can exceed that. So the PDF is written to a short temporary path and copied.
      tmp_pdf <- tempfile(fileext = ".pdf")
      ggplot2::ggsave(tmp_pdf, plot = full_plot, width = fig_width, height = fig_height, units = "in",
                      device = "pdf", useDingbats = FALSE, limitsize = FALSE, bg = palette_save_bg(bg))
      if (!file.copy(tmp_pdf, f, overwrite = TRUE)) stop("could not copy the PDF to ", f)
      unlink(tmp_pdf)
    } else if (has_ragg) {
      ggplot2::ggsave(f, plot = full_plot, width = fig_width, height = fig_height, units = "in",
                      dpi = PNG_DPI, device = ragg::agg_png, limitsize = FALSE, bg = palette_save_bg(bg))
    } else {
      ggplot2::ggsave(f, plot = full_plot, width = fig_width, height = fig_height, units = "in", dpi = PNG_DPI, limitsize = FALSE, bg = palette_save_bg(bg))
    }
    TRUE
  }, error = function(e) { note_warning(paste0("Could not write ", basename(f), ": ", conditionMessage(e))); FALSE })
  if (ok) out_files <- c(out_files, f)
}
f_pd  <- file.path(Stage_WD, "Plotted_data.csv");     data.table::fwrite(plotted_data,   f_pd)
f_pp  <- file.path(Stage_WD, "Plotted_points.csv");   data.table::fwrite(plotted_points, f_pp)
f_ax  <- file.path(Stage_WD, "Plotted_axes.csv");     data.table::fwrite(data.table::rbindlist(axes_list), f_ax)
f_sty <- file.path(Stage_WD, "Plotted_style.csv");    data.table::fwrite(data.table::data.table(
  Group = group_levels, Color = unname(group_colors), Linetype = unname(group_linetypes), Shape = unname(group_shapes), PlotBy = plot_by, Table = STAT_TABLE,
  Background = BACKGROUND, Transparent_On = TRANSPARENT_ON), f_sty)
f_leg <- file.path(Stage_WD, "Figure_legend.txt")
writeLines(c(strwrap(legend_text, width = 100),
             "", paste0("Statistics folder: ", basename(stat_dir)), paste0("P-value column: ", P_COLUMN),
             paste0("Table: ", basename(table_path))), f_leg)
out_files <- c(out_files, f_pd, f_pp, f_ax, f_sty, f_leg)
if (any(!file.exists(out_files))) {
  note_warning(paste0("Not found after saving (check the folder): ", paste(basename(out_files[!file.exists(out_files)]), collapse = ", ")))
}
cat("[INFO] Written to:\n       ", Stage_WD, "\n", sep = "")



#### Run record ################################################################
Inputs <- list(
  Script_Version = SCRIPT_VERSION,
  Completed_At   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER   = stat_dir,
  OUTPUT_FOLDER  = Stage_WD,
  Settings = list(
    P_COLUMN = P_COLUMN, PLOT_BY = plot_by, SD_STYLE = SD_STYLE, ERRORBAR_MIN_N = ERRORBAR_MIN_N, MIN_N_FOR_STARS = MIN_N_FOR_STARS, Y_LIMITS = as.numeric(Y_LIMITS), STAR_STEP = if (is.na(STAR_STEP)) "auto" else STAR_STEP,
    NO_STARS_TREATMENTS = as.list(NO_STARS_TREATMENTS), MARKERS = as.list(markers_plot), TREATMENTS = as.list(treatments_plot),
    CLUSTERS = as.list(clusters_plot), CATEGORY = if (is.null(PLOT_CATEGORY)) "none" else PLOT_CATEGORY, OUTPUT_FORMATS = as.list(OUTPUT_FORMATS), PNG_DPI = PNG_DPI,
    PANEL_COLUMNS = PANEL_COLUMNS, PANEL_SIZE_IN = PANEL_SIZE_IN, CAPTION_ON_FIGURE = CAPTION_ON_FIGURE,
    BACKGROUND = BACKGROUND, TRANSPARENT_ON = TRANSPARENT_ON,
    Palette_File = basename(palette_file()), Palette_Colors = as.list(group_colors), Palette_Linetypes = as.list(group_linetypes), POINT_SHAPES = if (is.character(POINT_SHAPES)) POINT_SHAPES else as.list(POINT_SHAPES), Point_Shapes = as.list(group_shapes),
    Palette_Adjusted_For_Contrast = as.list(palette_res$lifted)),
  Statistics_Folder      = basename(stat_dir),
  Statistics_Source      = if (is_timecourse) "04_TimecourseStatistics" else "03_Statistics",
  Statistics_File_MD5    = unname(tools::md5sum(stat_file)),
  Input_Table            = basename(table_path),
  Input_Table_MD5        = unname(tools::md5sum(table_path)),
  Panels                 = length(panels),
  Starred_Group_TimePoints = n_stars,
  N_Mismatch_With_Statistics = n_mismatch,
  Legend_Text            = legend_text,
  PNG_Device             = if ("png" %in% OUTPUT_FORMATS) (if (has_ragg) "ragg::agg_png" else "grDevices::png") else "not used",
  Random_Seed            = "not used in this script",
  Warnings               = as.list(warnings_log),
  Output_MD5             = file_md5(out_files),
  Session                = session_info_list(c("data.table", "jsonlite", "ggplot2", "gridExtra", "viridis")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = stat_dir)

cat("\n========================================================\n")
cat("[DONE] ", length(panels), " panels, ", n_stars, " starred time points (p-value column: ", P_COLUMN, ")\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("       Legend text: Figure_legend.txt\n")
cat("========================================================\n")
