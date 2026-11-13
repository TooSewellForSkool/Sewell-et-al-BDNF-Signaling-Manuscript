################################################################################
# 04_PlotCategoryHistograms.R -- histogram of each marker with the cut values drawn on it
################################################################################
# Layer:    Visualization (no statistics; reads the output of 02_AssignThresholdCategories.R)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-28
# Edited:   2026-10-28 (1.1: neighbouring level labels staggered)
# Version:  1.1   (NEW; replaces the histogram half of 13.1_ExtractByMarkerExpression_HiLo_v5.R)
#
# PURPOSE
#   Shows where the cut values that define the expression categories (None / Lo / Med / Hi)
#   sit on each marker's distribution, with the percent of all cells in every level. Made for the
#   supplement of the revision (Reviewer question: how were the receptor tiers chosen?). Run it once
#   for the published cuts and, for the sensitivity check, once per shifted run.
#
# HOW TO USE
#   1. Run 02_AssignThresholdCategories.R (v2 or later). 2. Source this script, pick that folder.
#
# WHAT YOU GET (in a new folder next to the category folder)
#   CategoryHistograms.pdf / .png   one panel per marker, dashed lines = cut values
#   Plotted_data.csv                the numbers drawn (bins, cut values, level percents)
#   Figure_legend.txt               legend sentence with the cut values, scale and cell count
#   run_parameters.json             settings and checksums
#
# Look: grey bars (density), dashed black cut lines, level label and percent above each level.
#   The published 13.1 histograms also drew a magenta density curve; this version draws bars only.
#   Bin width: 100 equal bins over the marker's range (13.1 used 175 bins); the shape is the same.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
CATEGORY_FOLDER <- ""     # a folder made by 02_AssignThresholdCategories.R ("" = choose in a pop-up)

# Written on the axis and in the legend. Check it against panel.csv (columns Transform, Factor).
SCALE_NOTE      <- "arcsinh(x / 5)"

# Names for the levels, per axis. Leave c() to print the level numbers only.
# Example for the receptor tiers:  LEVEL_NAMES <- list(x = c("None", "Lo", "Med", "Hi"), y = c("None", "Lo", "Med", "Hi"))
LEVEL_NAMES     <- list()

OUTPUT_FORMATS  <- c("pdf", "png")
PNG_DPI         <- 300
PANEL_WIDTH_IN  <- 6
PANEL_HEIGHT_IN <- 4
OUTPUT_LABEL    <- ""     # optional text added to the output folder name

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "04_PlotCategoryHistograms.R"
SCRIPT_VERSION           <- "1.1"
MODULE_BASE_NAME         <- "04_PlotCategoryHistograms"
REQUIRED_HELPERS_VERSION <- 7
UPSTREAM_MODULE          <- "02_AssignThresholdCategories"
DEFS_NAME                <- "Category_Definitions.csv"
HIST_NAME                <- "Histogram_Data.csv"


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


#### Set up environment ########################################################
pkg_versions <- dc_load_packages(cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite", "ggplot2", "gridExtra"))
has_ragg <- requireNamespace("ragg", quietly = TRUE)

cat("\n========================================================\n")
cat("[START] Histograms with cut values:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

warnings_log <- character(0)
note_warning <- function(txt) {
  warnings_log <<- c(warnings_log, txt)
  cat("[WARNING] ", txt, "\n", sep = "")
}


#### Choose the category folder ################################################
if (nzchar(CATEGORY_FOLDER)) {
  if (!dir.exists(CATEGORY_FOLDER)) stop("CATEGORY_FOLDER not found:\n  ", CATEGORY_FOLDER, call. = FALSE)
  cat_dir <- normalizePath(CATEGORY_FOLDER, winslash = "/", mustWork = TRUE)
} else {
  Original_WD <- get_working_dir("project")
  dirs <- list.dirs(Original_WD, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, DEFS_NAME)) & file.exists(file.path(dirs, HIST_NAME))]
  dirs <- sort(dirs, decreasing = TRUE)
  if (length(dirs) == 0) {
    stop("No folder with ", DEFS_NAME, " and ", HIST_NAME, " was found under:\n  ", Original_WD,
         "\nRun 02_AssignThresholdCategories.R (version 2 or later) first, or set CATEGORY_FOLDER.", call. = FALSE)
  }
  pick <- ask_choice("Which category folder should be plotted?", basename(dirs), preselect = basename(dirs)[1])
  cat_dir <- normalizePath(dirs[basename(dirs) == pick][1], winslash = "/", mustWork = TRUE)
}
cat("[INFO] Category folder:\n       ", cat_dir, "\n", sep = "")

up_record <- read_run_record(cat_dir)
if (is.null(up_record) || is.null(up_record[[UPSTREAM_MODULE]])) {
  note_warning(paste0("No ", UPSTREAM_MODULE, " section in the run record of this folder; provenance cannot be checked."))
}
defs <- data.table::fread(file.path(cat_dir, DEFS_NAME), data.table = FALSE)
hist_data <- data.table::fread(file.path(cat_dir, HIST_NAME), data.table = FALSE)
axes <- unique(defs$Axis)
n_cells <- if (!is.null(up_record[[UPSTREAM_MODULE]]$Cells)) as.integer(up_record[[UPSTREAM_MODULE]]$Cells) else sum(hist_data$Cells[hist_data$Axis == axes[1]])
shift_txt <- as.character(up_record[[UPSTREAM_MODULE]]$Settings$CUT_SHIFT)
cuts_by <- unique(defs$Cuts_Chosen_By)


#### One panel per marker ######################################################
make_panel <- function(ax) {
  h <- hist_data[hist_data$Axis == ax, , drop = FALSE]
  d <- defs[defs$Axis == ax, , drop = FALSE]
  mk <- h$Marker[1]
  width <- h$Bin_Upper - h$Bin_Lower
  h$Density <- h$Cells / (sum(h$Cells) * width)
  cuts <- stats::na.omit(unique(c(d$Lower_Cut, d$Upper_Cut)))
  rng <- c(min(h$Bin_Lower), max(h$Bin_Upper))
  # label position: middle of each level (clipped to the data range)
  lo <- ifelse(is.na(d$Lower_Cut), rng[1], d$Lower_Cut)
  hi <- ifelse(is.na(d$Upper_Cut), rng[2], d$Upper_Cut)
  nm <- LEVEL_NAMES[[ax]]
  lab_name <- if (length(nm) == nrow(d)) paste0(d$Level, " (", nm, ")") else as.character(d$Level)
  labs_df <- data.frame(x = (lo + hi) / 2, Label = paste0(lab_name, "\n", round(d$Percent_Of_Cells, 1), "%"))
  ymax <- max(h$Density)
  # stagger neighbouring labels so narrow levels (e.g. TrkB_i levels 1 and 2) do not overprint
  labs_df$y <- ymax * 1.15 - ymax * 0.17 * ((seq_len(nrow(labs_df)) %% 2) == 0)
  p <- ggplot2::ggplot(h, ggplot2::aes(xmin = Bin_Lower, xmax = Bin_Upper, ymin = 0, ymax = Density)) +
    ggplot2::geom_rect(fill = "grey70", colour = NA) +
    ggplot2::geom_vline(xintercept = cuts, linetype = "longdash") +
    ggplot2::geom_text(data = labs_df, ggplot2::aes(x = x, y = y, label = Label), inherit.aes = FALSE, size = 3.5, vjust = 1) +
    ggplot2::scale_y_continuous(limits = c(0, ymax * 1.2), expand = c(0, 0)) +
    ggplot2::labs(x = paste0(mk, "  (", SCALE_NOTE, ")"), y = "Density", title = mk) +
    ggplot2::theme_bw()
  list(plot = p, data = data.frame(Marker = mk, Axis = ax, h[, c("Bin_Lower", "Bin_Upper", "Cells", "Density")], stringsAsFactors = FALSE),
       cuts = data.frame(Marker = mk, Axis = ax, Cut = cuts), levels = d)
}
panels <- lapply(axes, make_panel)
plot_list <- lapply(panels, function(x) x$plot)


#### Legend text ###############################################################
cut_txt <- paste(vapply(panels, function(x) {
  paste0(x$cuts$Marker[1], " at ", paste(signif(x$cuts$Cut, 4), collapse = ", "))
}, character(1)), collapse = "; ")
legend_text <- paste0(
  "Distribution of ", paste(vapply(panels, function(x) x$cuts$Marker[1], character(1)), collapse = " and "),
  " in all ", format(n_cells, big.mark = ","), " analysed cells (all treatments and time points; ", SCALE_NOTE,
  " values after batch correction). Dashed lines mark the cut values that define the expression levels (", cut_txt,
  "; ", paste(cuts_by, collapse = "; "), "). A cell exactly on a cut value is assigned to the higher level. ",
  "Labels give each level and the percent of all cells in it.")


#### Save ######################################################################
Stage_WD <- make_stage_dir(dirname(cat_dir), paste0("CategoryHistograms", if (nzchar(OUTPUT_LABEL)) paste0("_", OUTPUT_LABEL) else ""))
out_files <- character(0)
ncol_f <- length(plot_list)
full_plot <- gridExtra::arrangeGrob(grobs = plot_list, nrow = 1, ncol = ncol_f)
for (fmt in OUTPUT_FORMATS) {
  f <- file.path(Stage_WD, paste0("CategoryHistograms.", fmt))
  ok <- tryCatch({
    if (fmt == "pdf") {
      tmp_pdf <- tempfile(fileext = ".pdf")   # R's pdf() cuts paths longer than 511 characters: write short, then copy
      ggplot2::ggsave(tmp_pdf, plot = full_plot, width = PANEL_WIDTH_IN * ncol_f, height = PANEL_HEIGHT_IN, units = "in",
                      device = "pdf", useDingbats = FALSE, limitsize = FALSE)
      if (!file.copy(tmp_pdf, f, overwrite = TRUE)) stop("could not copy the PDF to ", f)
      unlink(tmp_pdf)
    } else if (has_ragg) {
      ggplot2::ggsave(f, plot = full_plot, width = PANEL_WIDTH_IN * ncol_f, height = PANEL_HEIGHT_IN, units = "in",
                      dpi = PNG_DPI, device = ragg::agg_png, limitsize = FALSE, bg = "white")
    } else {
      ggplot2::ggsave(f, plot = full_plot, width = PANEL_WIDTH_IN * ncol_f, height = PANEL_HEIGHT_IN, units = "in",
                      dpi = PNG_DPI, limitsize = FALSE, bg = "white")
    }
    TRUE
  }, error = function(e) { note_warning(paste0("Could not write ", basename(f), ": ", conditionMessage(e))); FALSE })
  if (ok) out_files <- c(out_files, f)
}
f_pd <- file.path(Stage_WD, "Plotted_data.csv")
plotted <- do.call(rbind, lapply(panels, function(x) x$data))
data.table::fwrite(plotted, f_pd)
f_lv <- file.path(Stage_WD, "Plotted_levels.csv")
data.table::fwrite(do.call(rbind, lapply(panels, function(x) x$levels)), f_lv)
f_leg <- file.path(Stage_WD, "Figure_legend.txt")
writeLines(c(strwrap(legend_text, width = 100), "", paste0("Category folder: ", basename(cat_dir))), f_leg)
out_files <- c(out_files, f_pd, f_lv, f_leg)
cat("[INFO] Written to:\n       ", Stage_WD, "\n", sep = "")


#### Run record ################################################################
Inputs <- list(
  Script_Version = SCRIPT_VERSION,
  Completed_At   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER   = cat_dir,
  OUTPUT_FOLDER  = Stage_WD,
  Settings = list(SCALE_NOTE = SCALE_NOTE, LEVEL_NAMES = LEVEL_NAMES, OUTPUT_FORMATS = as.list(OUTPUT_FORMATS),
                  PNG_DPI = PNG_DPI, PANEL_WIDTH_IN = PANEL_WIDTH_IN, PANEL_HEIGHT_IN = PANEL_HEIGHT_IN),
  Category_Folder = basename(cat_dir),
  Cut_Shift_In_Category_Run = shift_txt,
  Cells = n_cells,
  Legend_Text = legend_text,
  Random_Seed = "not used in this script",
  Warnings = as.list(warnings_log),
  Output_MD5 = file_md5(out_files),
  Session = session_info_list(c("data.table", "jsonlite", "ggplot2", "gridExtra")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = cat_dir)

cat("\n========================================================\n")
cat("[DONE] ", length(plot_list), " histogram panel(s)\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("       Legend text: Figure_legend.txt\n")
cat("========================================================\n")
