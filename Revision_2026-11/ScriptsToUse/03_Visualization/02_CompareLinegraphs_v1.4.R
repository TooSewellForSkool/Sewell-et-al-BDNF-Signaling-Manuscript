################################################################################
# 02_CompareLinegraphs.R
################################################################################
# Layer:    Visualization -- OPTIONAL comparison of two versions of the same figure (not part of the normal workflow)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-26
# Version:  1   (NEW; no published counterpart)
#
# Purpose:
# * Put two runs of 01_PlotLinegraphs.R side by side, marker by marker: A (e.g. the published statistics) on the
#   left, B (e.g. the revised statistics) on the right. The lines, dots and axes must be the same in both
#   (the script checks this); only the stars differ.
# * Write Comparison_changes.csv: every plotted group/time point whose stars differ between A and B
#   (gained, lost, stronger or weaker), plus a summary.
# * Use it for the revision (published vs revised statistics) or to compare any two p-value choices. It is not run
#   by 01_PlotLinegraphs.R and nothing depends on it.
#
# What you need: two output folders of 01_PlotLinegraphs.R made from the SAME table, selection and markers
# (each holds Plotted_data.csv, Plotted_points.csv, Plotted_axes.csv, Plotted_style.csv).
#
# CHANGE LOG: v1 new script. v1.1 figure saved with limitsize = FALSE (tall figures >50 in were refused, so no PDF/PNG was written); prints an ERROR if no figure is written. v1.2 star layer no longer drawn in the legend (removes the large "a" key glyph). v1.3 line types and background follow the 01_PlotLinegraphs.R run (Plotted_style.csv); a 01 run made before v1.10 has none and is drawn as before (solid lines, published grey panel). v1.4 draws the mean points with the shapes recorded in Plotted_style.csv by 01_PlotLinegraphs.R v1.12 (older runs: round dots); nothing else changes.
#  A1 [BEHAVIOR] Stops if the two runs do not describe the same lines (same markers, groups, time points, n,
#     mean and SD within 1e-9): the comparison is only meaningful when only the p-values differ.
#  A2 [BEHAVIOR] The star positions follow 01_PlotLinegraphs.R (same rule). Look and layout follow 01 as well.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
FOLDER_A   <- ""               # an output folder of 01_PlotLinegraphs.R: the REFERENCE, drawn on the left ("" = pop-up)
FOLDER_B   <- ""               # an output folder of 01_PlotLinegraphs.R: the NEW version, drawn on the right ("" = pop-up)
LABEL_A    <- ""               # text above the left panels  ("" = taken from the run: the p-value column used)
LABEL_B    <- ""               # text above the right panels ("" = taken from the run)

MARKERS             <- c()     # c() = all markers the two runs share
ONLY_CHANGED_MARKERS <- FALSE  # TRUE: draw only the markers where at least one star differs
PAIRS_PER_ROW       <- 2       # marker pairs (A | B) per row of the figure

OUTPUT_FORMATS <- c("pdf", "png")
PNG_DPI        <- 300
PANEL_SIZE_IN  <- 5
OUTPUT_LABEL   <- ""

#### Fixed names -- do not change ###############################################
SCRIPT_NAME              <- "02_CompareLinegraphs.R"
SCRIPT_VERSION           <- "1.4"
MODULE_BASE_NAME         <- "02_CompareLinegraphs"
REQUIRED_HELPERS_VERSION <- 7



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
source(file.path(HELPERS_DIR, "01_Palette.R"), local = FALSE)   # v1.3: background, theme and line types as in 01_PlotLinegraphs.R



#### Set up environment ########################################################
pkg_versions <- dc_load_packages(cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite", "ggplot2", "gridExtra"))
has_ragg <- requireNamespace("ragg", quietly = TRUE)
cat("\n========================================================\n")
cat("[START] Side-by-side comparison of two line-graph runs:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")
warnings_log <- character(0)
note_warning <- function(txt) { warnings_log <<- c(warnings_log, txt); cat("[WARNING] ", txt, "\n", sep = "") }
one_of <- function(v, ok, nm) if (!all(v %in% ok) || length(v) == 0) stop(nm, " must be: ", paste(ok, collapse = " / "), call. = FALSE)
one_of(OUTPUT_FORMATS, c("pdf", "png"), "OUTPUT_FORMATS")
for (nm in c("PNG_DPI", "PANEL_SIZE_IN", "PAIRS_PER_ROW")) {
  v <- get(nm); if (!is.numeric(v) || length(v) != 1 || v <= 0) stop(nm, " must be one positive number.", call. = FALSE)
}
if (!is.logical(ONLY_CHANGED_MARKERS) || length(ONLY_CHANGED_MARKERS) != 1 || is.na(ONLY_CHANGED_MARKERS)) {
  stop("ONLY_CHANGED_MARKERS must be TRUE or FALSE.", call. = FALSE)
}



#### Choose the two folders ####################################################
needed <- c("Plotted_data.csv", "Plotted_points.csv", "Plotted_axes.csv", "Plotted_style.csv")
is_run_folder <- function(d) all(file.exists(file.path(d, needed)))
choose_folder <- function(path, what, exclude = NULL) {
  if (nzchar(path)) {
    if (!dir.exists(path)) stop(what, " not found:\n  ", path, call. = FALSE)
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  Original_WD <- get_working_dir("project")
  dirs <- list.dirs(Original_WD, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[!grepl("/[_.][^/]*", sub(Original_WD, "", dirs, fixed = TRUE))]
  dirs <- sort(dirs[vapply(dirs, is_run_folder, logical(1))], decreasing = TRUE)
  dirs <- setdiff(dirs, exclude)
  if (length(dirs) == 0) {
    stop("No line-graph output folder (with ", paste(needed, collapse = ", "), ") found in:\n  ", Original_WD,
         "\nRun 01_PlotLinegraphs.R (the current version) first, or set ", what, ".", call. = FALSE)
  }
  pick <- ask_choice(paste0("Which folder is ", what, "?"), basename(dirs), preselect = basename(dirs)[1])
  normalizePath(dirs[match(pick, basename(dirs))], winslash = "/", mustWork = TRUE)
}
dir_a <- choose_folder(FOLDER_A, "FOLDER_A (reference, left)")
dir_b <- choose_folder(FOLDER_B, "FOLDER_B (new version, right)", exclude = dir_a)
for (d in c(dir_a, dir_b)) if (!is_run_folder(d)) stop("Missing one of ", paste(needed, collapse = ", "), " in:\n  ", d, call. = FALSE)
if (identical(dir_a, dir_b)) stop("FOLDER_A and FOLDER_B are the same folder.", call. = FALSE)
cat("[INFO] A (left):  ", dir_a, "\n[INFO] B (right): ", dir_b, "\n", sep = "")

rd <- function(d, f) data.table::fread(file.path(d, f), na.strings = c("", "NA"), data.table = TRUE)
pa <- rd(dir_a, "Plotted_data.csv"); pb <- rd(dir_b, "Plotted_data.csv")
pts <- rd(dir_a, "Plotted_points.csv"); axes <- rd(dir_a, "Plotted_axes.csv"); sty <- rd(dir_a, "Plotted_style.csv")
rec_a <- read_run_record(dir_a); rec_b <- read_run_record(dir_b)
mod_a <- rec_a[["01_PlotLinegraphs"]]; mod_b <- rec_b[["01_PlotLinegraphs"]]
if (is.null(mod_a) || is.null(mod_b)) stop("Both folders need a 01_PlotLinegraphs section in run_parameters.json.", call. = FALSE)
for (d in list(pa, pb)) { d[, Group := as.character(Group)]; d[, Stars := ifelse(is.na(Stars), "", Stars)] }
label_a <- if (nzchar(LABEL_A)) LABEL_A else paste0("A: ", mod_a$Statistics_Source, " / ", mod_a$Settings$P_COLUMN)
label_b <- if (nzchar(LABEL_B)) LABEL_B else paste0("B: ", mod_b$Statistics_Source, " / ", mod_b$Settings$P_COLUMN)



#### The two runs must describe the same lines (A1) ############################
key <- c("Marker", "Group", "TimePoint", "Is_Appended_0hr")
setkeyv(pa, key); setkeyv(pb, key)
if (nrow(pa) != nrow(pb) || !identical(pa[, ..key], pb[, ..key])) {
  stop("The two runs do not plot the same markers / groups / time points.\n",
       "A: ", nrow(pa), " rows; B: ", nrow(pb), " rows. Make both runs from the same table with the same\n",
       "MARKERS, TREATMENTS and CLUSTERS (only P_COLUMN / the Statistics folder may differ).", call. = FALSE)
}
for (cn in c("N", "Mean", "SD")) {
  d <- max(abs(pa[[cn]] - pb[[cn]]), na.rm = TRUE)
  if (is.finite(d) && d > 1e-9) stop("Column ", cn, " differs between the runs (max difference ", signif(d, 3),
                                    "): they were not made from the same table/selection.", call. = FALSE)
}
if (!identical(unique(as.character(mod_a$Input_Table_MD5)), unique(as.character(mod_b$Input_Table_MD5)))) {
  note_warning("The two runs used tables with different checksums.")
}
cat("[INFO] Same lines in A and B: ", nrow(pa), " rows, ", data.table::uniqueN(pa$Marker), " markers\n", sep = "")



#### What changed ###############################################################
rank_star <- function(x) nchar(x)
chg <- merge(pa[Is_Appended_0hr == FALSE, .(Marker, Group, TimePoint, N, Mean, p_A = p, Stars_A = Stars, P_Column_A = P_Column)],
             pb[Is_Appended_0hr == FALSE, .(Marker, Group, TimePoint, p_B = p, Stars_B = Stars, P_Column_B = P_Column)],
             by = c("Marker", "Group", "TimePoint"), sort = FALSE)
chg[, Change := ifelse(Stars_A == Stars_B, "same",
                ifelse(Stars_A == "", "gained",
                ifelse(Stars_B == "", "lost",
                ifelse(rank_star(Stars_B) > rank_star(Stars_A), "stronger", "weaker"))))]
ref_rows <- chg$TimePoint == 0 & chg$Stars_A == "" & chg$Stars_B == ""
n_cmp <- sum(!is.na(chg$p_A) | !is.na(chg$p_B))
summ_tab <- chg[, .N, by = Change][order(-N)]
cat("[INFO] Group/time points compared: ", nrow(chg), "\n", sep = "")
for (i in seq_len(nrow(summ_tab))) cat("       ", format(summ_tab$Change[i], width = 9), summ_tab$N[i], "\n")
chg_out <- chg[order(Change == "same", Marker, Group, TimePoint)]
markers_changed <- unique(chg[Change != "same", Marker])
markers_all <- unique(pa$Marker)
markers_draw <- if (length(MARKERS) > 0) intersect(MARKERS, markers_all) else markers_all
if (length(MARKERS) > 0 && length(setdiff(MARKERS, markers_all)) > 0) {
  note_warning(paste0("MARKERS not in the runs: ", paste(setdiff(MARKERS, markers_all), collapse = ", ")))
}
if (ONLY_CHANGED_MARKERS) markers_draw <- intersect(markers_draw, markers_changed)
if (length(markers_draw) == 0) stop("No markers to draw (check MARKERS / ONLY_CHANGED_MARKERS).", call. = FALSE)
cat("[INFO] Markers drawn: ", length(markers_draw), " (", length(markers_changed), " with at least one changed star)\n", sep = "")



#### Panels (same rules as 01_PlotLinegraphs.R) ################################
group_levels <- sty$Group
group_colors <- stats::setNames(sty$Color, sty$Group)
# v1.3: line types and background are the ones recorded by 01_PlotLinegraphs.R (older runs without them: solid, published grey panel)
group_linetypes <- stats::setNames(if ("Linetype" %in% names(sty)) sty$Linetype else rep("solid", nrow(sty)), sty$Group)
# v1.4: point shapes recorded by 01_PlotLinegraphs.R v1.12 (older runs without them: round dots)
group_shapes <- stats::setNames(if ("Shape" %in% names(sty)) as.integer(sty$Shape) else rep(19L, nrow(sty)), sty$Group)
shapes_vary <- length(unique(group_shapes)) > 1
bg <- palette_background(if ("Background" %in% names(sty)) sty$Background[1] else "published",
                         if ("Transparent_On" %in% names(sty)) sty$Transparent_On[1] else "white")
plot_by <- sty$PlotBy[1]
tp_levels <- as.character(sort(unique(pa$TimePoint)))
star_label_rows <- function(d) d[Stars != "" & Is_Appended_0hr == FALSE]
make_panel <- function(m, stars_tbl, panel_label) {
  d <- pa[Marker == m]
  d_df <- as.data.frame(d)
  d_df$Group <- factor(as.character(d_df$Group), levels = group_levels)
  d_df$TimePoint_f <- factor(as.character(d_df$TimePoint), levels = tp_levels)
  pt <- as.data.frame(pts[Marker == m])
  pt$Group <- factor(as.character(pt$Group), levels = group_levels)
  pt$TimePoint_f <- factor(as.character(pt$TimePoint), levels = tp_levels)
  without0 <- unique(d$Treatment[d$Is_Appended_0hr]); with0 <- setdiff(unique(d$Treatment), without0)
  d_without0 <- d_df[d_df$Treatment %in% without0, ]; d_with0 <- d_df[d_df$Treatment %in% with0, ]
  ax <- axes[Marker == m]
  st <- as.data.frame(star_label_rows(stars_tbl[Marker == m]))
  p <- ggplot2::ggplot(pt, ggplot2::aes(x = TimePoint_f, y = Value, color = Group)) +
    ggplot2::geom_point(size = 1, alpha = 0.5, show.legend = !shapes_vary) +
    ggplot2::geom_ribbon(data = d_df, ggplot2::aes(x = TimePoint_f, y = Mean, ymin = Mean - SD, ymax = Mean + SD,
                                                  fill = Group, group = Group),
                         alpha = 0.2, color = NA, inherit.aes = FALSE) +
    ggplot2::geom_point(data = if (shapes_vary) d_without0[!d_without0$Is_Appended_0hr, ] else d_without0,
                        ggplot2::aes(y = Mean, shape = Group), size = 2) +
    ggplot2::geom_line(data = d_without0, ggplot2::aes(y = Mean, group = Group, linetype = Group), linewidth = 0.5) +
    ggplot2::geom_point(data = d_with0, ggplot2::aes(y = Mean, shape = Group), size = 2) +
    ggplot2::labs(title = m, subtitle = panel_label, color = plot_by, fill = plot_by, linetype = plot_by, shape = plot_by, x = "Time (hr)", y = ax$Y_Label[1]) +
    ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = group_linetypes, drop = FALSE) +
    ggplot2::scale_shape_manual(values = group_shapes, drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(ax$y_lower[1], ax$y_upper[1] * 1.12)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)) +
    palette_theme(bg)
  if (nrow(st) > 0) {
    st$Group <- factor(as.character(st$Group), levels = group_levels)
    st$TimePoint_f <- factor(as.character(st$TimePoint), levels = tp_levels)
    st$y <- (ax$y_upper[1] * 1.06) - (3 * (as.integer(st$Group) - 1))
    p <- p + ggplot2::geom_text(data = st, ggplot2::aes(x = TimePoint_f, y = y, label = Stars, color = Group),
                                inherit.aes = FALSE, size = 10, show.legend = FALSE)
  }
  p
}
panel_list <- list()
for (m in markers_draw) {
  panel_list[[length(panel_list) + 1]] <- make_panel(m, pa, label_a)
  panel_list[[length(panel_list) + 1]] <- make_panel(m, pb, label_b)
}
n_cols <- 2 * PAIRS_PER_ROW
n_slots <- n_cols * ceiling(length(panel_list) / n_cols)
if (length(panel_list) < n_slots) {
  panel_list <- c(panel_list, replicate(n_slots - length(panel_list), ggplot2::ggplot() + ggplot2::theme_void(), simplify = FALSE))
}
n_rows <- n_slots / n_cols
grid_panels <- gridExtra::arrangeGrob(grobs = panel_list, nrow = n_rows, ncol = n_cols)
title_grob <- grid::textGrob(paste0("Comparison: ", label_a, "   |   ", label_b), gp = grid::gpar(fontsize = 18, fontface = "bold", col = bg$ink))
full_plot <- gridExtra::arrangeGrob(title_grob, grid_panels, ncol = 1, heights = c(0.05, 0.95))



#### Save ######################################################################
Stage_WD <- make_stage_dir(dir_b, paste0("Comparison_A-vs-B", if (nzchar(OUTPUT_LABEL)) paste0("_", OUTPUT_LABEL) else ""))
out_files <- character(0)
for (fmt in OUTPUT_FORMATS) {
  f <- file.path(Stage_WD, paste0("Linegraph_Comparison.", fmt))
  ok <- tryCatch({
    if (fmt == "pdf") {
      tmp_pdf <- tempfile(fileext = ".pdf")      # R's pdf() cuts paths longer than 511 characters: write short, then copy
      ggplot2::ggsave(tmp_pdf, plot = full_plot, width = PANEL_SIZE_IN * n_cols, height = PANEL_SIZE_IN * n_rows,
                      units = "in", device = "pdf", useDingbats = FALSE, limitsize = FALSE, bg = palette_save_bg(bg))
      if (!file.copy(tmp_pdf, f, overwrite = TRUE)) stop("could not copy the PDF to ", f)
      unlink(tmp_pdf)
    } else if (has_ragg) {
      ggplot2::ggsave(f, plot = full_plot, width = PANEL_SIZE_IN * n_cols, height = PANEL_SIZE_IN * n_rows,
                      units = "in", dpi = PNG_DPI, device = ragg::agg_png, limitsize = FALSE, bg = palette_save_bg(bg))
    } else {
      ggplot2::ggsave(f, plot = full_plot, width = PANEL_SIZE_IN * n_cols, height = PANEL_SIZE_IN * n_rows,
                      units = "in", dpi = PNG_DPI, limitsize = FALSE, bg = palette_save_bg(bg))
    }
    TRUE
  }, error = function(e) { note_warning(paste0("Could not write ", basename(f), ": ", conditionMessage(e))); FALSE })
  if (ok) out_files <- c(out_files, f)
}
f_chg <- file.path(Stage_WD, "Comparison_changes.csv"); data.table::fwrite(chg_out, f_chg)
f_sum <- file.path(Stage_WD, "Comparison_summary.txt")
writeLines(c(paste0("A (left):  ", basename(dir_a), "   [", label_a, "]"),
             paste0("B (right): ", basename(dir_b), "   [", label_b, "]"),
             paste0("Group/time points compared: ", nrow(chg)),
             paste0(summ_tab$Change, ": ", summ_tab$N),
             paste0("Markers with at least one changed star: ", paste(markers_changed, collapse = ", "))), f_sum)
out_files <- c(out_files, f_chg, f_sum)
if (!any(grepl("\\.(pdf|png)$", out_files))) {
  cat("\n[ERROR] No figure was written. Reason(s):\n"); for (w in warnings_log) cat("        ", w, "\n", sep = "")
}
cat("[INFO] Written to:\n       ", Stage_WD, "\n", sep = "")

Inputs <- list(
  Script_Version = SCRIPT_VERSION, Completed_At = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER = dir_b, OUTPUT_FOLDER = Stage_WD,
  Settings = list(FOLDER_A = dir_a, FOLDER_B = dir_b, LABEL_A = label_a, LABEL_B = label_b,
                  MARKERS = as.list(MARKERS), ONLY_CHANGED_MARKERS = ONLY_CHANGED_MARKERS,
                  PAIRS_PER_ROW = PAIRS_PER_ROW, OUTPUT_FORMATS = as.list(OUTPUT_FORMATS), PNG_DPI = PNG_DPI),
  Plotted_Data_A_MD5 = unname(tools::md5sum(file.path(dir_a, "Plotted_data.csv"))),
  Plotted_Data_B_MD5 = unname(tools::md5sum(file.path(dir_b, "Plotted_data.csv"))),
  Group_TimePoints_Compared = nrow(chg),
  Changes = as.list(stats::setNames(summ_tab$N, summ_tab$Change)),
  Random_Seed = "not used in this script",
  Warnings = as.list(warnings_log),
  Output_MD5 = file_md5(out_files),
  Session = session_info_list(c("data.table", "jsonlite", "ggplot2", "gridExtra")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = dir_b)

cat("\n========================================================\n")
cat("[DONE] ", length(markers_draw), " marker pairs; ", sum(chg$Change != "same"), " changed group/time points\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("========================================================\n")
