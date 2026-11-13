################################################################################
# 03_PlotBargraphs.R
################################################################################
# Layer:    Visualization (no statistics are calculated here; p-values are READ from the Analysis layer)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-26
# Edited:   2026-10-26
# Version:  2   (NEW; replaces 16_Barplot_v7.r)
#
# Purpose:
# * Draw the condition bar graphs of the paper (e.g. Cont. Depr. / BDNF test / BDNF + K252a test at 1 hr):
#   one panel per marker (pseudobulk) or one panel per cluster/category (one figure per marker), bars = mean,
#   error bars = SD, dots = the individual samples, brackets with stars for the three pairwise comparisons.
# * The stars come from the statistics folder you choose (made by 05_ConditionStatistics.R) and you choose
#   WHICH p-value draws them (P_COLUMN). ALL p-values of every drawn comparison are written to
#   Plotted_data.csv, and the figure legend text (test, adjustment, family size, n) to Figure_legend.txt.
#
# Look: the published look is kept on purpose (grey ggplot theme, viridis "plasma" colours, bars with alpha 0.5,
#   5-inch panels, 45-degree axis labels, bracket order bottom to top = 1 vs 2, 1 vs 3, 2 vs 3).
#   The published brackets printed "p = 0.0076" (unadjusted Welch); the brackets now print stars (or "ns").
#
# ==============================================================================
# CHANGE LOG: 16_Barplot_v7.r (published) -> v1
# ==============================================================================
#  A1 [RESULT, declared] Brackets show stars (and "ns") of the p-value column you choose (default p_primary:
#     adjusted over the 3 pairs) instead of "p = <unadjusted Welch p>". The published p-values are in
#     Plotted_data.csv (column p_unadjusted).
#  A2 [RESULT] No star is drawn if either group of a comparison has fewer than MIN_N_FOR_STARS samples (default 3;
#     decision 2026-10-26). The p-value stays in Plotted_data.csv.
#  A3 [BEHAVIOR] Means and SDs are read from Conditions_groupstats.csv (computed by 05 from the same table);
#     the dots from Conditions_input_values.csv. Nothing is recalculated here.
#  A4 [BEHAVIOR] The published script placed the first bracket at 1.1 x the highest DOT and the next ones one
#     step up (STEP_INCREASE x highest dot); the same rule is used here, but the axis top is raised if a bracket
#     would otherwise be cut off.
#  A5 [BEHAVIOR] The published y axis label was found by looking for a mean of 1; here it says
#     "FC relative to NoTreat 0 hr" (the table's reference).
#
# B. Efficiency / readability / reproducibility
#  B1 One run record (run_parameters.json), PDF (editable) and PNG in one go, no random numbers.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
STATS_FOLDER   <- ""            # a folder made by 05_ConditionStatistics.R ("" = choose in a pop-up)

# Which p-value draws the stars?
#   "p_primary"        (default) adjusted over the 3 pairs of the marker (default of 05 v2: Welch + Holm = p_Holm)
#   "p_Tukey"          Tukey HSD (classic ANOVA; assumes equal spread)
#   "p_Holm"           Welch t-tests + Holm over the 3 pairs
#   "p_unadjusted"     Welch t-test, no adjustment (what the original submission showed as p = ...)
#   "q_BH_all_markers" Benjamini-Hochberg over all markers x pairs of the cluster/category
P_COLUMN       <- "p_primary"

MARKERS        <- c("pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La")           
# e.g. c("pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La"); c() = every marker in the statistics folder
CLUSTERS       <- c()           # e.g. c(1, 2, 4, 3, 17, 18); c() = all (pseudobulk tables have none)
CLUSTER_LABELS <- c()           # optional names for the panels, e.g. c("1" = "Neuronal", "2" = "Non-neuronal")
PANEL_TITLE_PREFIX <- "1hr Treatment: "

MIN_N_FOR_STARS <- 3            # no star if a group of the comparison has fewer samples (0 = always draw)
SHOW_NS         <- TRUE         # TRUE: "ns" above brackets that are not significant; FALSE: bracket without text
STEP_INCREASE   <- 0.15         # separation of the brackets (fraction of the highest dot), as published
COLOR_OPTION    <- "C"          # viridis option: A magma, B inferno, C plasma (published), D viridis, E cividis

# Output
OUTPUT_FORMATS  <- c("pdf", "png")
PNG_DPI         <- 300
PANEL_COLUMNS   <- 1            # panels per row
PANEL_SIZE_IN   <- 5            # width and height of one panel (inches)
OUTPUT_LABEL    <- "Fig2B"           # optional text added to the output folder name

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "03_PlotBargraphs.R"
SCRIPT_VERSION           <- "2"
MODULE_BASE_NAME         <- "03_PlotBargraphs"
REQUIRED_HELPERS_VERSION <- 7
UPSTREAM_MODULE          <- "05_ConditionStatistics"



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
pkg_versions <- dc_load_packages(cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite",
                                          "ggplot2", "ggsignif", "gridExtra", "viridis"))
has_ragg <- requireNamespace("ragg", quietly = TRUE)

cat("\n========================================================\n")
cat("[START] Bar graphs with stars from the Statistics folder:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

warnings_log <- character(0)
note_warning <- function(txt) {
  warnings_log <<- c(warnings_log, txt)
  cat("[WARNING] ", txt, "\n", sep = "")
}
if (!P_COLUMN %in% c("p_primary", "p_Holm", "p_Tukey", "p_unadjusted", "q_BH_all_markers")) {
  stop("P_COLUMN must be one of \"p_primary\", \"p_Holm\", \"p_Tukey\", \"p_unadjusted\", \"q_BH_all_markers\".", call. = FALSE)
}
if (!is.numeric(MIN_N_FOR_STARS) || length(MIN_N_FOR_STARS) != 1 || MIN_N_FOR_STARS < 0) {
  stop("MIN_N_FOR_STARS must be one number, 0 or larger.", call. = FALSE)
}



#### Choose and read the statistics folder #####################################
is_stat_folder <- function(d) file.exists(file.path(d, "Conditions_long.csv")) &&
                              file.exists(file.path(d, "Conditions_groupstats.csv")) &&
                              file.exists(file.path(d, "Conditions_input_values.csv"))
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
  if (length(dirs) == 0) stop("No 05_ConditionStatistics.R output folder was found in:\n  ", Original_WD,
                              "\nRun 05_ConditionStatistics.R first, or set STATS_FOLDER.", call. = FALSE)
  labels <- basename(dirs)
  pick <- ask_choice("Which Statistics folder should be plotted?", labels, preselect = labels[1])
  stat_dir <- normalizePath(dirs[match(pick, labels)], winslash = "/", mustWork = TRUE)
}
cat("[INFO] Statistics folder:\n       ", stat_dir, "\n", sep = "")
record <- require_modules(stat_dir, UPSTREAM_MODULE, how_to_fix = "Run 05_ConditionStatistics.R first.")

rd <- function(f) {
  x <- data.table::fread(file.path(stat_dir, f), na.strings = c("", "NA"), data.table = TRUE)
  if ("Cluster" %in% names(x)) x[, Cluster := as.character(Cluster)]
  if ("ThreshCat" %in% names(x)) { x[, ThreshCat := as.character(ThreshCat)]; x[is.na(ThreshCat), ThreshCat := ""] }
  x
}
long  <- rd("Conditions_long.csv")
fams  <- rd("Conditions_families.csv")
gstat <- rd("Conditions_groupstats.csv")
vals  <- rd("Conditions_input_values.csv")
if (!"p_Tukey" %in% names(long)) long[, p_Tukey := NA_real_]   # folders made by 05 version 1 have no p_Tukey
if (!P_COLUMN %in% names(long)) stop("Conditions_long.csv has no column '", P_COLUMN, "'.", call. = FALSE)
conds <- unique(gstat$Condition)
is_pseudobulk <- all(is.na(long$Cluster))
if (any(long$ThreshCat != "")) note_warning("Expression categories (ThreshCat) are not separated in the panels; they would be mixed. Select one category before plotting.")

markers_plot <- if (length(MARKERS) == 0) unique(long$Marker) else MARKERS
if (length(setdiff(markers_plot, long$Marker)) > 0) {
  stop("MARKERS not in the Statistics folder: ", paste(setdiff(markers_plot, long$Marker), collapse = ", "), call. = FALSE)
}
clusters_all <- unique(long$Cluster)
if (!is_pseudobulk) clusters_all <- clusters_all[order(suppressWarnings(as.numeric(clusters_all)), clusters_all)]
clusters_plot <- if (is_pseudobulk) NA_character_ else if (length(CLUSTERS) == 0) clusters_all else as.character(CLUSTERS)
if (!is_pseudobulk && length(setdiff(clusters_plot, clusters_all)) > 0) {
  stop("CLUSTERS not in the Statistics folder: ", paste(setdiff(clusters_plot, clusters_all), collapse = ", "), call. = FALSE)
}
group_colors <- stats::setNames(viridis::viridis(n = length(conds), option = COLOR_OPTION), conds)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
p_is_q <- P_COLUMN == "q_BH_all_markers"

# Axis breaks as in the published script
nice_breaks <- function(mx) {
  by <- if (mx < 2) 0.2 else if (mx < 5) 0.5 else if (mx < 10) 1 else if (mx < 15) 2 else if (mx < 20) 5 else NA
  if (is.na(by)) seq(0, ceiling(mx * 1.5)) else seq(0, ceiling(mx * 1.5), by = by)
}



#### One panel ##################################################################
plotted_list <- list(); panels_by_marker <- list(); n_stars <- 0L
make_panel <- function(marker, cluster) {
  sel <- function(d) if (is.na(cluster)) d[Marker == marker & is.na(Cluster)] else d[Marker == marker & Cluster == cluster]
  gs <- sel(gstat); lg <- sel(long)
  pts <- if (is.na(cluster)) vals[is.na(Cluster)] else vals[Cluster == cluster]
  pts <- data.table::data.table(Condition = pts$Condition, Value = pts[[marker]])[is.finite(Value)]
  if (nrow(gs) == 0 || nrow(pts) == 0) return(NULL)
  gs[, Condition := factor(Condition, levels = conds)]
  pts[, Condition := factor(Condition, levels = conds)]
  lg <- lg[order(Pair)]
  lg[, P_Used := get(P_COLUMN)]
  lg[, Suppressed_Low_N := MIN_N_FOR_STARS > 0 & (is.na(N_Group1) | is.na(N_Group2) | pmin(N_Group1, N_Group2) < MIN_N_FOR_STARS)]
  lg[, Stars := ifelse(Suppressed_Low_N, "", stars(P_Used))]
  lg[, Label := ifelse(Suppressed_Low_N | is.na(P_Used), "", ifelse(Stars == "", if (SHOW_NS) "ns" else "", Stars))]
  n_stars <<- n_stars + sum(lg$Stars != "")
  maxv <- max(pts$Value, na.rm = TRUE)
  y1 <- 1.1 * maxv
  lg[, y_position := y1 + (Pair - 1) * STEP_INCREASE * maxv]
  top_needed <- max(lg$y_position) + 0.12 * maxv
  y_top <- max(round(maxv * 1.5, 1), top_needed)
  ann <- data.frame(start = lg$Group1, end = lg$Group2, label = lg$Label, y = lg$y_position, stringsAsFactors = FALSE)
  title_txt <- if (is.na(cluster)) paste0(PANEL_TITLE_PREFIX, "Pseudobulk - ", marker) else
    paste0(PANEL_TITLE_PREFIX, if (!is.null(names(CLUSTER_LABELS)) && cluster %in% names(CLUSTER_LABELS)) CLUSTER_LABELS[[cluster]] else paste0("Cluster ", cluster))
  p <- ggplot2::ggplot(pts, ggplot2::aes(x = Condition, y = Value)) +
    ggplot2::geom_col(data = gs, ggplot2::aes(x = Condition, y = Mean, fill = Condition), color = "black", alpha = 0.5, width = 0.9) +
    ggplot2::geom_point(ggplot2::aes(fill = Condition), color = "black", size = 2, shape = 21, stroke = 0.5) +
    ggplot2::geom_errorbar(data = gs[is.finite(SD)], ggplot2::aes(x = Condition, ymin = Mean - SD, ymax = Mean + SD),
                           width = 0.2, color = "black", inherit.aes = FALSE) +
    ggsignif::geom_signif(data = ann, ggplot2::aes(xmin = start, xmax = end, annotations = label, y_position = y),
                          manual = TRUE, color = "black", textsize = 3.88, vjust = 0.4, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(breaks = nice_breaks(maxv), limits = c(0, y_top)) +
    ggplot2::labs(title = title_txt, x = "Treatment", y = paste0(marker, ": FC relative to NoTreat 0 hr"), fill = "Treatment") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  key <- paste(marker, cluster)
  plotted_list[[key]] <<- cbind(
    data.table::data.table(Marker = marker, Cluster = cluster, P_Column_Drawn = P_COLUMN),
    lg[, c("Pair", "Group1", "Group2", "N_Group1", "N_Group2", "Mean_Group1", "Mean_Group2", "Omnibus_Test", "Omnibus_p",
           "PostHoc", "p_unadjusted", "p_Holm", "p_primary", "q_BH_all_markers", "Family_Size_BH", "p_Tukey",
           "Stars", "Label", "Suppressed_Low_N"), with = FALSE])
  p
}

if (is_pseudobulk) {
  panels <- lapply(markers_plot, function(m) make_panel(m, NA_character_))
  panels <- Filter(Negate(is.null), panels)
  figures <- list(list(name = "PseudoBulk", title = "PseudoBulk Plots", panels = panels))
} else {
  cl_tag <- if (length(clusters_plot) <= 6) paste(clusters_plot, collapse = "-") else paste0(length(clusters_plot), "clusters")
  figures <- lapply(markers_plot, function(m) {
    panels <- Filter(Negate(is.null), lapply(clusters_plot, function(cl) make_panel(m, cl)))
    list(name = paste0("Clusters_", cl_tag, "_", sub("_[^_]*$", "", m)),
         title = paste0("Cluster", paste(clusters_plot, collapse = ","), " Plots"), panels = panels)
  })
}
figures <- Filter(function(f) length(f$panels) > 0, figures)
if (length(figures) == 0) stop("Nothing could be drawn (no data for the chosen markers/clusters).", call. = FALSE)
plotted_data <- data.table::rbindlist(plotted_list, fill = TRUE)



#### Legend text ################################################################
pc_txt <- switch(P_COLUMN,
  p_primary = paste0("p-values adjusted for the three pairwise comparisons of each marker (",
                     if (all(plotted_data$PostHoc == "welch_holm")) "Welch's unpaired two-sided t-tests with Holm adjustment" else
                       "see column PostHoc in Plotted_data.csv for the test used per marker", ")"),
  p_Tukey = "Tukey HSD p-values from a one-way ANOVA (assumes equal spread of the groups)",
  p_Holm = "p-values of Welch's t-tests (unpaired, two-sided) with Holm adjustment for the three pairwise comparisons of each marker",
  p_unadjusted = "unadjusted p-values of Welch's unpaired two-sided t-tests (no correction for multiple comparisons)",
  q_BH_all_markers = paste0("Benjamini-Hochberg adjusted p-values (false discovery rate) over all marker x pair comparisons of the ",
                            "cluster/category (family size ", paste(sort(unique(plotted_data$Family_Size_BH)), collapse = "/"), "), ",
                            "calculated from unadjusted Welch's t-tests"))
n_by <- range(c(plotted_data$N_Group1, plotted_data$N_Group2), na.rm = TRUE)
legend_text <- paste0(
  "Bars show the mean and error bars the SD of the samples in each condition; dots are individual samples (n = ",
  if (n_by[1] == n_by[2]) n_by[1] else paste0(n_by[1], "-", n_by[2]), " per condition). Brackets compare conditions; stars give ",
  pc_txt, ": * ", if (p_is_q) "q" else "p", " < 0.05, ** < 0.01, *** < 0.001; ns = not significant. ",
  "Every p-value of every drawn comparison (unadjusted, Holm, BH and the one drawn) is listed in Plotted_data.csv.",
  if (any(plotted_data$Suppressed_Low_N)) paste0(" Comparisons with a group of fewer than ", MIN_N_FOR_STARS, " samples are drawn without statistics.") else "")



#### Save ######################################################################
Stage_WD <- make_stage_dir(stat_dir, paste0("Bargraphs_", P_COLUMN, if (nzchar(OUTPUT_LABEL)) paste0("_", OUTPUT_LABEL) else ""))
out_files <- character(0)
for (fg in figures) {
  panel_list <- unname(fg$panels)
  ncol_f <- min(PANEL_COLUMNS, length(panel_list)); nrow_f <- ceiling(length(panel_list) / ncol_f)
  n_slots <- ncol_f * nrow_f
  if (length(panel_list) < n_slots) panel_list <- c(panel_list, replicate(n_slots - length(panel_list), ggplot2::ggplot() + ggplot2::theme_void(), simplify = FALSE))
  grid_panels <- gridExtra::arrangeGrob(grobs = panel_list, nrow = nrow_f, ncol = ncol_f)
  title_grob <- grid::textGrob(fg$title, gp = grid::gpar(fontsize = 18, fontface = "bold"))
  full_plot <- gridExtra::arrangeGrob(title_grob, grid_panels, ncol = 1, heights = c(0.05, 0.95))
  base_name <- paste0("BarPlot_", if (is_pseudobulk) "PseudoBulk_" else "", fg$name, "-Grid")
  for (fmt in OUTPUT_FORMATS) {
    f <- file.path(Stage_WD, paste0(base_name, ".", fmt))
    ok <- tryCatch({
      if (fmt == "pdf") {
        tmp_pdf <- tempfile(fileext = ".pdf")    # R's pdf() cuts paths longer than 511 characters: write short, then copy
        ggplot2::ggsave(tmp_pdf, plot = full_plot, width = PANEL_SIZE_IN * ncol_f, height = PANEL_SIZE_IN * nrow_f, units = "in",
                        device = "pdf", useDingbats = FALSE, limitsize = FALSE)
        if (!file.copy(tmp_pdf, f, overwrite = TRUE)) stop("could not copy the PDF to ", f)
        unlink(tmp_pdf)
      } else if (has_ragg) {
        ggplot2::ggsave(f, plot = full_plot, width = PANEL_SIZE_IN * ncol_f, height = PANEL_SIZE_IN * nrow_f, units = "in",
                        dpi = PNG_DPI, device = ragg::agg_png, limitsize = FALSE)
      } else {
        ggplot2::ggsave(f, plot = full_plot, width = PANEL_SIZE_IN * ncol_f, height = PANEL_SIZE_IN * nrow_f, units = "in",
                        dpi = PNG_DPI, limitsize = FALSE)
      }
      TRUE
    }, error = function(e) { note_warning(paste0("Could not write ", basename(f), ": ", conditionMessage(e))); FALSE })
    if (ok) out_files <- c(out_files, f)
  }
}
if (length(out_files) == 0) cat("[ERROR] No figure was written; see the warnings above.\n")
f_pd  <- file.path(Stage_WD, "Plotted_data.csv"); data.table::fwrite(plotted_data, f_pd)
f_leg <- file.path(Stage_WD, "Figure_legend.txt")
writeLines(c(strwrap(legend_text, width = 100), "", paste0("Statistics folder: ", basename(stat_dir)),
             paste0("P-value column: ", P_COLUMN)), f_leg)
out_files <- c(out_files, f_pd, f_leg)
cat("[INFO] Written to:\n       ", Stage_WD, "\n", sep = "")



#### Run record ################################################################
Inputs <- list(
  Script_Version = SCRIPT_VERSION,
  Completed_At   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  INPUT_FOLDER   = stat_dir,
  OUTPUT_FOLDER  = Stage_WD,
  Settings = list(P_COLUMN = P_COLUMN, MARKERS = as.list(markers_plot), CLUSTERS = as.list(as.character(clusters_plot)),
                  MIN_N_FOR_STARS = MIN_N_FOR_STARS, SHOW_NS = SHOW_NS, STEP_INCREASE = STEP_INCREASE,
                  COLOR_OPTION = COLOR_OPTION, OUTPUT_FORMATS = as.list(OUTPUT_FORMATS), PNG_DPI = PNG_DPI,
                  PANEL_COLUMNS = PANEL_COLUMNS, PANEL_SIZE_IN = PANEL_SIZE_IN),
  Statistics_Folder = basename(stat_dir),
  Panels = sum(vapply(figures, function(f) length(f$panels), integer(1))),
  Figures = length(figures),
  Starred_Comparisons = n_stars,
  Legend_Text = legend_text,
  PNG_Device = if ("png" %in% OUTPUT_FORMATS) (if (has_ragg) "ragg::agg_png" else "grDevices::png") else "not used",
  Random_Seed = "not used in this script",
  Warnings = as.list(warnings_log),
  Output_MD5 = file_md5(out_files),
  Session = session_info_list(c("data.table", "jsonlite", "ggplot2", "ggsignif", "gridExtra", "viridis")))
record_params(MODULE_BASE_NAME, Inputs, Stage_WD, parent_dir = stat_dir)

cat("\n========================================================\n")
cat("[DONE] ", Inputs$Panels, " panels in ", length(figures), " figure(s), ", n_stars, " starred comparisons (p-value column: ", P_COLUMN, ")\n", sep = "")
cat("       Warnings: ", length(warnings_log), "\n", sep = "")
cat("       Legend text: Figure_legend.txt;  every p-value: Plotted_data.csv\n")
cat("========================================================\n")
