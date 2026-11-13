################################################################################
# ColourMap.R  --  draws the palette so you can check the colours before (re)making figures
################################################################################
# Layer:    Visualization / Palette (no data are read; nothing here changes a figure)
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-28
# Version:  1.1
#
# Purpose:
# * One line per colour, labelled with what it stands for (cluster, treatment, category ...), drawn in its line type on
#   a rectangle of the background colour, with the hex code and the contrast ratio with that background.
# * A second page shows every co-plotted group as seen with normal vision and with three kinds of colour blindness
#   (simulated: deuteranopia, protanopia, tritanopia).
# * Palette_Check.csv: for each co-plotted group, the smallest contrast between neighbouring lines, the smallest colour
#   difference (OKLab, x100) between any two lines of the group under every vision type, and the contrast with the background.
#
# How to use:
# 1. Choose BACKGROUND below (same words as in 01_PlotLinegraphs.R). 2. Source the script.
# 3. Look at ColourMap_<background>.png / .pdf and Palette_Check.csv (written next to this script).
# Colours are changed in Helpers/Palette_Colors.csv (one line per colour), then run this script again.
#
# ==============================================================================
# CHANGE LOG
# ==============================================================================
# v1.1 colours that are too faint for the background are adjusted per co-plotted group (needs 01_Palette.R version 2).
# v1 NEW (2026-10-28).
################################################################################

#### USER SETTINGS -- edit this section only ###################################
BACKGROUND     <- "white"      # "white", "transparent", "published" (grey panel) or any colour, e.g. "#1E1E1E"
TRANSPARENT_ON <- "white"      # only for "transparent": the colour it will be placed on, "white" or "black"
SHOW_LINE_TYPES <- TRUE        # draw each sample line in the line type it gets in line graphs
PNG_DPI        <- 200

#### Fixed names -- do not change ###############################################
SCRIPT_NAME <- "ColourMap.R"

#### Find the helper files (do not edit this block) ############################
local({
  start_dir <- tryCatch({
    p <- rstudioapi::getSourceEditorContext()$path
    if (is.character(p) && nzchar(p)) dirname(p) else getwd()
  }, error = function(e) getwd())
  search_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  pal_helper <- NA_character_
  for (i in 1:5) {
    cand <- c(file.path(search_dir, "Helpers", "01_Palette.R"), file.path(search_dir, "01_Palette.R"))
    cand <- cand[file.exists(cand)]
    if (length(cand) > 0) { pal_helper <- cand[1]; break }
    search_dir <- dirname(search_dir)
  }
  if (is.na(pal_helper)) stop("Could not find Helpers/01_Palette.R. Save this script inside your ScriptsToUse folder ",
                              "(looked upward from: ", start_dir, ").", call. = FALSE)
  assign("HELPERS_DIR", dirname(pal_helper), envir = .GlobalEnv)
  assign("OUT_DIR", normalizePath(start_dir, winslash = "/", mustWork = FALSE), envir = .GlobalEnv)
  source(pal_helper, local = FALSE)
})

bg  <- palette_background(BACKGROUND, TRANSPARENT_ON)
pal <- palette_read()
cat("[START] Colour map; background = ", bg$name, " (reference colour ", bg$reference, "; minimum contrast ", bg$min_ratio, ")\n", sep = "")

#### Colour-blindness simulation and colour difference ########################
MACHADO <- list(                       # Machado et al. 2009, severity 1.0, applied to linear RGB
  deutan = matrix(c(0.367322, 0.860646, -0.227968, 0.280085, 0.672501, 0.047413, -0.011820, 0.042940, 0.968881), 3, byrow = TRUE),
  protan = matrix(c(0.152286, 1.052583, -0.204868, 0.114503, 0.786281, 0.099216, -0.003882, -0.048116, 1.051998), 3, byrow = TRUE),
  tritan = matrix(c(1.255528, -0.076749, -0.178779, -0.078411, 0.930809, 0.147602, 0.004733, 0.691367, 0.303900), 3, byrow = TRUE))
to_lin   <- function(c01) ifelse(c01 <= 0.04045, c01 / 12.92, ((c01 + 0.055) / 1.055)^2.4)
from_lin <- function(l) ifelse(l <= 0.0031308, 12.92 * l, 1.055 * l^(1 / 2.4) - 0.055)
simulate_cvd <- function(cols, type) {
  if (type == "normal") return(unname(cols))
  vapply(cols, function(col) {
    lin <- to_lin(as.numeric(grDevices::col2rgb(col)) / 255)
    out <- pmin(pmax(as.numeric(MACHADO[[type]] %*% lin), 0), 1)
    grDevices::rgb(from_lin(out[1]), from_lin(out[2]), from_lin(out[3]))
  }, character(1), USE.NAMES = FALSE)
}
oklab <- function(col) {
  c <- to_lin(as.numeric(grDevices::col2rgb(col)) / 255)
  l <- (0.4122214708 * c[1] + 0.5363325363 * c[2] + 0.0514459929 * c[3])^(1 / 3)
  m <- (0.2119034982 * c[1] + 0.6806995451 * c[2] + 0.1073969566 * c[3])^(1 / 3)
  s <- (0.0883024619 * c[1] + 0.2817188376 * c[2] + 0.6299787005 * c[3])^(1 / 3)
  c(0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
    1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
    0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)
}
delta_e <- function(a, b) 100 * sqrt(sum((oklab(a) - oklab(b))^2))
worst_difference <- function(cols) {
  if (length(cols) < 2) return(NA_real_)
  w <- Inf
  for (type in c("normal", "deutan", "protan", "tritan")) {
    s <- simulate_cvd(cols, type)
    for (i in 1:(length(s) - 1)) for (j in (i + 1):length(s)) w <- min(w, delta_e(s[i], s[j]))
  }
  w
}

#### Colours after the automatic contrast adjustment ##########################
show <- pal[pal$Set %in% c("Cluster", "MergedGroup", "Treatment", "Category", "Tier", "Merge"), , drop = FALSE]
# co-plotted groups are adjusted together (clusters, merged groups, tiers); every other colour on its own (same rule as the plots)
adjust_group <- ifelse(show$Set %in% c("Cluster", "MergedGroup", "Tier"), paste(show$Set, show$CoPlotGroup), paste(show$Set, show$ID))
show$Final <- palette_adjust(show$Color, adjust_group, bg$reference, bg$min_ratio)
show$Contrast <- round(palette_contrast(show$Final, bg$reference), 2)
show$Adjusted <- show$Final != show$Color
grp_pos <- ave(seq_len(nrow(show)), paste(show$Set, show$CoPlotGroup), FUN = seq_along)
show$Linetype <- ifelse(show$Set == "Cluster" & SHOW_LINE_TYPES, PALETTE_LINETYPES[((grp_pos - 1) %% length(PALETTE_LINETYPES)) + 1], "solid")
show$Name <- paste0(show$Set, " ", show$ID, ": ", show$Label)

#### Page 1: one line per colour ###############################################
tag <- gsub("[^A-Za-z0-9]+", "", bg$name)
png_file <- file.path(OUT_DIR, paste0("ColourMap_", tag, ".png"))
pdf_file <- file.path(OUT_DIR, paste0("ColourMap_", tag, ".pdf"))
draw_map <- function() {
  n <- nrow(show)
  op <- graphics::par(mar = c(0, 0, 0, 0), bg = bg$page); on.exit(graphics::par(op))
  graphics::plot.new(); graphics::plot.window(xlim = c(0, 10), ylim = c(0, n + 2), xaxs = "i", yaxs = "i")
  graphics::rect(0, 0, 10, n + 2, col = bg$reference, border = NA)             # the background rectangle (reference colour)
  graphics::text(0.3, n + 1.1, paste0("Colour map  -  background: ", bg$name,
                 if (bg$transparent) paste0(" (shown on ", bg$reference, ")") else "", "   -   contrast ratio with the background"),
                 adj = c(0, 0.5), col = bg$ink, font = 2, cex = 0.9)
  for (k in seq_len(n)) {
    y <- n + 0.4 - k
    graphics::segments(0.3, y, 3.2, y, col = show$Final[k], lwd = 5, lty = show$Linetype[k], lend = 1)
    graphics::text(3.5, y, show$Name[k], adj = c(0, 0.5), col = bg$ink, cex = 0.75)
    graphics::text(9.9, y, paste0(show$Final[k], "  ", sprintf("%.1f", show$Contrast[k]), ":1",
                   if (show$Contrast[k] < bg$min_ratio) "  LOW" else "", if (show$Adjusted[k]) "  (adjusted)" else ""),
                   adj = c(1, 0.5), col = bg$ink, cex = 0.65, family = "mono")
  }
}
png_dev <- function(f, w, h) {
  if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(f, width = w, height = h, units = "in", res = PNG_DPI, background = bg$page)
  else grDevices::png(f, width = w, height = h, units = "in", res = PNG_DPI, bg = bg$page)
}
h_map <- 0.2 * (nrow(show) + 2) + 0.3
png_dev(png_file, 8.5, h_map); draw_map(); grDevices::dev.off()
grDevices::pdf(pdf_file, width = 8.5, height = h_map, bg = bg$page); draw_map(); grDevices::dev.off()

#### Page 2: co-plotted groups under colour-blindness simulation ###############
cl <- show[show$Set == "Cluster", , drop = FALSE]
types <- c("normal", "deutan", "protan", "tritan")
cvd_png <- file.path(OUT_DIR, paste0("ColourMap_CVD_", tag, ".png"))
png_dev(cvd_png, 9, 4)
op <- graphics::par(mar = c(0.2, 0.2, 0.2, 0.2), mfrow = c(4, 1), bg = bg$page)
for (type in types) {
  graphics::plot.new(); graphics::plot.window(xlim = c(0, nrow(cl)), ylim = c(0, 1.25), xaxs = "i", yaxs = "i")
  graphics::rect(0, 0, nrow(cl), 1.25, col = bg$reference, border = NA)
  sim <- simulate_cvd(cl$Final, type)
  graphics::rect(seq_len(nrow(cl)) - 0.95, 0.05, seq_len(nrow(cl)) - 0.05, 0.85, col = sim, border = NA)
  graphics::text(seq_len(nrow(cl)) - 0.5, 0.45, cl$ID,
                 col = ifelse(palette_luminance(sim) < 0.25, "white", "black"), cex = 0.7)
  graphics::text(0.05, 1.05, type, adj = c(0, 0.5), col = bg$ink, cex = 0.8, font = 2)
}
graphics::par(op); grDevices::dev.off()

#### Check table per co-plotted group #########################################
grp_rows <- split(cl, cl$CoPlotGroup)
check <- do.call(rbind, lapply(names(grp_rows), function(g) {
  x <- grp_rows[[g]]
  data.frame(CoPlotGroup = g, Clusters = paste(x$ID, collapse = ","), N = nrow(x),
             Min_Contrast_With_Background = min(x$Contrast),
             Min_Neighbour_Contrast = if (nrow(x) > 1) round(min(palette_contrast(x$Final[-nrow(x)], x$Final[-1])), 2) else NA_real_,
             Worst_Colour_Difference_Any_Vision = round(worst_difference(x$Final), 1),
             Colours_Adjusted_For_Contrast = paste(x$ID[x$Adjusted], collapse = ","),
             stringsAsFactors = FALSE)
}))
utils::write.csv(check, file.path(OUT_DIR, paste0("Palette_Check_", tag, ".csv")), row.names = FALSE)
cat("\nCo-plotted groups (colour difference: OKLab x100, smallest over normal/deutan/protan/tritan; about 8 or more is easy to tell apart):\n")
# Write summary file
utils::write.csv(show, file.path(OUT_DIR, paste0("Palette_Summary_", tag, ".csv")), row.names = FALSE)
print(check, row.names = FALSE)
cat("\n[DONE] Written next to this script:\n  ", basename(png_file), "\n  ", basename(pdf_file), "\n  ", basename(cvd_png), "\n  ",
    paste0("Palette_Check_", tag, ".csv"), "\n", sep = "")
