################################################################################
# 01_Palette.R   (shared by all plotting scripts; you do not run this file)
################################################################################
# Layer:    Helpers
# Author:   Jonathon Sewell (with Claude)
# Created:  2026-10-28
# Version:  1.2   (PALETTE_HELPERS_VERSION below)
#
# Purpose:
# * ONE place for the colours of every figure. The colours themselves are in Palette_Colors.csv (same folder);
#   change a colour THERE (one line per cluster / treatment / category), not in the plotting scripts.
# * Each cluster belongs to a "co-plotted group" (column CoPlotGroup): clusters that are drawn together in one
#   line graph. Inside a group the colours are chosen to be told apart, also by colour-blind readers, and the
#   LINE TYPE also differs (solid / long dash / dash / dots / dot-dash), so lines remain separable
#   even where two colours are close. A cluster always gets the same colour AND line type, in every figure.
# * BACKGROUND: "white", "transparent", "published" (the grey panel of the original figures) or any colour
#   name / hex code. Colours whose contrast with the background is too low are lightened or darkened
#   automatically (see MIN_RATIO_*), so a figure stays legible on a dark slide or a transparent figure.
#
# ==============================================================================
# CHANGE LOG
# ==============================================================================
# v1.2 [RESULT: how the mean points look only] POINT SHAPES (PALETTE_HELPERS_VERSION 3). palette_clusters() and palette_treatments() also
#    return $shapes (R point codes). Clusters: the shape follows the cluster's position in its co-plotted group, exactly like the line
#    type (a group of one = round dot). Treatments: round dots, unless Palette_Colors.csv has a column "Shape" with a value for the
#    treatment (optional; the column does not have to exist). palette_shapes_custom() checks a shape list given in the plotting script.
#    No colour, line type, statistic or number changes.
# v1.1 FIX: on a dark or transparent-on-black background the first version raised every too-dark colour to exactly the
#    minimum contrast, so the steps of a gradient collapsed into (nearly) one colour (smallest colour difference 0.4-2.2).
#    Now a co-plotted group is adjusted TOGETHER: its colours keep their order of lightness and are re-spaced evenly between
#    the minimum contrast and a lighter top (dark background), or capped below the maximum lightness (light background).
#    Colours that already pass are left unchanged. On white nothing changes.
# v1 [RESULT: colours only] NEW. Replaces rainbow() cluster colours (neighbouring clusters had neighbouring hues,
#    and yellows) and the per-script colour rules. Treatment colours NoTreat, ContDepr, CompMed, BDNF are unchanged
#    (BDNF+K252a = the former magma orange). No statistic, no number and no layout is affected.
################################################################################

PALETTE_HELPERS_VERSION <- 3

# Minimum contrast ratio (WCAG relative luminance) between a line colour and the background it is drawn on.
#   Light backgrounds: 2.0  (the primary blue gradient starts at about 2.1: light blue for the progenitors;
#                            a higher floor would crowd its steps together, see Decisions_Log.md 2026-10-28)
#   Dark backgrounds:  3.0  (the usual rule for graphical objects)
MIN_RATIO_LIGHT_BG <- 1.0
MIN_RATIO_DARK_BG  <- 3.0

# Line types used inside a co-plotted group (in the order of the group in Palette_Colors.csv)
PALETTE_LINETYPES <- c("solid", "longdash", "dashed", "dotted", "dotdash")

# Point shapes of the mean points, used the same way as the line types (R point codes):
#   19 = filled circle (the default point of every figure before v1.2), 17 = filled triangle, 15 = filled square,
#   18 = filled diamond, 4 = cross (x; the only one that is not filled, kept as the 5th so shapes and line types cycle together)
PALETTE_SHAPES <- c(19L, 17L, 15L, 18L, 4L)

# Checks a shape list typed in a plotting script (POINT_SHAPES): whole numbers 0-25, one per group (or one value for all).
# Returns the shapes named by group, or stops with a message.
palette_shapes_custom <- function(shapes, group_levels, name = "POINT_SHAPES") {
  v <- suppressWarnings(as.integer(unlist(shapes)))
  if (anyNA(v) || any(v < 0L | v > 25L)) stop(name, " must be \"auto\", \"none\" or whole numbers from 0 to 25 (R point codes), e.g. c(19, 17, 15).", call. = FALSE)
  if (length(v) == 1L) v <- rep(v, length(group_levels))
  if (length(v) != length(group_levels)) stop(name, " has ", length(v), " value(s) but the figure has ", length(group_levels),
                                              " groups (", paste(group_levels, collapse = ", "), "). Give one per group, in plotting order.", call. = FALSE)
  if (any(v >= 21L)) stop(name, ": shapes 21-25 need a fill colour and are not supported; use 15, 17, 18, 19 or 4.", call. = FALSE)
  stats::setNames(v, group_levels)
}


#### Reading the palette #######################################################
palette_file <- function() {
  candidates <- c(file.path(HELPERS_DIR, "Palette_Colors.csv"),
                  file.path(HELPERS_DIR, "..", "Palette", "Palette_Colors.csv"))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) {
    stop("Palette_Colors.csv was not found next to 01_Palette.R (folder: ", HELPERS_DIR, ").", call. = FALSE)
  }
  normalizePath(candidates[1], winslash = "/", mustWork = TRUE)
}

palette_read <- function() {
  pal <- utils::read.csv(palette_file(), stringsAsFactors = FALSE, colClasses = "character")
  need <- c("Set", "ID", "Label", "Color", "CoPlotGroup")
  if (!all(need %in% names(pal))) stop("Palette_Colors.csv must have the columns: ", paste(need, collapse = ", "), call. = FALSE)
  pal
}


#### Colour maths ##############################################################
palette_luminance <- function(col) {
  rgb <- grDevices::col2rgb(col) / 255
  lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  as.numeric(0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ])
}

# Contrast ratio between two colours (1 = identical, 21 = black on white); vectorised over 'a'
palette_contrast <- function(a, b) {
  la <- palette_luminance(a); lb <- palette_luminance(b)
  (pmax(la, lb) + 0.05) / (pmin(la, lb) + 0.05)
}

palette_hex <- function(rgb01) grDevices::rgb(rgb01[1], rgb01[2], rgb01[3])

# Lighten (dark background) or darken (light background) a colour step by step until its contrast with
# 'ref' reaches min_ratio. A colour that already passes is returned unchanged.
palette_lift_one <- function(col, ref, min_ratio) {
  if (palette_contrast(col, ref) >= min_ratio) return(col)
  ref_dark <- palette_luminance(ref) < 0.18
  target <- if (ref_dark) c(1, 1, 1) else c(0, 0, 0)
  base <- as.numeric(grDevices::col2rgb(col)) / 255
  for (t in seq(0, 1, by = 0.01)) {
    cand <- palette_hex(base + (target - base) * t)
    if (palette_contrast(cand, ref) >= min_ratio) return(cand)
  }
  palette_hex(target)
}


# Move a colour to a given relative luminance by blending it towards white (to lighten) or black (to darken)
palette_set_luminance <- function(col, y_target) {
  base <- as.numeric(grDevices::col2rgb(col)) / 255
  y0 <- palette_luminance(col)
  if (abs(y0 - y_target) < 1e-4) return(col)
  up <- y_target > y0
  target <- if (up) c(1, 1, 1) else c(0, 0, 0)
  lo <- 0; hi <- 1
  for (i in 1:40) {
    mid <- (lo + hi) / 2
    y <- palette_luminance(palette_hex(base + (target - base) * mid))
    if ((up && y < y_target) || (!up && y > y_target)) lo <- mid else hi <- mid
  }
  palette_hex(base + (target - base) * hi)
}

# Adjust the colours of ONE co-plotted group for a background (reference colour 'ref') so that every colour reaches
# 'min_ratio'. Colours that all pass are returned unchanged. Otherwise the group is re-spaced as a whole (same order of
# lightness, equal contrast steps), so a gradient stays a gradient.
palette_adjust_group <- function(cols, ref, min_ratio) {
  n <- length(cols)
  if (n == 0 || all(palette_contrast(cols, ref) >= min_ratio)) return(cols)
  y <- palette_luminance(cols); yref <- palette_luminance(ref)
  dark <- yref < 0.18
  if (n == 1) {
    y_limit <- if (dark) min_ratio * (yref + 0.05) - 0.05 else (yref + 0.05) / min_ratio - 0.05
    return(palette_set_luminance(cols, if (dark) max(y, y_limit) else min(y, y_limit)))
  }
  if (dark) {
    y_floor <- min_ratio * (yref + 0.05) - 0.05
    bottom <- y_floor
    top <- max(max(y), min(1, y_floor * (1 + 0.8 * (n - 1))))
  } else {
    y_ceil <- (yref + 0.05) / min_ratio - 0.05
    top <- min(max(y), y_ceil)
    bottom <- min(min(y), max(0, (y_ceil + 0.05) / (1 + 0.8 * (n - 1)) - 0.05))
  }
  targets <- exp(seq(log(bottom + 0.05), log(top + 0.05), length.out = n)) - 0.05    # equal contrast steps
  out <- cols
  ord <- order(y)                                                                  # keep the order of lightness
  for (k in seq_len(n)) out[ord[k]] <- palette_set_luminance(cols[ord[k]], targets[k])
  out
}

# Adjust a whole set of colours: 'groups' says which colours are adjusted together (same value = same group)
palette_adjust <- function(cols, groups, ref, min_ratio) {
  out <- cols
  for (g in unique(groups)) {
    idx <- which(groups == g)
    out[idx] <- palette_adjust_group(cols[idx], ref, min_ratio)
  }
  out
}


#### Background ################################################################
# background: "white", "transparent", "published" or any colour.
# transparent_reference: the colour a transparent figure will be placed on ("white" or "black"); it decides the
#   contrast check and the colour of text and axes.
palette_background <- function(background = "white", transparent_reference = "white") {
  background <- tolower(trimws(background))
  if (background == "published") {
    page <- "white"; panel <- "grey92"; ref <- "grey92"; transparent <- FALSE
  } else if (background == "transparent") {
    page <- "transparent"; panel <- "transparent"; ref <- transparent_reference; transparent <- TRUE
  } else {
    ok <- !inherits(try(grDevices::col2rgb(background), silent = TRUE), "try-error")
    if (!ok) stop("BACKGROUND = \"", background, "\" is not a colour R knows. Use \"white\", \"transparent\", \"published\", ",
                  "a colour name or a hex code such as \"#1E1E1E\".", call. = FALSE)
    page <- background; panel <- background; ref <- background; transparent <- FALSE
  }
  dark <- palette_luminance(ref) < 0.18
  list(name = background, page = page, panel = panel, reference = ref, transparent = transparent,
       dark = dark, ink = if (dark) "white" else "black",
       grid = grDevices::adjustcolor(if (dark) "white" else "black", alpha.f = 0.15),
       min_ratio = if (dark) MIN_RATIO_DARK_BG else MIN_RATIO_LIGHT_BG)
}

# ggplot2 theme for the chosen background (NULL for "published": the default grey theme stays as it was)
palette_theme <- function(bg) {
  if (identical(bg$name, "published")) return(NULL)
  fill <- if (bg$transparent) NA else bg$panel
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = fill, colour = NA),
    plot.background  = ggplot2::element_rect(fill = fill, colour = NA),
    panel.grid.major = ggplot2::element_line(colour = bg$grid),
    panel.grid.minor = ggplot2::element_line(colour = bg$grid, linewidth = 0.25),
    panel.border     = ggplot2::element_rect(fill = NA, colour = bg$ink, linewidth = 0.4),
    axis.text        = ggplot2::element_text(colour = bg$ink),
    axis.title       = ggplot2::element_text(colour = bg$ink),
    axis.ticks       = ggplot2::element_line(colour = bg$ink),
    plot.title       = ggplot2::element_text(colour = bg$ink, hjust = 0.5, face = "bold"),
    plot.subtitle    = ggplot2::element_text(colour = bg$ink, hjust = 0.5),
    legend.background = ggplot2::element_rect(fill = fill, colour = NA),
    legend.key       = ggplot2::element_rect(fill = fill, colour = NA),
    legend.text      = ggplot2::element_text(colour = bg$ink),
    legend.title     = ggplot2::element_text(colour = bg$ink))
}

# 'bg' argument of ggsave(): the page colour ("transparent" gives a figure without a background)
palette_save_bg <- function(bg) bg$page


#### Colours for the plots #####################################################
# Clusters: a list with colours and line types, named by cluster. Unknown cluster numbers (the palette has 1-20)
# get a neutral grey and a warning is returned in $missing.
#   merged = TRUE: the clusters are the 4 merged groups of the merged-cluster tables (R1_Merge), not R1 clusters
palette_clusters <- function(ids, bg, merged = FALSE) {
  pal <- palette_read()
  set_name <- if (merged) "MergedGroup" else "Cluster"
  p <- pal[pal$Set == set_name, , drop = FALSE]
  ids_chr <- as.character(ids)
  idx <- match(ids_chr, p$ID)
  cols <- ifelse(is.na(idx), "#808080", p$Color[idx])
  # line type = position of the cluster inside its co-plotted group (solid for a group of one)
  lty <- rep("solid", length(ids_chr))
  shp <- rep(PALETTE_SHAPES[1], length(ids_chr))
  if (!merged) {
    for (k in seq_along(ids_chr)) {
      if (is.na(idx[k])) next
      grp <- p$ID[p$CoPlotGroup == p$CoPlotGroup[idx[k]]]
      pos <- match(ids_chr[k], grp)
      lty[k] <- PALETTE_LINETYPES[((pos - 1) %% length(PALETTE_LINETYPES)) + 1]
      shp[k] <- PALETTE_SHAPES[((pos - 1) %% length(PALETTE_SHAPES)) + 1]
    }
  }
  # the whole group of the palette is adjusted together (so a figure that shows only some of its clusters still gets the
  # same colours as one that shows all of them)
  p$Adjusted <- palette_adjust(p$Color, p$CoPlotGroup, bg$reference, bg$min_ratio)
  lifted <- ifelse(is.na(idx), cols, p$Adjusted[idx])
  list(colors = stats::setNames(lifted, ids_chr), linetypes = stats::setNames(lty, ids_chr),
       shapes = stats::setNames(shp, ids_chr),
       lifted = ids_chr[lifted != cols], missing = ids_chr[is.na(idx)])
}

# Treatments: named colours, or NULL when a treatment is not in the palette (the script then keeps its old rule)
palette_treatments <- function(group_levels, bg) {
  pal <- palette_read()
  p <- pal[pal$Set == "Treatment", , drop = FALSE]
  idx <- match(group_levels, p$ID)
  if (anyNA(idx)) return(NULL)
  cols <- p$Color[idx]
  lifted <- palette_adjust(cols, p$ID[idx], bg$reference, bg$min_ratio)       # one group per treatment: each colour on its own
  shp <- rep(PALETTE_SHAPES[1], length(group_levels))
  if ("Shape" %in% names(p)) {                       # optional column: one R point code per treatment, blank = round dot
    sv <- suppressWarnings(as.integer(p$Shape[idx]))
    shp[!is.na(sv)] <- sv[!is.na(sv)]
  }
  list(colors = stats::setNames(lifted, group_levels),
       linetypes = stats::setNames(rep("solid", length(group_levels)), group_levels),
       shapes = stats::setNames(shp, group_levels),
       lifted = group_levels[lifted != cols], missing = character(0))
}
