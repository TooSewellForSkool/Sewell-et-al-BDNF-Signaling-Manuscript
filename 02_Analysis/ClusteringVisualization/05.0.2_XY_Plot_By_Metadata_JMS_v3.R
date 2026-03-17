
# Author: Jonathon Sewell, adapted from XY_Plot_By_Metadata.R (author Austin Keeler & Ashley Hirt - 15 July, 2020)
# Date: 2025-10-02
# Version: JMS_3
# Update:
# (1) Enable the ability to plot clusters


# Created by Austin Keeler & Ashley Hirt
# 15 July 2020
# Plot a single UMAP colored by a metadata variable

print("Start UMAP_ColorByMetadata")

rm(list = ls())
.libPaths( c( .libPaths(), "/project/zunderlab/R/4.0.3_Pipeline2") )

library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer) # Optional: depends on how color plot
library(viridis) # Optional: depends on how color plot



#### Input parameters ##########################################################
LAYOUT_ROUND <- 2
LAYOUT_GROUP <- 1 # Don't use group 0, not assigned XY coordinates

METADATA_WITHIN_METADATA <- TRUE # TRUE / FALSE
SELECTED_METADATA <- "Condition" # make sure this matches your metadata column name exactly
# "Expt" / "Condition" / "Treatment"
SELECTED_METADATA_VALUE <- c("BDNF_16hr") # See 'if (METADATA_TO_PLOT == "Condition")' line for options
# c("NoTreat", "BDNF")

CLUSTER_ROUND <- 1
CLUSTERS_MERGED <- FALSE # TRUE / FALSE
SELECTED_CLUSTER <- c("Cluster_R1") # c() if not filtering by cluster / "Cluster_R1" / "Cluster_R2"
PlotByIndCluster <- FALSE # TRUE / FALSE
SELECTED_CLUSTER_VALUE <- c(18)
# c(1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 15, 17, 18, 20) 

METADATA_TO_PLOT <- "TimePoint" # make sure this matches your metadata column name exactly
# "Condition" / "Treatment" / "TimePoint"

if (METADATA_TO_PLOT == "Condition") {
  METADATA_ORDER <-c("NoTreat_0hr",
                     "ContDepr_5min", "ContDepr_15min", "ContDepr_1hr", "ContDepr_4hr", "ContDepr_16hr",
                     "CompMed_5min", "CompMed_15min", "CompMed_1hr", "CompMed_4hr", "CompMed_16hr",
                     "BDNF_5min", "BDNF_15min", "BDNF_1hr", "BDNF_4hr", "BDNF_16hr",
                     "ContDepr_Test", "BDNFTest_1hr", "BDNF+K252aTest_1hr"
                     )
  } else if (METADATA_TO_PLOT == "Treatment") {
    METADATA_ORDER <- c("NoTreat", 
                        "ContDepr", "CompMed", "BDNF", 
                        "BDNFTest", "BDNF+K252aTest"
                        )
    } else if (METADATA_TO_PLOT == "TimePoint") {
      METADATA_ORDER <- c("0", "0.083","0.25","1","4","16"
                          #,"72"
                          )
      } else {
        METADATA_ORDER <- NULL
        }
# ^^ This needs to match with the values in metadata.csv in the METADATA_TO_PLOT column
# Otherwise, set METADATA_ORDER to NULL and metadata will be sorted alphabetically/numerically



METADATA_FILENAME <- "metadata.csv"
FILENUMS_FILENAME <- "filenums.csv"
LAYOUT_BASENAME <- "umap_xy_RX.csv"
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
FOLDER_BASENAME <- paste0("05.0.2-", # Script info
                          format(Sys.time(), "%Y-%m-%d"), # Date info 
                          " PlotByMetadata_VAR_XY")
output_dir <- sub("VAR", paste0("L", LAYOUT_ROUND),
                  FOLDER_BASENAME)

POINT_SIZE <- 0.1
LABEL_SIZE <- 5
METADATA_SCALE <- "continuous" # "continuous" or "discrete"
# CONTINUOUS_COLORSCALE <- "plasma" # "plasma", "viridis", etc.
OUTPUT_FORMAT <- "png" # "png" or "pdf"
OUTPUT_HEIGHT <- 7
OUTPUT_WIDTH <- 7

OmittedColor <- "lightgray"
ColorPaletteToUse <- "viridis" # viridis / RColorBrewer
ColorToUse <- c("D")
# With viridis, A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo
# With RColorBrewer, use display.brewer.all(colorblindFriendly = TRUE) to identify palette codes
# use c(# , name), where '#' is the number of colors from 'name' palette to be used for colorRampPalette interpolation
# Ex. c(8, "Set2") / c(12, "Paired") / c(8, "Dark2") / c(11, "RdYlBu") / c(11, "PiYG") / c(11, "BrBG")

OUTPUT_FILENAME <- paste0("Colored_by_", METADATA_TO_PLOT, ".png")

if (LAYOUT_ROUND > 1) {
  output_filename <- paste0("Colored_by_", METADATA_TO_PLOT, "_L",
                            LAYOUT_ROUND, "_G", LAYOUT_GROUP, ".png")
  layout_cols <- c("group", "umap_x", "umap_y")
} else {
  output_filename <- paste0("Colored_by_", METADATA_TO_PLOT, "_L",
                            LAYOUT_ROUND, ".png")
  layout_cols <- c("umap_x", "umap_y")
}

if (METADATA_WITHIN_METADATA) {
  output_filename <- sub(".png", 
                         paste0("_", paste0(SELECTED_METADATA_VALUE, collapse = ","), ".png"), 
                         output_filename)
}

if (CLUSTERS_MERGED) {
  output_filename <- sub(".png", 
                         "_merged.png", 
                         output_filename)
}


#### Read in necessary files ###################################################
# Read in metadata
metadata <- read.table(METADATA_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)

# Read in filenum labels
filenums <- read.table(FILENUMS_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)[,]

# Read UMAP layout, has two columns, x and y, where each row is a cell
layout_filename <- sub("RX", paste0("R", LAYOUT_ROUND), LAYOUT_BASENAME)
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)
colnames(layout_in) <- layout_cols



if (!is.null(SELECTED_CLUSTER)) { # Read in clusters if user selects a cluster round
  clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
  if (CLUSTERS_MERGED) {
    clusters_filename <- sub(".csv", "_merged.csv", clusters_filename)
  }
  clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                            check.names=FALSE)
  
  if (CLUSTER_ROUND > 1) { # Combine cluster rounds into clusters_in
    clusters_filename.1 <- sub("RX", paste0("R", CLUSTER_ROUND-1), CLUSTERS_BASENAME)
    clusters_in.1 <- read.table(clusters_filename.1, header=TRUE, sep=",",
                                check.names=FALSE)
    
    clusters_in <- cbind(clusters_in.1, clusters_in)
    }
}

# Generate a dataframe to build on
if (!is.null(SELECTED_CLUSTER)) {
  plot_df <- cbind(layout_in, 
                   clusters_in) # Added clusters_in
} else {
  plot_df <- layout_in
}












#### Organize input files ######################################################

if (!METADATA_TO_PLOT %in% colnames(metadata)) {
  # Get metadata order (important for plotting with discrete and especially continuous scales!)
  # If no order is specified by METADATA_ORDER, then calculate alphabetical order.
  # Otherwise, use order specified by METADATA_ORDER
  if (is.null(METADATA_ORDER)) { 
    unique_sorted <- as.character(sort(unique(plot_df[ , METADATA_TO_PLOT])))
    metadata_order <- 1:length(unique_sorted)
    names(metadata_order) <- as.character(sort(unique(plot_df[ , METADATA_TO_PLOT])))
    md_converted_for_plot <- as.character(plot_df[ , METADATA_TO_PLOT])
    } else {
      metadata_order <- 1:length(METADATA_ORDER)
      names(metadata_order) <- METADATA_ORDER
    }
  
  # Combine initialized plot_df with METADATA_TO_PLOT info 
  plot_df <- cbind(plot_df,
                    md_converted_for_plot)
  
  # This comes from Stack Overflow.  I don't really understand it,
  # but it's what let's us rename legend labels how we want them to look.
  # https://stackoverflow.com/questions/12075037/ggplot-legends-change-labels-order-and-title
  plot_df$md_level <- factor(plot_df$md_converted_for_plot,
                             levels = 1:length(metadata_order),
                              labels = names(metadata_order))
  } else {
    # Get metadata order (important for plotting with discrete and especially continuous scales!)
    # If no order is specified by METADATA_ORDER, then calculate alphabetical order.
    # Otherwise, use order specified by METADATA_ORDER
    if (is.null(METADATA_ORDER)) {
      unique_sorted <- sort(unique(as.character(metadata[ , METADATA_TO_PLOT])))
      metadata_order <- 1:length(unique_sorted)
      names(metadata_order) <- sort(unique(as.character(metadata[ , METADATA_TO_PLOT])))
      } else {
        metadata_order <- 1:length(METADATA_ORDER)
        names(metadata_order) <- METADATA_ORDER
        }
    
    # Make metadata key, which will be used to convert metadata and file_nums
    metadata_key <- as.character(metadata[ , METADATA_TO_PLOT])
    names(metadata_key) <- 1:length(metadata[ , METADATA_TO_PLOT])
    
    # Make filenum key, which will be used to convert File_Num to plotting level
    filenum_key <- metadata_order[metadata_key]
    names(filenum_key) <- 1:length(metadata[ , METADATA_TO_PLOT])
    
    # Convert filenames for each cell event to metadata values
    md_converted_for_plot <- filenum_key[filenums]
      
    # Combine initialized plot_df with METADATA_TO_PLOT info 
    plot_df <- cbind(plot_df,
                      md_converted_for_plot)
    
    # This comes from Stack Overflow.  I don't really understand it,
    # but it's what let's us rename legend labels how we want them to look.
    # https://stackoverflow.com/questions/12075037/ggplot-legends-change-labels-order-and-title
    plot_df$md_level <- factor(plot_df$md_converted_for_plot,
                               levels = 1:length(metadata_order),
                               labels = names(metadata_order))
    }



# 
# if (is.null(METADATA_ORDER)) {
#   if (!METADATA_TO_PLOT %in% colnames(metadata)) {
#     unique_sorted <- as.character(sort(unique(plot_df[ , METADATA_TO_PLOT])))
#     metadata_order <- 1:length(unique_sorted)
#     names(metadata_order) <- as.character(sort(unique(plot_df[ , METADATA_TO_PLOT])))
#     } else {
#       unique_sorted <- sort(unique(as.character(metadata[ , METADATA_TO_PLOT])))
#       metadata_order <- 1:length(unique_sorted)
#       names(metadata_order) <- sort(unique(as.character(metadata[ , METADATA_TO_PLOT])))
#       }
#   } else {
#   metadata_order <- 1:length(METADATA_ORDER)
#   names(metadata_order) <- METADATA_ORDER
# }
# 
# 
# # Make metadata key, which will be used to convert metadata and file_nums
# if (!METADATA_TO_PLOT %in% colnames(metadata)) {
#   metadata_key <- as.character(plot_df[ , METADATA_TO_PLOT])
#   names(metadata_key) <- 1:length(plot_df[ , METADATA_TO_PLOT])
#   } else {
#     metadata_key <- as.character(metadata[ , METADATA_TO_PLOT])
#     names(metadata_key) <- 1:length(metadata[ , METADATA_TO_PLOT])
#     }
# 
# # Make filenum key, which will be used to convert File_Num to plotting level
# filenum_key <- metadata_order[metadata_key]
# names(filenum_key) <- 1:length(metadata[ , METADATA_TO_PLOT])
# 
# # Convert filenames for each cell event to metadata values
# if (!METADATA_TO_PLOT %in% colnames(metadata)) {
#   md_converted_for_plot <- metadata_key
#   
#   
# } else {
#   # Make filenum key, which will be used to convert File_Num to plotting level
#   filenum_key <- metadata_order[metadata_key]
#   names(filenum_key) <- 1:length(metadata[ , METADATA_TO_PLOT])
#   
#   md_converted_for_plot <- filenum_key[filenums]
# }
# 
# plot_df <- cbind(plot_df,
#                  md_converted_for_plot)
# 
# 
# 
# # This comes from Stack Overflow.  I don't really understand it,
# # but it's what let's us rename legend labels how we want them to look.
# # https://stackoverflow.com/questions/12075037/ggplot-legends-change-labels-order-and-title
# plot_df$md_level <- factor(plot_df$md_converted_for_plot,
#                            levels = 1:length(metadata_order),
#                            labels = names(metadata_order))

# Generate second set of above values if plotting METADATA_WITHIN_METADATA
if (METADATA_WITHIN_METADATA) {
  if (!SELECTED_METADATA %in% colnames(metadata)) {
    # Get metadata order (important for plotting with discrete and especially continuous scales!)
    # If no order is specified by METADATA_ORDER, then calculate alphabetical order.
    # Otherwise, use order specified by METADATA_ORDER
  
      unique_sorted_2 <- as.character(sort(unique(plot_df[ , SELECTED_METADATA])))
      #metadata_order <- 1:length(unique_sorted)
      metadata_order_2 <- 1:length(unique_sorted_2)
      names(metadata_order_2) <- as.character(sort(unique(plot_df[ , SELECTED_METADATA])))
      md_converted_for_plot_2 <- as.character(plot_df[ , SELECTED_METADATA])
    
    
    # Combine initialized plot_df with METADATA_TO_PLOT info 
    plot_df <- cbind(plot_df,
                     md_converted_for_plot_2)
    
    # This comes from Stack Overflow.  I don't really understand it,
    # but it's what let's us rename legend labels how we want them to look.
    # https://stackoverflow.com/questions/12075037/ggplot-legends-change-labels-order-and-title
    plot_df$md_level_2 <- factor(plot_df$md_converted_for_plot_2,
                               levels = 1:length(metadata_order_2),
                               labels = names(metadata_order_2))
    
    } else {
      # Get metadata order (important for plotting with discrete and especially continuous scales!)
      # If no order is specified by METADATA_ORDER, then calculate alphabetical order.
      # Otherwise, use order specified by METADATA_ORDER
      unique_sorted_2 <- sort(unique(as.character(metadata[ , SELECTED_METADATA])))
      metadata_order_2 <- 1:length(unique_sorted_2)
      names(metadata_order_2) <- sort(unique(as.character(metadata[ , SELECTED_METADATA])))
      
      # Make metadata key, which will be used to convert metadata and file_nums
      metadata_key_2 <- as.character(metadata[ , SELECTED_METADATA])
      names(metadata_key_2) <- 1:length(metadata[ , SELECTED_METADATA])
    
      # Make filenum key, which will be used to convert File_Num to plotting level
      filenum_key_2 <- metadata_order_2[metadata_key_2]
      names(filenum_key_2) <- 1:length(metadata[ , SELECTED_METADATA])
    
      # Convert filenames for each cell event to metadata values
      md_converted_for_plot_2 <- filenum_key_2[filenums]
      
      # Combine initialized plot_df with METADATA_TO_PLOT info 
      plot_df <- cbind(plot_df,
                       md_converted_for_plot_2)
    
      # This comes from Stack Overflow.  I don't really understand it,
      # but it's what let's us rename legend labels how we want them to look.
      # https://stackoverflow.com/questions/12075037/ggplot-legends-change-labels-order-and-title
      plot_df$md_level_2 <- factor(plot_df$md_converted_for_plot_2,
                                   levels = 1:length(metadata_order_2),
                                   labels = names(metadata_order_2))
      }
  }




















# # Generate second set of above values if plotting METADATA_WITHIN_METADATA
#     metadata_key_2 <- as.character(metadata[ , SELECTED_METADATA])
#     names(metadata_key_2) <- 1:length(metadata[ , SELECTED_METADATA])
#     
#     metadata_order_2 <- 1:length(unique(metadata[ , SELECTED_METADATA]))
#     names(metadata_order_2) <- unique(metadata[ , SELECTED_METADATA])
#     
#     filenum_key_2 <- metadata_order_2[metadata_key_2]
#     names(filenum_key_2) <- 1:length(metadata[ , SELECTED_METADATA])
#     
#     md_converted_for_plot_2 <- filenum_key_2[filenums]
#     
#     plot_df <- cbind(plot_df, md_converted_for_plot_2)
#     
#     plot_df$md_level_2 <- factor(plot_df$md_converted_for_plot_2,
#                                  levels = 1:length(metadata_order_2),
#                                  labels = names(metadata_order_2))


    
    
    
    
if (LAYOUT_ROUND > 1) {
  plot_df <- plot_df[which(plot_df[ , "group"] == LAYOUT_GROUP), ]
  }

    
# Randomize order for plotting, to avoid cells at the end (often
# the last file concatenated) going on top of everything else
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)), ]
row.names(plot_df) <- NULL #remove old indices, so everything has a new order



#### Generate output folder (if necessary) and set as working directory ########
if (dir.exists(output_dir)) { # If output directory exists, set output directory as new working directory
  print("output_dir exists")
  setwd(output_dir)
} else { # Create output directory if it does not exists
  dir.create(output_dir)
  setwd(output_dir)
}

# Define the colors to be used
if (ColorPaletteToUse == "viridis") {
  colors <- viridis(n = length(unique(plot_df[["md_level"]])),
                    option = ColorToUse)
  
} else if (ColorPaletteToUse == "RColorBrewer") {
  colors <- colorRampPalette(color = brewer.pal(n = as.numeric(ColorToUse[1]),
                                                name = ColorToUse[2]))(length(unique(plot_df[["md_level"]])))
}

#### Generate and save plot ####################################################
Plotting_Function <- function(PLOT_CLUSTER) {
  
if (METADATA_WITHIN_METADATA) {
  plot_df <- plot_df %>% 
    mutate(grouped_colors = ifelse(test = plot_df[["md_level_2"]] %in% SELECTED_METADATA_VALUE &
                                     plot_df[ , sub("Cluster_", "", SELECTED_CLUSTER)] %in% PLOT_CLUSTER, # PLOT_CLUSTER = SELECTED_CLUSTER_VALUE
                                   yes = plot_df[["md_level"]], 
                                   no = OmittedColor))
  
  # color_mapping <- c(setNames(scales::hue_pal()(length(unique(plot_df[["md_level"]]))), sort(unique(plot_df[["md_converted_for_plot"]]))), "gray" = "gray")
  
  p <- ggplot(data = plot_df, 
              aes(x = umap_x, y = umap_y)) +
    
    geom_point(data = plot_df[plot_df$grouped_colors == OmittedColor, ],
               size = POINT_SIZE, color = OmittedColor, alpha = 0.5) +

    geom_point(data = plot_df[plot_df$grouped_colors != OmittedColor, ],
               aes(color = plot_df[plot_df$grouped_colors != OmittedColor, "md_level"]),
               size = POINT_SIZE, alpha = 0.75) + 
    
    labs(color = METADATA_TO_PLOT) +
    guides(color = guide_legend(override.aes = list(size = LABEL_SIZE))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          #legend.key.size = unit(3,"line")
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = NA))
  
  } else {
    p <- ggplot(data = plot_df, aes(x = umap_x, y = umap_y,
                                  group = md_level, color = md_level)) +
    geom_point(size = POINT_SIZE) +
    labs(color = METADATA_TO_PLOT) +
    guides(color = guide_legend(override.aes = list(size = LABEL_SIZE))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          #legend.key.size = unit(3,"line")
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = NA))
    }

if (METADATA_SCALE == "discrete") {
  p <- p + scale_color_discrete()
}
if (METADATA_SCALE == "continuous") {
  p <- p + scale_color_manual(values = colors)
}
  
  if (!is.null(SELECTED_CLUSTER)) {
    output_filename <- sub(".png", 
                           paste0("_", SELECTED_CLUSTER, 
                                  "-", paste0(PLOT_CLUSTER, collapse = ","), ".png"), 
                           output_filename)
  }

  ggsave(filename = output_filename, 
         plot = p, 
         device = OUTPUT_FORMAT,
         height = OUTPUT_HEIGHT, width = OUTPUT_WIDTH)
  
}


if (!is.null(SELECTED_CLUSTER) & PlotByIndCluster) {
  for (i in SELECTED_CLUSTER_VALUE) {
    Plotting_Function(i)
    }
  } else {
    Plotting_Function(SELECTED_CLUSTER_VALUE)
    }







#### Reset working directory ###################################################
setwd("..")

print("End UMAP_ColorByMetadata")
