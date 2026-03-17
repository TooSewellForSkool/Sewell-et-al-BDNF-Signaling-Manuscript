
# Author: Jonathon Sewell, adapted from XY_Plot_By_Metadata.R (author Austin Keeler & Ashley Hirt - 15 July, 2020)
# Date: 2025-08-22
# Version: JMS_2
# Update:
# (1) Enables coloring only user-defined clusters


# Corey Williams, University of Virginia
#15 Jul, 2019
#Plot colored by expression of markers

print("Start Plot_XY_By_Cluster")

rm(list = ls())
.libPaths( c( .libPaths(), "/project/zunderlab/R/4.0.3_Pipeline2") )

library(ggfortify)

print("libraries loaded")



## Variable input parameters #################################################
LAYOUT_ROUND <- 2
LAYOUT_GROUP <- 1 # Don't use group 0, not assigned XY coordinates

CLUSTER_ROUND <- 1
CLUSTERS_MERGED <- FALSE # TRUE / FALSE
CLUSTERS_TO_PLOT <- c(1)  
    # c() to color all / c(1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 15, 17, 18, 20)
CLUSTERS_TO_PLOT_RED <- TRUE # TRUE / FALSE

OUTPUT_DEVICE <- "png" # "png", "pdf", "jpeg", etc
POINT_SIZE <- 0.05 # 0.01 for all cells / 0.05 for subsets

OmittedColor <- "lightgray"
ColorPaletteToUse <- "RColorBrewer" # viridis / RColorBrewer / Rainbow
if (CLUSTER_ROUND == 1) {
  ColorPaletteToUse <- "Rainbow"
  }
ColorToUse <- c(12, "Paired")
# With viridis, A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo
# With RColorBrewer, use display.brewer.all(colorblindFriendly = TRUE) to identify palette codes
# use c(# , name), where '#' is the number of colors from 'name' palette to be used for colorRampPalette interpolation
# Ex. c(8, "Set2") / c(12, "Paired") / c(8, "Dark2") / c(11, "RdYlBu") / c(11, "PiYG") / c(11, "BrBG")



#### Fixed parameters ##########################################################
LAYOUT_BASENAME <- "umap_xy_RX.csv"
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
OUTPUT_BASENAME <- "ClustersVAR_xy"

if (LAYOUT_ROUND > 1) {
  output_basename <- sub("VAR", paste0("_L", LAYOUT_ROUND,
                                       "_G", LAYOUT_GROUP,
                                       "_C", CLUSTER_ROUND),
                         OUTPUT_BASENAME)
} else {
  output_basename <- sub("VAR", paste0("_L", LAYOUT_ROUND,
                                       "_C", CLUSTER_ROUND),
                         OUTPUT_BASENAME)
}

if (!is.null(CLUSTERS_TO_PLOT)) {
  output_basename <- paste0(output_basename,
                            "_Cluster-", paste0(CLUSTERS_TO_PLOT, collapse = ","))
}

if (CLUSTERS_TO_PLOT_RED) {
  output_basename <- paste0(output_basename, "_red")
}

if (CLUSTERS_MERGED) {
  output_basename <- paste0(output_basename, "_merged")
                         
}

output_dir <- paste0("05.0.3-", # Script info
                     format(Sys.time(), "%Y-%m-%d"), # Date info
                     " Plot_XY_By_Cluster")

print("Input parameters loaded, reading needed files")



#### Read needed files #########################################################
layout_filename <- sub("RX", paste0("R", LAYOUT_ROUND), LAYOUT_BASENAME)
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)
clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
if (CLUSTERS_MERGED) {
  clusters_filename <- sub(".csv", "_merged.csv", clusters_filename)
}
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)

# Update CLUSTERS_TO_PLOT if all clusters are to be plotted
if (is.null(CLUSTERS_TO_PLOT)) {
  CLUSTERS_TO_PLOT <- sort(unique(clusters_in[[paste0("R", CLUSTER_ROUND)]]))
  }

print("Prepping data to plot")



#### Prep dataframe for plotting ###############################################
plot_df <- as.data.frame(cbind(layout_in, clusters_in))

if (LAYOUT_ROUND > 1) {
  colnames(plot_df) <- c("group", "umap_x", "umap_y", "cluster")
  plot_df <- plot_df[which(plot_df[,"group"]==LAYOUT_GROUP),]
} else {
  colnames(plot_df) <- c("umap_x", "umap_y", "cluster")
}



#### Randomize order for plotting ##############################################
# to avoid cells at the end (often the last file concatenated) going on top of everything else
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)),]
row.names(plot_df) <- NULL #remove old indices, so everything has a new order

print("Data ready to plot, plotting")



#### Generate output folder (if necessary) and set as working directory ########
if (dir.exists(output_dir)) { # If output directory exists, set output directory as new working directory
  print("output_dir exists")
  setwd(output_dir)
  } else { # Create output directory if it does not exists
    dir.create(output_dir)
    setwd(output_dir)
    }



#### Save plots colored by each marker #########################################
# plot_df_ColorColumn <- plot_df %>% 
#   mutate(grouped_colors = ifelse(test = plot_df[ , "cluster"] %in% CLUSTERS_TO_PLOT, 
#                                  yes = plot_df[["cluster"]], 
#                                  no = "gray"))

# Define the colors to be used
if (ColorPaletteToUse == "viridis") {
  colors.final <- viridis(n = length(unique(plot_df[["cluster"]])),
                          option = ColorToUse)
  
  colors.final[-which(sort(unique(plot_df[["cluster"]])) %in% CLUSTERS_TO_PLOT)] <- OmittedColor
  if (CLUSTERS_TO_PLOT_RED) {
    colors.final[which(sort(unique(plot_df[["cluster"]])) %in% CLUSTERS_TO_PLOT)] <- "red"
  }
  
  } else if (ColorPaletteToUse == "RColorBrewer") {
    colors.final <- colorRampPalette(color = brewer.pal(n = as.numeric(ColorToUse[1]),
                                                        name = ColorToUse[2]))(length(unique(plot_df[["cluster"]])))
  
    colors.final[-which(sort(unique(plot_df[["cluster"]])) %in% CLUSTERS_TO_PLOT)] <- OmittedColor
    if (CLUSTERS_TO_PLOT_RED) {
      colors.final[which(sort(unique(plot_df[["cluster"]])) %in% which(sort(unique(plot_df[["cluster"]])) %in% CLUSTERS_TO_PLOT))] <- "red"
      }
    } else if (ColorPaletteToUse == "Rainbow") {
      numclusters <- length(unique(plot_df[["cluster"]])) # Calculate number of clusters to color
      colors.final <- rainbow(numclusters, alpha = 1) # Generate color palette: any RColorBrewer or viridis palette
      colors.final[-which(sort(unique(plot_df[["cluster"]])) %in% CLUSTERS_TO_PLOT)] <- OmittedColor
      if (CLUSTERS_TO_PLOT_RED) {
        colors.final[which(sort(unique(plot_df[["cluster"]])) %in% CLUSTERS_TO_PLOT)] <- "red"
        }
      }

colors.final.df <- data.frame(cluster = sort(unique(plot_df[["cluster"]])), 
                              grouped_colors = colors.final)
names(colors.final) <- colors.final.df$cluster

plot_df_ColorColumn <- merge(plot_df, colors.final.df, 
                             by = "cluster")

OmittedColor_Rows <- plot_df_ColorColumn[["grouped_colors"]] == OmittedColor

p <- ggplot(data = plot_df_ColorColumn, 
            aes(x = umap_x, y = umap_y, 
                color = factor(cluster))) +
  
  geom_point(data = plot_df_ColorColumn[OmittedColor_Rows, ],
             size = POINT_SIZE, color = OmittedColor, alpha = 0.5) +
  
  geom_point(data = plot_df_ColorColumn[!OmittedColor_Rows, ],
             aes(color = factor(plot_df_ColorColumn[!OmittedColor_Rows, "cluster"])),
             size = POINT_SIZE, alpha = 0.75) +
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  
  scale_color_manual(values = colors.final) +
  
  guides(color = guide_legend(title = "Cluster",
                              override.aes = list(shape = 15, size = 8)))
#^^should work for changing size/shape of legend elements... might have to tweak size per preference

output_filename <- paste0(output_basename, ".", OUTPUT_DEVICE)
ggsave(output_filename, p,  device = OUTPUT_DEVICE, height = 21,width = 21)



#### Reset working directory ###################################################
setwd("..")

print("End UMAP_ColorByMetadata")
print("File outputted")
print("End Plot_XY_By_Cluster")
