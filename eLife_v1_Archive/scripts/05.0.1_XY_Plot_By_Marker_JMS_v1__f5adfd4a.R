
# Author: Jonathon Sewell, adapted from 06_Plot_Expression_Heatmap.R (author Corey Williams - 15 July, 2019)
# Date: 2025-08-19
# Version: 1
# Update:
# (1) Input expression matrices using fread from data.table package for quicker import
# (2) Updated aes_string() in ggplot to aes() - aes_string() "was deprecated in ggplot2 3.0.0."
# (3) Updated output folder name (FOLDER_BASENAME)


#Corey Williams, University of Virginia
#15 Jul, 2019
#Plot colored by expression of markers

print("Start Plot_XY_By_Marker.R")

rm(list = ls())
.libPaths( c( .libPaths(), "/project/zunderlab/R/4.0.3_Pipeline2") )

library(ggfortify)
library(data.table)

## Input parameters
LAYOUT_ROUND <- 2
LAYOUT_GROUP <- 1 # Don't use group 0, not assigned XY coordinates
FOLDER_BASENAME <- paste0("05.0.1-", # Script info
                          format(Sys.time(), "%Y-%m-%d"), # Date info 
                          " Marker_Expression_VAR_XY")
EXPRS_MAT_INFILE <- "expression_matrix_analysis.csv"
EXPRS_MAT_OTHER <- "expression_matrix_other.csv"
LAYOUT_BASENAME <- "umap_xy_RX.csv"
PANEL_FILENAME <- "panel.csv"
POINT_SIZE <- 0.1
PLOT_HEIGHT <- 7
PLOT_WIDTH <- 7

if (LAYOUT_ROUND > 1) {
    output_dir <- sub("VAR", paste0("L", LAYOUT_ROUND,
                                    "_G", LAYOUT_GROUP),
                      FOLDER_BASENAME)
    layout_cols <- c("group", "umap_x", "umap_y")
} else {
    output_dir <- sub("VAR", paste0("L", LAYOUT_ROUND),
                      FOLDER_BASENAME)
    layout_cols <- c("umap_x", "umap_y")
}

## Read needed files
exprs_mat <- cbind(fread(EXPRS_MAT_INFILE, 
                         stringsAsFactors = FALSE, header=TRUE, data.table = FALSE),
                   fread(EXPRS_MAT_OTHER, 
                         stringsAsFactors = FALSE, header=TRUE, data.table = FALSE))

fread(EXPRS_MAT_INFILE, 
      stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)

layout_filename <- sub("RX", paste0("R", LAYOUT_ROUND), LAYOUT_BASENAME)
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)
colnames(layout_in) <- layout_cols

panel <- read.table(PANEL_FILENAME, header=TRUE, sep=",", check.names=FALSE)

## Prep dataframe for plotting
plot_vars <- panel[which(panel[,"Plot"]==1), "Fixed_Param"]
plot_df <- as.data.frame(cbind(exprs_mat[,plot_vars], layout_in))

if (LAYOUT_ROUND > 1) {
    plot_df <- plot_df[which(plot_df[,"group"]==LAYOUT_GROUP),]
}

## Randomize order for plotting, to avoid cells at the end (often the last file concatenated) going on top of everything else
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)),]
row.names(plot_df) <- NULL #remove old indices, so everything has a new order

## Save plots colored by each marker

dir.create(output_dir)
setwd(output_dir)
#loop through variables to plot
for (v in plot_vars){
  ggsave(paste0(v,".png"),
         plot=ggplot(plot_df, aes(x=.data[["umap_x"]], y=.data[["umap_y"]], color=.data[[v]])) +
           geom_point(size = POINT_SIZE) +
           scale_color_viridis_c() + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")),
         height = PLOT_HEIGHT, width = PLOT_WIDTH)
}
setwd("..")

print("End Plot_XY_By_Marker.R")
