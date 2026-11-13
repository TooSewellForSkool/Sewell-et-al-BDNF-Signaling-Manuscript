#Corey Williams, University of Virginia
#15 Jul, 2019
#Plot colored by expression of markers


print("Start Violin_Plots.R")

rm(list = ls())
.libPaths( c("/project/zunderlab/R/4.0.3_Pipeline2", .libPaths()) )

library(ggfortify)
library(ggstance)
library(ggpubr)
library(forcats)

## Input parameters
CLUSTER_ROUND <- 1
EXPRS_MAT_INFILE <- "expression_matrix_analysis.csv"
EXPRS_MAT_OTHER <- "expression_matrix_other.csv"
PANEL_FILENAME <- "panel.csv"
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
VIOLIN_HEIGHT_FACTOR <- 5

## Read needed files
exprs_mat <- cbind(read.table(EXPRS_MAT_INFILE, header=TRUE, sep=",",
                              check.names=FALSE),
                   read.table(EXPRS_MAT_OTHER, header=TRUE, sep=",",
                              check.names=FALSE))
panel <- read.table(PANEL_FILENAME, header=TRUE, sep=",", check.names=FALSE)

clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)

## Prep dataframe for plotting
plot_vars <- panel[which(panel[,"Plot"]==1), "Fixed_Param"]
plot_df <- as.data.frame(cbind(exprs_mat[,plot_vars],clusters_in))
colnames(plot_df)[ncol(plot_df)] <- "Cluster"

#Make list of violin plots by cluster
plist = sapply(plot_vars, function(marker_plot) {
  if (marker_plot == plot_vars[1]){
    ggplot(plot_df, 
           aes(x = plot_df[,marker_plot], y = fct_rev(factor(Cluster)), 
               fill = factor(Cluster))) + 
      geom_violinh(trim=FALSE,scale = "width") + 
      xlab(marker_plot) + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            legend.position = "none",axis.text.x = element_blank(),axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8))
  }
  else {
    ggplot(plot_df, 
           aes(x = plot_df[,marker_plot], y = fct_rev(factor(Cluster)), 
               fill = factor(Cluster))) + 
      geom_violinh(trim=FALSE,scale = "width") + 
      xlab(marker_plot) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            legend.position = "none",axis.text = element_blank(),axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8))
  }
}, simplify=FALSE)
#save Violin plots
ggsave("ViolinPlots.png",
       annotate_figure(ggarrange(plotlist = plist,ncol=length(plist)),
                       left = text_grob("Cluster",rot=90)),
       height=max(clusters_in)/VIOLIN_HEIGHT_FACTOR,width=length(plot_vars))

print("End 12_ViolinPlots.R")
