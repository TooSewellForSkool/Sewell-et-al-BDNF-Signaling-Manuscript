## Plot clusters by their UMAP layout coordinates

rm(list = ls())
.libPaths( c("/project/zunderlab/R/4.0.3_Pipeline2", .libPaths()) )

library(data.table)

## Input parameters

CLUSTER_ROUND <- 1 # Need to change this for the level of subsetting/subclustering!
CLUSTERS_MERGED <- FALSE # TRUE / FALSE

CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
GROUPS_BASENAME <- "cluster_RX_groups.csv"
LAYOUT_BASENAME <- "umap_xy_RX.csv"
FULL_OUTPUT_BASENAME <- "cluster_RX_xy_All_Full"
IND_OUTPUT_BASENAME <- "cluster_RX_xy_All_Ind"
ALL_GROUP_OUTPUT_BASENAME <- "cluster_RX_xy_GX_All"
IND_GROUP_OUTPUT_BASENAME <- "cluster_RX_xy_GX_Ind"
PLOT_FMT <- "png" #png or pdf ---note that label background is broken in pdf
PLOT_WTH <- 1000 #1000 for png, 10 for pdf
PLOT_HGT <- 1000 #1000 for png, 10 for pdf
IND_COL <- "red"
#POINT_SIZE <- 20

# Try and calculate the best point size --may want to try and find
# a better function that gives good plots for all datasets.
point.size <- function(n) {
  #pt.size <- 0.5+10/(nrow(um_R1)^0.5) # min=0.5, max=100
  #pt.size <- 1 + 20*(0.999^(n*10))
  if (n <= 70) {
    pt.size <- 25
  } else if (n <= 100000000) {
    pt.size <- 20/log10(0.1*n) - 2.5  # meets at 23.66589
  } else {
    pt.size <- 2.5/log10(0.1*n) # meets at 0.3571429
  }
  return(pt.size)
}

# PNG output function - will be used multiple times if R > 1
plot_all_clusters <- function(xy, cl, outfile, plt.fmt, plt.wth, plt.hgt, pt.size) {
  #get cluster centers to use for plotting cluster labels
  clust_centers <- c()
  for (c in sort(unique(cl))) {
    clust_centers <- rbind(clust_centers, colMeans(xy[which(cl==c),]))
  }
  
  # Write out the plot as PNG or PDF
  if (plt.fmt == "png") {
    png(paste0(outfile,".png"), width=plt.wth, height=plt.hgt)
  } else if (plt.fmt == "pdf") {
    pdf(paste0(outfile,".pdf"), width=plt.wth, height=plt.hgt)
  }
  # Plot cells
  plot(x=xy[,1], y=xy[,2], pch=".", col=cl, cex=pt.size,
       xaxt="n", yaxt="n", ann=FALSE)
  # Get center and plot label for each cluster
  for (c in sort(unique(cl))) {
    # Calculate center as mean_X and mean_Y
    clust_center <- colMeans(xy[which(cl==c),])
    # Make "background box" text
    text_bg = paste(rep("\U2588", nchar(c)), collapse="")
    # Plot semi-transparent background box and cluster labels
    text(clust_center[1], clust_center[2], text_bg, cex=2,
         col=rgb(red=1, green=1, blue=1, alpha=0.75))
    text(clust_center[1], clust_center[2], c, cex=2)
  }
  dev.off()
}

# PNG output function - will be used multiple times if R > 1
plot_ind_clusters <- function(xy, cl, outfile, plt.fmt,
                              plt.wth, plt.hgt, pt.size,
                              cl.col="red") {
  
  unique_cl <- sort(unique(cl))
  
  pr <- ceiling(sqrt(length(unique_cl)))
  
  adj.wth <- plt.wth*pr/2
  adj.hgt <- plt.hgt*pr/2
  
  # Write out the plot as PNG or PDF
  if (plt.fmt == "png") {
    png(paste0(outfile,".png"), width=adj.wth, height=adj.hgt)
  } else if (plt.fmt == "pdf") {
    pdf(paste0(outfile,".pdf"), width=adj.wth, height=adj.hgt)
  }
  par(mfrow=c(pr,pr))
  for (c in unique_cl) {
    tmp_cl <- cl
    tmp_cl[which(cl!=c)] <- 0
    tmp_cl[which(cl==c)] <- 1
    if (any(cl!=c)) {
      plot(x=xy[which(cl!=c),1], y=xy[which(cl!=c),2], pch=".", cex=pt.size, col="grey",
           main=paste("Cluster",c), axes=FALSE, ylab="", xlab="", cex.main=5,
           xlim=c(min(xy[,1]), max(xy[,1])), ylim=c(min(xy[,2]), max(xy[,2])))
      points(x=xy[which(cl==c),1], y=xy[which(cl==c),2], pch=".", cex=pt.size*2, col=cl.col)
    } else {
      plot(x=xy[which(cl==c),1], y=xy[which(cl==c),2], pch=".", cex=pt.size*2, col=cl.col,
           main=paste("Cluster",c), axes=FALSE, ylab="", xlab="", cex.main=5)
    }
    
    box()
    
  }
  dev.off()
}


# Read in cluster_assigns file
clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
if (CLUSTERS_MERGED) {
  clusters_filename <- sub(".csv", "_merged.csv", clusters_filename)
}

clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)
# take only clusters from the round of interest
cl <- clusters_in[,paste0("R", CLUSTER_ROUND)]

# input umap xy coordinates from full dataset (always R1)
layout_filename <- sub("RX", "R1", LAYOUT_BASENAME)
layout_R1 <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)

point_size <- point.size(nrow(layout_R1))

# First, plot the current clusters on the full dataset
full_outname <- sub("RX", paste0("R", CLUSTER_ROUND), FULL_OUTPUT_BASENAME)
if (CLUSTERS_MERGED) {
  full_outname <- paste0(full_outname, "_merged")
}

plot_all_clusters(layout_R1, cl, full_outname, PLOT_FMT, PLOT_WTH,
                  PLOT_HGT, point_size)

ind_outname <- sub("RX", paste0("R", CLUSTER_ROUND), IND_OUTPUT_BASENAME)
if (CLUSTERS_MERGED) {
  ind_outname <- paste0(ind_outname, "_merged")
}
plot_ind_clusters(layout_R1, cl, ind_outname, PLOT_FMT, PLOT_WTH,
                  PLOT_HGT, point_size)

# Input umap xy coordinates from subsets and plots these
# individually if R > 1
if (CLUSTER_ROUND > 1) {
  # Read in current subset cluster grouping
  layout_filename <- sub("RX", paste0("R", CLUSTER_ROUND), LAYOUT_BASENAME)
  layout_RX <- read.table(layout_filename, header=TRUE, sep=",",
                          check.names=FALSE)
  # Get all the unique groups
  groups <- sort(unique(layout_RX[,"Group"]))
  nonzero_groups = groups[groups > 0]
  # For each group, pull out clusters and XY coordinates, then plot
  for (g in nonzero_groups) {
    # Get basename for output file
    all_basename <- sub("RX", paste0("R", CLUSTER_ROUND), ALL_GROUP_OUTPUT_BASENAME)
    all_basename <- sub("GX", paste0("G", g), all_basename)
    if (CLUSTERS_MERGED) {
      all_basename <- paste0(all_basename, "_merged")
    }
    # Get indices for all group members
    sub_idxs <- which(layout_RX[,"Group"] == g)
    # Calculate what will hopefully be a point size for this plot
    point_size <- point.size(length(sub_idxs))
    # Plot by UMAP_XY for all group members, colored by cluster assigns
    plot_all_clusters(layout_RX[sub_idxs,c(paste0("R", CLUSTER_ROUND, "_X"),
                                           paste0("R", CLUSTER_ROUND, "_Y"))],
                      cl[sub_idxs], all_basename, PLOT_FMT, PLOT_WTH,
                      PLOT_HGT, point_size)
    
    ind_basename <- sub("RX", paste0("R", CLUSTER_ROUND), IND_GROUP_OUTPUT_BASENAME)
    ind_basename <- sub("GX", paste0("G", g), all_basename)
    plot_ind_clusters(layout_RX[sub_idxs,c(paste0("R", CLUSTER_ROUND, "_X"),
                                           paste0("R", CLUSTER_ROUND, "_Y"))],
                      cl[sub_idxs], ind_basename, PLOT_FMT, PLOT_WTH,
                      PLOT_HGT, point_size, IND_COL)
  }
}

