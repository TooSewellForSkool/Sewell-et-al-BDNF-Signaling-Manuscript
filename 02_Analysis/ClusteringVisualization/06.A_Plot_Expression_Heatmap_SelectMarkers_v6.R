
# Author: Jonathon Sewell, adapted from 06_Plot_Expression_Heatmap.R
# Date: 2025-08-24
# Version: 6
# Update:
# (1) Enabled user-defined selection of markers
# (2) Enabled choice of mean or median for plotted values


# Install ComplexHeatmap with BiocManager
# https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
# Instructions/guide for ComplexHeatmap:
# https://jokergoo.github.io/ComplexHeatmap-reference/book/

rm(list = ls())
.libPaths( c("/project/zunderlab/R/4.0.3_Pipeline2", .libPaths()) )

library(ComplexHeatmap)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
library(dendsort)
library(dendextend) #not used here, but can use to reorder/flip dendrogram branches
library(viridis)
library(stringr)
library(data.table)


Original_WD <- getwd()


#### Variables to Change ######################################################
# Define cluster variables
CLUSTER_ROUND <- 1 # Need to change this for the level of subsetting/subclustering!
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
CLUSTERS_MERGED <- FALSE # TRUE / FALSE

ClustersToPlot <- c(1, 7, 2, 13, 4) # c(1, 7, 2, 13, 4) / "All"

ClusterOrder <- c()
  # c(2, 11, 4, 1, 14, 15, 5, # Gradient A
  #                 12, 10, 6, # Gradient B
  #                 19, 21, # Astrocytic
  #                 8, 18, # Microglial
  #                 22, # Oligodendrocytic
  #                 16, 20, # Other Progenitor-like
  #                 17, # Misc. Neuron
  #                 7, 13, # Misc. Glia
  #                 3, 9) # Unknown
    # c() if you want to allow heatmap to group clusters independent of user

# Define marker variables
UserSelectsMarkers <- TRUE

if (UserSelectsMarkers) {
  UserDefinedMarkers <- c("Tuj1_89Y", "NeuN_150Nd", "Sox2_160Gd", "DCX_163Dy", "MAP2_164Dy", "BLBP_172Yb")
  # UserDefinedMarkers <- c(# "pSrc_154Sm", "p-cJUN_162Dy", "pPLCG-2_144Nd",
  #                         "pPLK1_113In", "H3K9ac_115In", "pGSK3b_139La", 
  #                         # "p-p38_142Nd", "p-p53_151Eu", "pSTAT3_153Eu", 
  #                         "pATM_159Tb", "yH2AX_165Ho","pAKt_168Er", "pERK_171Yb", "pS6_175Lu", "pCreb_176Lu")
  # "Tuj1_89Y", "pPLK1_113In", "H3K9ac_115In", "pGSK3b_139La", "CXC3CR1_141Pr", "p-p38_142Nd", "Vimentin_143Nd", "pPLCG-2_144Nd", "OligoO4_145Nd",  
  # "s100b_146Nd", "pSTAT5_147Sm", "HB9_148Nd", "CGRP_149Sm", "NeuN_150Nd", "p-p53_151Eu", "CHAT_152Sm", "pSTAT3_153Eu", "pSrc_154Sm", "nNOS_155Gd"     
  # "TrkB_s_158Gd", "pATM_159Tb", "Sox2_160Gd", "p75NTR_161Dy", "p-cJUN_162Dy", "DCX_163Dy", "MAP2_164Dy", "yH2AX_165Ho", "pNFKB_166Er", "Islet1_167Er", 
  # "pAKt_168Er", "GFAP_169Tm", "pTrkB_170Er", "pERK_171Yb", "BLBP_172Yb", "CC3_173Yb", "TrkB_i_174Yb", "pS6_175Lu", "pCreb_176Lu", "DNA_191Ir", "DNA_193Ir"
  }

# Define output selections and naming
MeanOrMedian <- "median" # "mean" / "median"
HeatmapCellLabels <- TRUE # TRUE / FALSE
OUTPUT.TYPE <- "pdf" # choose png or pdf

OUTPUT_BASENAME <- paste0("cluster_RX_exprs_heatmap_SelectMarkers_", 
                          if (length(ClustersToPlot) <= 1) {
                                "All"
                            } else {
                              paste0("Clusters", min(ClustersToPlot), "-", max(ClustersToPlot))
                              },
                          if (HeatmapCellLabels) {
                            paste0("-labeled")
                          },
                          if (length(ClusterOrder) > 1) {
                            paste0("-ordered")
                          }
                          )

# Export csv file of data being plotted in heatmaps
WRITE_MEANS <- TRUE # TRUE / FALSE


#### Variables that may need to be changed #####################################
FILE_NUM <- paste0("filenums", ".csv")
EXPRS_MAT_INFILE1 <- "expression_matrix_analysis.csv"
EXPRS_MAT_INFILE2 <- "expression_matrix_other.csv"
PANEL_INFILE <- paste0("panel", ".csv")
MEANS_BASENAME <- paste0("cluster_RX_exprs_means", ".csv")
METADATA_FILENAME <- paste0("metadata", ".csv")
clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
if (CLUSTERS_MERGED) {
  clusters_filename <- sub(".csv", "_merged.csv", clusters_filename)
  OUTPUT_BASENAME <- paste0(OUTPUT_BASENAME, "_merged")
  MEANS_BASENAME <- paste0(MEANS_BASENAME, "_merged")
}


Original_WD <- getwd()



#### Read in necessary files ###################################################
# Read in file numbers associated with the expression file
FileNum <- read.table(FILE_NUM, header=TRUE, sep=",",
                      check.names=FALSE)

# Read in panel (will show which expression markers will be included)
panel <- read.table(PANEL_INFILE, header=TRUE, sep=",",
                    check.names=FALSE)

# Read in metadata
metadata <- read.table(METADATA_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)

# Read in cluster file
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)
# Take only clusters from the round of interest
clusts <- clusters_in[,paste0("R", CLUSTER_ROUND)]

# Read in expression matrix and perform initial clean up of data
exprs_mat_IN1 <- fread(paste0(Original_WD, "/", EXPRS_MAT_INFILE1), 
                       stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)
exprs_mat_IN2 <- fread(paste0(Original_WD, "/", EXPRS_MAT_INFILE2), 
                       stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)
exprs_mat_IN_ALL <- cbind(exprs_mat_IN1, exprs_mat_IN2)

exprs_mat_IN <- cbind(FileNum, #Append FileNum for indexing
                      clusters_in, #Append clusters_in for indexing
                      exprs_mat_IN_ALL[ , panel$Fixed_Param[which(panel$Plot==1)]]) #Only selects markers from panel$Plot==1
colnames(exprs_mat_IN) <- c("File_Num", "Cluster", colnames(exprs_mat_IN[3:length(colnames(exprs_mat_IN))]))

# Adjust ClustersToPlot, if necessary
if (length(ClustersToPlot) <= 1) {
  ClustersToPlot <- sort(unique(clusts))
  } else {
    ClustersToPlot <- ClustersToPlot
    }

exprs_mat_IN <- exprs_mat_IN[exprs_mat_IN$Cluster %in% ClustersToPlot, ]



#### Create variables to be used later #########################################
# Markers used to Cluster = "Markers_Analysis"
  Markers_Analysis <- panel$Fixed_Param[panel$Analysis == 1]
# All markers in panel = "Markers_Plot"
  Markers_Plot <- panel$Fixed_Param[panel$Plot == 1] 
# ID Only = "ID_Markers"
  ID_Markers  <-  c(Markers_Plot[!stringr::str_starts(Markers_Plot, "p") & # Remove all markers starting with "p" (most signaling)
                                   !Markers_Plot %in% c("H3K9ac_115In", "yH2AX_165Ho", "CC3_173Yb", # Remove all other non-ID markers
                                                        "DNA_191Ir", "DNA_193Ir")],
                    "p75NTR_161Dy") # Add back ID markers starting with "p"
# Signaling markers = "Signaling_Markers"
  Signaling_Markers <- Markers_Plot[!Markers_Plot %in% ID_Markers &
                                      !Markers_Plot == "DNA_191Ir" &
                                      !Markers_Plot == "DNA_193Ir"]

# Combine marker groupings into a list
  MarkerSets <- list(All = Markers_Plot, 
                     IDmarkers = ID_Markers, 
                     Clustering = Markers_Analysis, 
                     Signaling = Signaling_Markers)

  # Append UserDefinedMarkers if applicable 
  if (UserSelectsMarkers) {
    MarkerSets <- append(MarkerSets, list(SelectMarkers = UserDefinedMarkers))
  }

# Create an output folder for resulting images
output_dir <- paste0("06.A-", # Script ID
                     format(Sys.time(), "%Y-%m-%d"), # Date info
                     "_ExprHeatmaps_", # Description
                     "cluster-R", CLUSTER_ROUND, # Cluster round
                     "_", MeanOrMedian, # Values displayed
                     if (CLUSTERS_MERGED) { # Clusters merged or not
                       "_merged"
                       }
                     )
dir.create(output_dir)
setwd(paste0(Original_WD, "/", output_dir))


#### Generate expression matrix and use to generate the heatmaps ###############
# Create a function to generate heatmaps of the selected marker set
ExpressionHeatmap <- function(MarkerIDSet) {

  exprs_mat_ID <- exprs_mat_IN[MarkerSets[[MarkerIDSet]]]
  
  clust_SummaryValue <- c()
  
  apply(exprs_mat_ID[which(exprs_mat_IN$Cluster == 1), ], 2, median, na.rm = TRUE)
  
  unique_clusts <- sort(unique(exprs_mat_IN$Cluster))
  for (i in unique_clusts) {
    if (MeanOrMedian == "mean") {
      clust_SummaryValue <- rbind(clust_SummaryValue,
                                  colMeans(exprs_mat_ID[which(exprs_mat_IN$Cluster == i), ]))
      } else if (MeanOrMedian == "median") {
        clust_SummaryValue <- rbind(clust_SummaryValue,
                                    apply(X = exprs_mat_ID[which(exprs_mat_IN$Cluster == i), ], 
                                          MARGIN = 2, FUN = median, na.rm = TRUE))
        }
    }
  rownames(clust_SummaryValue) <- unique_clusts
  
  # Modify the output filename based on the marker set selected
  if (MarkerIDSet == 1) {
    OUTPUT_BASENAME <- paste0(gsub("Select", "All", OUTPUT_BASENAME), "_", MeanOrMedian)
    } else if (MarkerIDSet == 2) {
      OUTPUT_BASENAME <- paste0(gsub("Select", "ID", OUTPUT_BASENAME), "_", MeanOrMedian)
      } else if (MarkerIDSet == 3) {
        OUTPUT_BASENAME <- paste0(gsub("Select", "Cluster", OUTPUT_BASENAME), "_", MeanOrMedian)
        } else if (MarkerIDSet == 4) {
          OUTPUT_BASENAME <- paste0(gsub("Select", "Signaling", OUTPUT_BASENAME), "_", MeanOrMedian)
          } else if (MarkerIDSet == 5) {
            OUTPUT_BASENAME <- paste0(gsub("Select", paste0("Selected", length(UserDefinedMarkers)), OUTPUT_BASENAME), "_", MeanOrMedian)
            }
  
  # Generate a csv file, if desired
  if (WRITE_MEANS) {
    write.csv(clust_SummaryValue, paste0(sub("RX", paste0("R", CLUSTER_ROUND),
                                    OUTPUT_BASENAME), ".csv"))
    }
  
  # Set up and apply a function to normalize on 0-1 scale
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  norm_means <- apply(clust_SummaryValue,2,range01)
  rownames(norm_means) <- unique_clusts
  
  # Generate dendrogram for columns and rows
  col_dend <- dendsort(hclust(dist(t(norm_means))))
  
  if (length(ClusterOrder) <= 1) {
    row_dend <- dendsort(hclust(dist(norm_means)))
    clust_order <- NULL
  } else { 
    row_dend <- FALSE
    clust_order <- ClusterOrder
      }
  
  # Write label for heatmap column
  hm_cols <- sub("_\\d*[A-Za-z]*$","",colnames(norm_means))
  
  # Define the output file type
  if (OUTPUT.TYPE == "png") {
    png(paste0(sub("RX", paste0("R", CLUSTER_ROUND), OUTPUT_BASENAME), ".png"),
        width=1600, height=800)
    } else if (OUTPUT.TYPE == "pdf") {
      pdf(paste0(sub("RX", paste0("R", CLUSTER_ROUND), OUTPUT_BASENAME), ".pdf"),
          width=16, height=8)
      }

  # Generate the heatmap
  draw(Heatmap(norm_means, 
               col = viridis(100, option="viridis"),
               cluster_rows = row_dend,
               row_dend_width = unit(50, "mm"),
               cluster_columns = col_dend,
               column_dend_height = unit(50, "mm"),
               column_labels = hm_cols,
               row_labels = unique_clusts,
               row_order = clust_order,
               
               # function to add values to tiles
               cell_fun = if (HeatmapCellLabels) {
                 function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.1f", norm_means[i, j]), x, y, gp = gpar(fontsize = 8))
                   }
                 }
               )
       )
  dev.off()
  }

for (i in 1:length(MarkerSets)) {
  ExpressionHeatmap(i)
}



### Reset directory
setwd(Original_WD)

print("Complete! Have fun staring at your new plot!")


