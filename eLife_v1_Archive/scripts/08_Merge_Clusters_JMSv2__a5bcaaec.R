
# Updated script:
# Jonathon Sewell
# Version: 2
# 2026-02-24

# Source script: 08_Merge_Clusters.R


## Merge Clusters


rm(list = ls())
.libPaths( c("/project/zunderlab/R/4.0.3_Pipeline2", .libPaths()) )

library(data.table)


INPUT.FOLDER <- getwd() # setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) to set source file location as working directory


### Input parameters ###########################################################
## User-defined
  CLUSTER_ROUND <- 1 # Need to change this for the level of subsetting/subclustering!

  ScriptsToTransfer <- c("05.1_3D_UMAP_v4.R",
                         "12.1_PlotSignalingHeatmap_ByCluster_BySample_PercAboveThreshold_NormToCtrl_v6.R",
                         "16_Barplot_v7.r",
                         "LineGraphAnalysis_v10.R")

## Naming variables
  CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
  clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
  MERGE_KEY <- "cluster_RX_merge.csv"
  merge_key_filename <- sub("RX", paste0("R", CLUSTER_ROUND), MERGE_KEY)
    # Note: 
    # - User must make this MERGE_KEY csv file separately
    # - Row numbers become the new merged clusters, which are filled with old clusters to be merged
  
  OUTPUT.FOLDER.NAME <- paste0("08-", # Script ID
                               format(Sys.time(), "%Y-%m-%d"), # Date info
                               "_MergeClusterAnalysis")

  
  
## Filenames for files that will be transferred to output folder
  LAYOUT_BASENAME <- "/umap_xy_RX.csv"
  layout_filename <- list.files(INPUT.FOLDER, pattern = ".csv") %>% # finds all files in Original_WD of the input_file_type
    (\(x) x[grepl("umap_xy", x)])() # pass the list of file names to an anonymous function that filters for input_file_transformation
  PANEL_FILENAME <- "/panel.csv"
  METADATA_FILENAME <- "/metadata.csv"
  FILENUMS_FILENAME <- "/filenums.csv"
  EXPRS_ANALYSIS_FILENAME <- "/expression_matrix_analysis.csv"
  EXPRS_OTHER_FILENAME <- "/expression_matrix_other.csv"
  
  
  
### Read in necessary files ####################################################
## Read in files necessary for merging
  clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                            check.names=FALSE)

  merge_key <- read.table(merge_key_filename, header=FALSE, sep=",",
                          check.names=FALSE)

print("Necessary files read in")




### Generate the merged cluster file ###########################################
# Take only clusters from the round of interest
cl_in <- clusters_in[,paste0("R", CLUSTER_ROUND)]

cl_out <- cl_in

for (i in 1:nrow(merge_key)) {
  to_merge <- merge_key[i,][!is.na(merge_key[i,])] #remove NAs from uneven merge nums
  merge_min <- min(to_merge)
  cl_out[which(cl_out %in% to_merge)] <- merge_min
}

tmp_clust_names <- sort(unique(cl_out))

remove_gaps_key <- cbind(tmp_clust_names, 1:length(tmp_clust_names))

for (i in 1:nrow(remove_gaps_key)) {
  if (remove_gaps_key[i,1] != remove_gaps_key[i,2]) {
    cl_out[which(cl_out == remove_gaps_key[i,1])] = remove_gaps_key[i,2]
  }
}

clusters_outfile <- sub(".csv", "_merged.csv", clusters_filename)

df_out <- data.frame(cl_out)
colnames(df_out) <- paste0("R", CLUSTER_ROUND)



### Write output files #########################################################
## Create and set ouput directory
dir.create(OUTPUT.FOLDER.NAME) #Create output directory
OUTPUT.FOLDER <- paste0(INPUT.FOLDER, "/", OUTPUT.FOLDER.NAME)
setwd(OUTPUT.FOLDER) #Set directory to output folder



### Save files into OUTPUT.FOLDER ##############################################
## Output new files for downstream analysis
fwrite(df_out, 
       file = clusters_outfile, 
       row.names = FALSE)


## Output original files into OUTPUT.FOLDER so they are available for downstream analysis
fwrite(merge_key, 
       file = merge_key_filename, 
       row.names = FALSE)

fwrite(clusters_in,
       file = sub("/", "", sub("/", "", clusters_filename)), # Remove slash at beginning of clusters_filename
       row.names = FALSE)

file.copy(from = paste0(INPUT.FOLDER, EXPRS_ANALYSIS_FILENAME), 
          to = paste0(OUTPUT.FOLDER, EXPRS_ANALYSIS_FILENAME))

file.copy(from = paste0(INPUT.FOLDER, EXPRS_OTHER_FILENAME), 
          to = paste0(OUTPUT.FOLDER, EXPRS_OTHER_FILENAME))

file.copy(from = paste0(INPUT.FOLDER, FILENUMS_FILENAME), 
          to = paste0(OUTPUT.FOLDER, FILENUMS_FILENAME))

file.copy(from = paste0(INPUT.FOLDER, METADATA_FILENAME), 
          to = paste0(OUTPUT.FOLDER, METADATA_FILENAME))

file.copy(from = paste0(INPUT.FOLDER, PANEL_FILENAME), 
          to = paste0(OUTPUT.FOLDER, PANEL_FILENAME))

sapply(layout_filename, function(NAME) {
  file.copy(from = paste0(INPUT.FOLDER, "/", NAME), 
            to = paste0(OUTPUT.FOLDER, "/", NAME))
})

sapply(ScriptsToTransfer, function(NAME) {
  file.copy(from = paste0(INPUT.FOLDER, "/", NAME), 
            to = paste0(OUTPUT.FOLDER, "/", NAME))
})


#### Return to original working directory ######################################
setwd(INPUT.FOLDER)



