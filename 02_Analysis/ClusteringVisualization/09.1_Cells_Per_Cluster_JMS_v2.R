
# Author: Updated by Jonathon Sewell (source: 09.1_Cells_Per_Cluster.R)
# Date: 2026-02-24
# Version: 2
# Update:
# (1) Cleaned up and fixed bugs



print("Start Cells_Per_Cluster.R")

rm(list = ls())
.libPaths( c("/project/zunderlab/R/4.0.3_Pipeline2", .libPaths()) )

library(data.table)
library(dplyr)
library(cluster)
library(ggfortify)
library(ggplot2)
library(gridExtra)


Original_WD <- getwd()


### Input parameters ###########################################################
## User-defined
# Select the cluster round being used
  CLUSTER_ROUND <- 1
# Do you want to incorporate defined groups from cluster_RX_groups.csv?
  USE_GROUPS <- TRUE # TRUE / FALSE
# By which variable do you want to count for the output plot?
  SummarizeBy <- c("Group_Labels") # "Clusters" "Expt", "Treatment", "Group_Assigns", "Group_Labels"
# What is the output plot file type?
  PERCENTAGES_OUTPUT_TYPE <- "pdf" # "pdf" / "png"


## Input file names
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
GROUPS_BASENAME <- "cluster_RX_groups.csv"
CELLS_PER_CLUSTER_OUT <- "cells_per_cluster_RX.csv"

fileNums_filename <- "filenums.csv"
clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
metadata_filename <- "metadata.csv"
groups_filename <- sub("RX", paste0("R", CLUSTER_ROUND), GROUPS_BASENAME)


### Read in necessary files and merge into working dataframe ###################
# Read in input files
fileNums_in <- read.table(fileNums_filename, header=TRUE, sep=",",
                          check.names=FALSE) %>% rename("File_Num" = "V1")

clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE) %>% rename("Clusters" = paste0("R", CLUSTER_ROUND))

metadata_in <- read.table(metadata_filename, header=TRUE, sep=",",
                          check.names=FALSE)

fileNums_metadata <- merge(fileNums_in, 
                           metadata_in[which(!colnames(metadata_in) %in% c("File_Name", "Reordered", "BarcodeSet", "SampleCollection", "Event_Count"))], 
                           by = "File_Num")

# Combine dataframes to be used for analysis
if (USE_GROUPS) {
  groups_in <- read.table(groups_filename, header=TRUE, sep=",",
                          check.names=FALSE) %>% rename("Clusters" = paste0("R", CLUSTER_ROUND))
  
  input_df <- merge(cbind(fileNums_metadata, clusters_in), groups_in[c("Clusters", "Group_Assigns", "Group_Labels")], by = "Clusters")
  
  } else {
    
    input_df <- cbind(fileNums_metadata, clusters_in)
    
    }



### Calculate percentages based on input data and user-defined variables #######
percents_df <- input_df %>% 
  group_by(across(all_of(SummarizeBy[SummarizeBy != "Replicate"]))) %>% 
  summarize(Count = table(.data[[SummarizeBy]])[1]) %>% 
  mutate(Percent = (Count / sum(Count)) * 100)



### Generate and save plot from calculated values ##############################
# Generate plot
PERCENTS_PLOT <- percents_df %>% 
  ggplot(aes(x = factor(as.character(.data[[SummarizeBy]]), levels = sort(unique(.data[[SummarizeBy]]))), 
             y = Percent)) + 
  geom_segment(aes(x = factor(as.character(.data[[SummarizeBy]]), levels = sort(unique(.data[[SummarizeBy]]))), 
                   xend = factor(as.character(.data[[SummarizeBy]]), levels = sort(unique(.data[[SummarizeBy]]))),
                   y = 0, yend = Percent)) +
  
  geom_point(aes(color = factor(as.character(.data[[SummarizeBy]]), levels = sort(unique(.data[[SummarizeBy]])))),
             size = 5) +
  geom_text(aes(label = round(Percent, digits = 2)),
            nudge_y = 0.5,
            nudge_x = 0.25,
            size = 4,
            angle = 30) +
  scale_y_continuous(breaks = seq(0, 
                                  max(percents_df["Percent"], na.rm = TRUE), 
                                  by = round(signif(max(percents_df["Percent"], na.rm = TRUE), digits = 1) / 10))
                     ) +
  labs(title = "Percent by Cluster",
       color = "Cluster") + 
  xlab("Cluster") + 
  annotation_custom(tableGrob(percents_df["Count"],
                              theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.1, fontsize = 10)),
                                                     rowhead = list(fg_params = list(fontface = "bold", fontsize = 10)))), 
                    xmin = nrow(percents_df) * 0.75, xmax = nrow(percents_df) * 1.1, 
                    ymin = 0, ymax = max(percents_df["Percent"], na.rm = TRUE))
  

# Create an output folder for resulting images
output_dir <- paste0("09.1-", # Script ID
                     format(Sys.time(), "%Y-%m-%d"), # Date info
                     "_ClusterPercentages_", # Description
                     "cluster-R", CLUSTER_ROUND, # Cluster round
                     if (USE_GROUPS) { # Clusters merged or not
                       paste0("_Grouped-", paste0(sort(unique(groups_in[["Group_Labels"]])), collapse = ","))
                       }
                     )
dir.create(output_dir)
setwd(paste0(Original_WD, "/", output_dir))


# Save and output plot
ggsave(filename = paste0("ClusterPercentages_",
                          "By", SummarizeBy,
                          ".", PERCENTAGES_OUTPUT_TYPE),
       plot = PERCENTS_PLOT,
       width = 7, height = 5,
       units = "in", dpi = 300
       )
  

### Write csv file with values from the output plot ############################
# Determine the count of each cluster
clusters_out <- table(clusters_in)

clusters_out <- cbind(as.numeric(names(clusters_out)), as.matrix(clusters_out))

colnames(clusters_out) <- c(paste0("Cluster_", "R", CLUSTER_ROUND), "Count")
print(clusters_out)

outname <- sub("RX", paste0("R", CLUSTER_ROUND), CELLS_PER_CLUSTER_OUT)
fwrite(clusters_out, file = outname)



### Reset directory ############################################################
setwd(Original_WD)


print("End Cells_Per_Cluster.R")
print("Complete! Have fun staring at your new plot!")
