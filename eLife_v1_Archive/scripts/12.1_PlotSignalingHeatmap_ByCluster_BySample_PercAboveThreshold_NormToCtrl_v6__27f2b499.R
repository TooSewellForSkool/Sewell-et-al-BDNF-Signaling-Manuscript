
# Author: Jonathon Sewell
# Date: 2025-09-19
# Version: 7
# Update:
# Update:
# (1) Enables selection of output file type
# (2) Removed unnecessary hashed-out code
# (3) Updated formatting

# PURPOSE: Generate a Heatmap of the nth percent of all samples above a user-defined threshold by cluster based on a user-defined parameter
# INPUT:  _ "filenums.csv"    _ "panel.csv"   _ "cluster_RX_exprs_means.csv" (X = clustering round)    
#         _ "expression_matrix_analysis.csv"  _ "expression_matrix_other.csv"
#         _ "metadata.csv" 
#     ***NOTE: make sure "metadata.csv" includes the following columns: 
#            _ File_Name   _ File_Num    _ Reordered       _ BarcodeSet  _ SampleCollection	
#            _ Expt	      _ Condition	  _ Condition_Only  _ TimePoint   _ Replicate

# OUTPUT: _ Complete dataset thresholded by the user-defined parameters ("[...]_RAWvalues.csv")
#         _ Fold Change data from the "RAWvalues" dataset relative to a user-defined condition(s) ("_FC.csv", _FC_RepsAvg.csv)
#         _ log2(Fold change) data of the "FC" dataset ("_log2FC.csv", _log2FC_RepsAvg.csv)
#         _ Heatmap of Fold Change data organized by cluster or by treatment condition (Heatmap_ByCluster_FC_RepsAvg.png, Heatmap_ByCondition_FC_RepsAvg.png)
#         _ Heatmap of log2(FC) data organized by cluster or by treatment condition (Heatmap_ByCluster_log2FC_RepsAvg.png, Heatmap_ByCondition_log2FC_RepsAvg.png)


# NOTES:  _ This script was originally modified from 06.3_Plot_Expression_HeatmapByCluster_BySample.R from Pipeline2
#         _ To normalize data to a specified condition, the exprs_matrix data undergoes a Linear Transform from the input Arcsinh Transform data


rm(list = ls())
.libPaths( c("/project/zunderlab/R/4.0.3_Pipeline2", .libPaths()) )

### Packages to call (or install) ###
# Install ComplexHeatmap with BiocManager
# https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
# install.packages('ComplexHeatmap')
# Instructions/guide for ComplexHeatmap:
# https://jokergoo.github.io/ComplexHeatmap-reference/book/
library(ComplexHeatmap)
library(dendsort)
library(dendextend) #not used here, but can use to reorder/flip dendrogram branches
library(viridis)
# install.packages("stringr")
library("stringr")
library(circlize)
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(gridExtra)
library(grid)
library(ggplot2)
library(cowplot)

Original_WD <- getwd() # setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) to set source file location as working directory

#### Variables to Change ######################################################
# Define how to initially subset the input dataframe based on larger-scale metadata
WorkingSampleSet_Identifier <- "Expt" # i.e. BarcodeSet, SampleCollection, Expt
WorkingSampleSet <- c(5) # Select a sample set from the complete data set to work with (i.e. metadata$SampleCollection == 1)
# Cheatsheet: NoTreat = "1" / ContDepr, Rescue, BNDF = "2, 3, 4" / K252a expt = "5" / Stock = "6" / Universal = "7"

# Define how to set the threshold and on which timepoint to base this threshold
Set_Threshold <- 0.95 # 1st Quartile = 0.25 / 2nd Quartile (Median) = 0.5 / 3rd Quartile = 0.75
Threshold_TimePoint <- 0 # Must match value from metadata$TimePoint

# Define how to subset further based on condition (Note: only need to include base condition, no timepoint information)
ConditionsToExclude <- c() # c() if none are to be excluded
  # Expt 1: "NoTreat"
  # Expt 2: "ContDepr"
  # Expt 3: "CompMed"
  # Expt 4: "BDNF"
  # Expt 5: "NoTreatTest" / BDNFTest" / "K252aTest" / "ContDeprTest"

# Define how to threshold the data (against time=0hr, or against another condition)
CtrlPerTimepoint <- FALSE
  # TRUE if comparison is to a control of the same time point (ex. Expt at 1hr compared to Ctrl at 1hr); 
  # FALSE if comparison is to a control at a different time point (ex. Expt at 1hr compared to Ctrl at 0hr)

# Output what condition(s) to compare to for fold change calculations
PostThreshold_NormCondition <- if (isTRUE(CtrlPerTimepoint)) {
  if (mean(WorkingSampleSet %in% c(1, 2, 3, 4)) > 0) {
    PostThreshold_NormCondition <- c("NoTreat_0hr", 
                                     "ContDepr_5min", "ContDepr_15min", "ContDepr_1hr", "ContDepr_4hr", "ContDepr_16hr", 
                                     # "ContDepr_72hr" # de-hash if using 72hr timepoint
    )
  } else if (mean(WorkingSampleSet %in% c(5)) == 1) {
    PostThreshold_NormCondition <- c("ContDepr_Test")
  } else {
    print("Make sure to include Expt 1")
    PostThreshold_NormCondition <- c("NoTreat_0hr")
  }
} else if (isFALSE(CtrlPerTimepoint)) {
  if (mean(WorkingSampleSet %in% c(1, 2, 3, 4)) > 0) {
    PostThreshold_NormCondition <- c("NoTreat_0hr")
  } else if (mean(WorkingSampleSet %in% c(5)) == 1) {
    PostThreshold_NormCondition <- c("NoTreat_0hr") # Make sure this matches metadata.csv
  } else {
    print("Make sure to include Expt 1")
    PostThreshold_NormCondition <- c("NoTreat_0hr")
  }
} else {
  print("Check CtrlPerTimepoint variable is TRUE/FALSE")
}

# Define whether THRESHOLD_CATEGORY should be considered (only if 13.1_ExtractByMarkerExpression_HiLo_v3.R already used)
IsThereThresholdCategory <- FALSE # TRUE / FALSE - TRUE if using output from 13.1_ExtractByMarkerExpression_HiLo_v3.R
WhichToPlot <- "lo" # hi / med / lo
if (IsThereThresholdCategory) { # Original_WD should refer to directory generated by 13.1_ExtractByMarkerExpression_HiLo_v3.R
  MarkerThresholdedOn <- str_extract(basename(Original_WD), 
                                     "(?<=ExtractByMarkerExpression_)[^_]+_[^_]+")
}

# Define Cluster information
CLUSTER_ROUND <- 1 # Need to change this for the level of subsetting/subclustering!
CLUSTERS_MERGED <- FALSE
ClustersToRemove <- c() # c() if not removing any
ClusterOrder <- c() # c() if ClusterOrder is to go un-defined

# Define if there are any File_Name values to exclude (outliers or otherwise)
Exclude_FileNames <- FALSE # TRUE if there are any samples that are considered to be outliers you'd like to remove
FileNamesToExclude <- c() # Add any File_Name value(s) - ex. c("Sample001_normalized_BDNF_15min_1_Cells_To_Analyze_157Gated_quantile.fcs")

# Plotting information
PLOT_HISTOGRAMS <- FALSE # TRUE / FALSE
plot_column_length <- 5 # Scalar for histogram plotting

OUTPUT.TYPE <- "pdf" # choose png or pdf - for histogram and heatmaps
OUTPUT_BASENAME <- paste0("Top", (1-Set_Threshold)*100, "percOfTimePoint", Threshold_TimePoint, "_" , "Thresholded_RepsAvg(minusOutlierReps)")

MarkerOrder <- c("Condition", "Cluster", "TimePoint", "Replicate", # Keep these the same
                 "pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La",
                 "pPLK1_113In", "H3K9ac_115In", "pSTAT3_153Eu", "pSTAT5_147Sm",
                 "pATM_159Tb", "yH2AX_165Ho", "p-p38_142Nd", "pNFKB_166Er", "p-cJUN_162Dy", 
                 "p-p53_151Eu", "CC3_173Yb",
                 "pPLCG-2_144Nd", "pSrc_154Sm"
                 ) 
            # c("Condition", "Cluster", "TimePoint", "Replicate", # Keep these the same
            # # "TrkB_s_158Gd", "TrkB_i_174Yb", "pTrkB_170Er", 
            # "pPLCG-2_144Nd", "pERK_171Yb", "pAKt_168Er", "pSrc_154Sm", "pCreb_176Lu",
            # "p-p38_142Nd", "pS6_175Lu", 
            # "pSTAT3_153Eu", "pSTAT5_147Sm",
            # "pGSK3b_139La", "pPLK1_113In", "H3K9ac_115In", "pNFKB_166Er",
            # # "p75NTR_161Dy", 
            # "p-cJUN_162Dy", "p-p53_151Eu",
            # "pATM_159Tb", "yH2AX_165Ho", "CC3_173Yb") 
        # or c() for un-ordered markers

# Define which markers to include in plots
SelectMarkersToPlot <- c("pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La",
                         "pPLK1_113In", "H3K9ac_115In", "pSTAT3_153Eu", "pSTAT5_147Sm",
                         "pATM_159Tb", "yH2AX_165Ho", "p-p38_142Nd", "pNFKB_166Er", "p-cJUN_162Dy", 
                         "p-p53_151Eu", "CC3_173Yb",
                         "pPLCG-2_144Nd", "pSrc_154Sm"
                         ) 


            # c(# "TrkB_s_158Gd",
            #                        # "TrkB_i_174Yb",
            #                        # "pTrkB_170Er",
            #                        "pPLCG-2_144Nd",
            #                        "pERK_171Yb",
            #                        "pAKt_168Er",
            #                        "pSrc_154Sm",
            #                        "pCreb_176Lu",
            #                        "p-p38_142Nd",
            #                        "pS6_175Lu",
            #                        "pSTAT3_153Eu",
            #                        "pSTAT5_147Sm",
            #                        "pGSK3b_139La",
            #                        "pPLK1_113In",
            #                        "H3K9ac_115In",
            #                        "pNFKB_166Er",
            #                        # "p75NTR_161Dy",
            #                        "p-cJUN_162Dy",
            #                        "p-p53_151Eu",
            #                        "pATM_159Tb",
            #                        "yH2AX_165Ho",
            #                        "CC3_173Yb"
            #                        )
        # "All" to plot all markers



#### Variables to that may need to be changed##################################
FILE_NUM <- paste0("filenums", ".csv")
PANEL_INFILE <- paste0("panel", ".csv")
METADATA_FILENAME <- paste0("metadata", ".csv")
EXPRS_MAT_INFILE1 <- "expression_matrix_analysis.csv"
EXPRS_MAT_INFILE2 <- "expression_matrix_other.csv"

CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
if (CLUSTERS_MERGED) {
  clusters_filename <- sub(".csv", "_merged.csv", clusters_filename)
}

THRESHOLD_CATEGORY <- "ThresholdCategory.csv"


#### Read in necessary files ##################################################
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

# Read in THRESHOLD_CATEGORY file
if (IsThereThresholdCategory) {
  ThresholdCategory <- read.table(THRESHOLD_CATEGORY, header=TRUE, sep=",",
                                  check.names=FALSE)
}

# Read in expression matrix and perform initial clean up of data
exprs_mat_IN1 <- fread(paste0(Original_WD, "/", EXPRS_MAT_INFILE1), 
                       stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)

exprs_mat_IN2 <- fread(paste0(Original_WD, "/", EXPRS_MAT_INFILE2), 
                       stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)
exprs_mat_IN_ALL <- cbind(exprs_mat_IN1, exprs_mat_IN2)


if (IsThereThresholdCategory) {
  AppendedFileNames <- c("File_Num", "Cluster", "ThresholdCategory")
  exprs_mat_IN <- cbind(FileNum, #Append FileNum for indexing
                        clusters_in, #Append clusters_in for indexing
                        ThresholdCategory, #Append ThresholdCategory for indexing
                        exprs_mat_IN_ALL[ , panel$Fixed_Param[which(panel$Plot==1)]]) #Only selects markers from panel$Plot==1
  colnames(exprs_mat_IN) <- c(AppendedFileNames, 
                              colnames(exprs_mat_IN[(length(AppendedFileNames)+1):length(colnames(exprs_mat_IN))]))
} else {
  AppendedFileNames <- c("File_Num", "Cluster")
  exprs_mat_IN <- cbind(FileNum, #Append FileNum for indexing
                        clusters_in, #Append clusters_in for indexing
                        exprs_mat_IN_ALL[ , panel$Fixed_Param[which(panel$Plot==1)]]) #Only selects markers from panel$Plot==1
  colnames(exprs_mat_IN) <- c(AppendedFileNames, 
                              colnames(exprs_mat_IN[(length(AppendedFileNames)+1):length(colnames(exprs_mat_IN))]))
}


# Define the FileNums to use based on WorkingSampleSet and subset exprs_mat_IN
WorkingSampleSet_FileNums <- metadata$File_Num[metadata[ , WorkingSampleSet_Identifier] %in% WorkingSampleSet &
                                                 !sub("_.*", "", metadata$Condition) %in% ConditionsToExclude]

exprs_mat <- exprs_mat_IN %>% 
  dplyr::filter(File_Num %in% WorkingSampleSet_FileNums) #Subset based on WorkingSampleSet



#### Prep the input data so that it can undergo normalization #################
# Change the linear data to arcsinh transformed data (factor = 5)
Transform_Factor <- 5
exprs_mat_Linear <- cbind(exprs_mat[ , AppendedFileNames], 
                          5 * sinh(exprs_mat[ , (length(AppendedFileNames)+1):ncol(exprs_mat)]))

# Remove any clusters that will not be necessary to analyze
if (length(ClustersToRemove) > 0) {
  exprs_mat_Linear <- exprs_mat_Linear[-which(exprs_mat_Linear$Cluster %in% ClustersToRemove), ]
  } else if (length(ClustersToRemove) == 0) {
    exprs_mat_Linear <- exprs_mat_Linear
    }



#### Create variables to be used later ########################################
Markers_Analysis <- panel$Fixed_Param[panel$Analysis == 1] #Markers used for clustering
Markers_Plot <- panel$Fixed_Param[panel$Plot == 1] #Markers used for plotting
ID_Markers  <-  panel$Fixed_Param[panel$Plot == 1 &
                                    #!panel$Fixed_Param %in% Markers_Analysis &
                                    !stringr::str_starts(panel$Fixed_Param, "p") & #Remove all markers starting with "p" (most signaling)
                                    !panel$Fixed_Param == "H3K9ac_115In" & #Remove all other non-ID markers
                                    !panel$Fixed_Param == "TrkB_s_158Gd" &
                                    !panel$Fixed_Param == "yH2AX_165Ho" &
                                    !panel$Fixed_Param == "CC3_173Yb" &
                                    !panel$Fixed_Param == "TrkB_i_174Yb" &
                                    !panel$Fixed_Param == "DNA_191Ir" &
                                    !panel$Fixed_Param == "DNA_193Ir"]
Signaling_Markers <- Markers_Plot[!Markers_Plot %in% ID_Markers &
                                    !Markers_Plot == "DNA_191Ir" &
                                    !Markers_Plot == "DNA_193Ir"]

# Define which markers will be excluded from the heatmaps
if (length(SelectMarkersToPlot) > 0) {
  Markers_OmitFromPlot <- Markers_Plot[!Markers_Plot %in% SelectMarkersToPlot | 
                                         !Markers_Plot %in% Signaling_Markers]
} else {
  Markers_OmitFromPlot <- Markers_Plot[!Markers_Plot %in% Signaling_Markers]
}

metadata_exprs_mat <- metadata[metadata$File_Num %in% WorkingSampleSet_FileNums, ]
NumOfChannels <- ncol(exprs_mat_Linear %>%  select(-all_of(AppendedFileNames)))
Sample_Conditions_Timepoints <- metadata_exprs_mat$TimePoint
Sample_Conditions_Timepoints_REORDERED <-  metadata_exprs_mat$TimePoint[order(metadata_exprs_mat$Reordered)]
Sample_Conditions_Replicate <- metadata_exprs_mat$Replicate
Sample_Conditions_Replicate_REORDERED <- metadata_exprs_mat$Replicate[order(metadata_exprs_mat$Reordered)]
Unique_Cluster <- sort(unique(exprs_mat_Linear$Cluster))
NumOfClusters <- length(Unique_Cluster)
Unique_FileNum <- sort(unique(exprs_mat_Linear$File_Num))
Unique_FileNum_REORDERED <- metadata_exprs_mat[order(metadata_exprs_mat$Reordered), 2]
Sample_Conditions <- metadata_exprs_mat$Condition
Sample_Conditions_UNIQUE <- unique(metadata_exprs_mat$Condition)
Sample_Conditions_REORDERED <- metadata_exprs_mat$Condition[order(metadata_exprs_mat$Reordered)]
Sample_Conditions_UNIQUE_REORDERED <- unique(metadata_exprs_mat$Condition[order(metadata_exprs_mat$Reordered)])

Threshold_Condition <- unique(metadata_exprs_mat$Condition[which(metadata_exprs_mat$TimePoint == Threshold_TimePoint &
                                                                   !grepl("Test", metadata_exprs_mat$Condition))])



#### Create folder for output files ###########################################
output_dir_All <- paste0("12.1-", # Script ID
                         format(Sys.time(), "%Y-%m-%d"), # Date info
                         paste0("_PercAboveThresh_", 
                                Set_Threshold*100, 
                                "thPerc",
                                "Of", sub("_.*", "", Threshold_Condition)), # Transformation info
                         "_Expt", paste(WorkingSampleSet, collapse = ","), # Expt ID
                         if (IsThereThresholdCategory) { # IsThereThresholdCategory
                           paste0("-", WhichToPlot)
                         },
                         if (length(ClustersToRemove) > 0) {
                           "_RemovedClusters"
                         }
)
dir.create(output_dir_All)
setwd(paste0(Original_WD, "/", output_dir_All))



#### Generate primary data frame to loop off of to generate heatmaps ##########
# Define the File_Num(s) to set the threshold against
Threshold_FileNum <- metadata$File_Num[metadata[ , WorkingSampleSet_Identifier] %in% WorkingSampleSet &
                                         is.na(metadata$TimePoint) == FALSE &
                                         metadata$TimePoint == Threshold_TimePoint]

# Take the average of the portion of the File_Num defined by Threshold_TimePoint
# Find the cutoff for each File_Num
if (IsThereThresholdCategory) {
  Threshold_CutOff_Avg <- exprs_mat_Linear %>% 
    dplyr::filter(Cluster %in% Unique_Cluster, 
                  File_Num %in% Threshold_FileNum, 
                  ThresholdCategory == WhichToPlot) %>% # Only from defined ThresholdCategory
    group_by(Cluster, File_Num) %>%
    summarize(across(where(is.numeric), 
                     ~ quantile(.x, probs = Set_Threshold)), 
              .groups = "drop") %>% 
  # Take the average of each column by cluster and rename File_Num, Cluster
    group_by(Cluster) %>%
    summarize(across(where(is.numeric), 
                     mean), 
              .groups = "drop") %>%
    mutate(File_Num = 0) %>%
    relocate(File_Num) # Move File_Num to the first column
  } else {
    Threshold_CutOff_Avg <- exprs_mat_Linear %>% 
      dplyr::filter(Cluster %in% Unique_Cluster, File_Num %in% Threshold_FileNum) %>%
      group_by(Cluster, File_Num) %>%
      summarize(across(where(is.numeric), 
                       ~ quantile(.x, probs = Set_Threshold)), 
                .groups = "drop") %>% 
      # Take the average of each column by cluster and rename File_Num, Cluster
      group_by(Cluster) %>%
      summarize(across(where(is.numeric), 
                       mean), 
                .groups = "drop") %>%
      mutate(File_Num = 0) %>%
      relocate(File_Num) # Move File_Num to the first column
  }



#### Identify the events in each FileNum that is above the defined threshold by CLUSTER and Marker (COLUMN) ####
# Create a named vector of thresholds for fast lookup
thresholds <- Threshold_CutOff_Avg %>%
  select(Cluster, everything()) %>%
  tibble::column_to_rownames("Cluster")
  
if (IsThereThresholdCategory) {
  exprs_mat_Linear_AboveThreshold <- exprs_mat_Linear %>%
    dplyr::filter(ThresholdCategory == WhichToPlot) %>% # Only from defined ThresholdCategory
    group_by(Cluster) %>% 
    mutate(across(
      where(is.numeric),
      ~ ifelse(.x > as.numeric(thresholds[as.character(unique(Cluster)), cur_column()]), .x, NA)
      ))
  } else {
    exprs_mat_Linear_AboveThreshold <- exprs_mat_Linear %>%
      group_by(Cluster) %>% 
      mutate(across(
        where(is.numeric),
        ~ ifelse(.x > as.numeric(thresholds[as.character(unique(Cluster)), cur_column()]), .x, NA)
      ))
  }


#### Calculate the percent of events in each FileNum that is above the defined threshold by CLUSTER and Marker (COLUMN) ####
# Calculate percent above threshold for each File_Num and Cluster
PercAboveThreshold <- exprs_mat_Linear %>%
  group_by(Cluster, File_Num) %>%
  summarise(across(
    where(is.numeric),
    ~ mean(.x > thresholds[as.character(unique(Cluster)), cur_column()]) * 100
  ), .groups = "drop")


# Add 'NA' values to any Cluster that is missing a FileNum - this will create a row for each Cluster represented across all File_Nums
PercAboveThreshold_PlusNA <- PercAboveThreshold %>%
  complete(Cluster, File_Num)


# Add a few more identifying columns
PercAboveThreshold_Indexed <- cbind(File_Name = rep(metadata_exprs_mat$File_Name, times = NumOfClusters), # Add File_Name column
                                    Sample_Condition = rep(Sample_Conditions, times = NumOfClusters), # Add Sample_Condition column
                                    Reordered = rep(metadata_exprs_mat$Reordered, times = NumOfClusters), # Add Reordered file numbers column
                                    File_Num = PercAboveThreshold_PlusNA$File_Num, # Add File_Num column
                                    Cluster = PercAboveThreshold_PlusNA$Cluster, # Add Cluster column
                                    TimePoints = rep(Sample_Conditions_Timepoints, times = NumOfClusters), # Add TimePoints column
                                    Replicate = rep(metadata_exprs_mat$Replicate, times = NumOfClusters), # Add Replicate column
                                    PercAboveThreshold_PlusNA[ , -which(colnames(PercAboveThreshold_PlusNA) %in% c("File_Num", "Cluster"))])# Reorganize PercAboveThreshold_PlusNA


# Remove any samples that are considered outliers
if (Exclude_FileNames) {
  PercAboveThreshold_Indexed_OutliersRemoved <- PercAboveThreshold_Indexed %>%
    dplyr::filter(!File_Name %in% FileNamesToExclude)
} else {
  PercAboveThreshold_Indexed_OutliersRemoved <- PercAboveThreshold_Indexed
}

# Reorder data frame
PercAboveThreshold_Indexed_OutliersRemoved_REORDERED <- PercAboveThreshold_Indexed_OutliersRemoved %>% 
  arrange(Cluster, Reordered)



### Plot each marker as a histogram and incorporate the respective Threshold_CutOff_Avg value
# Function to create histogram + density plot
create_histogram_plot <- function(data, marker, threshold, identifier, quantile_cutoff = NULL, inset = FALSE) {
  if (!is.null(quantile_cutoff)) {
    data <- data[data[[marker]] < quantile_cutoff, ]
  }
  
  p <- ggplot(data, aes(x = .data[[marker]])) +
    # Histogram of data
    geom_histogram(aes(y = after_stat(density), fill = Treatment),
                   color = "black", alpha = 0.25, bins = 100) +
    # Add density plot of data
    geom_density(aes(colour = Treatment), lwd = 1.2, linetype = 2, alpha = 0.25) +
    # Add line denoting respective threshold value from Threshold_CutOff_Avg
    geom_vline(xintercept = threshold, linetype = 'longdash') +
    # Adjust plot formating
    labs(x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.background = element_rect(color = "black", fill = "white", linewidth = 1))
  
  # Update formatting based on inset value
  if (inset) {
    p <- p + annotate("text", x = Inf, y = Inf, 
                      hjust = 1.05, vjust = 1.5,
                      size = 3, fontface = "bold", label = "95th Percentile Cutoff") +
      theme(legend.position = "none")
  } else {
    p <- p + ggtitle(marker)
  }
  
  # Update formatting based on Identifier
  if (identifier == "PercAboveThreshold") {
    p <- p + xlim(0, 100)
  }
  
  return(p)
}

plot_histograms_grid <- function(Input_Dataframe, Identifier, Extra_Columns, Include_Inset = TRUE) {
  # Loop through clusters to generate histograms for each marker
  for (j in Unique_Cluster) {
    # Preprocess data once per cluster
    cluster_data <- Input_Dataframe %>%
      filter(Cluster == j) %>%
      left_join(metadata_exprs_mat[, c("File_Num", "Treatment")], by = "File_Num")
    
    markers <- setdiff(colnames(Input_Dataframe), Extra_Columns)
    
    # Precompute 95th percentiles
    quantiles_to_plot <- sapply(markers, function(marker) unname(quantile(cluster_data[[marker]] %>% na.omit(), 0.95)))
    
    # Generate plots
    if (Include_Inset) {
      plot_list <- lapply(markers, function(marker) {
        full_plot <- create_histogram_plot(data = cluster_data %>% 
                                             select(Treatment, all_of(marker)) %>% 
                                             na.omit(), 
                                           marker = marker, 
                                           threshold = if (Identifier == "PercAboveThreshold") {
                                             (1-Set_Threshold) * 100
                                           } else {
                                             Threshold_CutOff_Avg[[j, marker]]
                                           },
                                           identifier = Identifier)
        partial_plot <- create_histogram_plot(data = cluster_data %>% 
                                                select(Treatment, all_of(marker)) %>% 
                                                na.omit(), 
                                              marker = marker, 
                                              threshold = if (Identifier == "PercAboveThreshold") {
                                                (1-Set_Threshold) * 100
                                              } else {
                                                Threshold_CutOff_Avg[[j, marker]]
                                              },
                                              identifier = Identifier,
                                              quantile_cutoff = quantiles_to_plot[[marker]], 
                                              inset = TRUE)
        
        ggdraw() +
          draw_plot(full_plot) +
          draw_plot(partial_plot, x = 0.39, y = 0.52, width = 0.4, height = 0.4)
      })
    } else {
      plot_list <- lapply(markers, function(marker) {
        full_plot <- create_histogram_plot(data = cluster_data, 
                                           marker = marker, 
                                           threshold = if (Identifier == "PercAboveThreshold") {
                                             (1-Set_Threshold) * 100
                                             } else {
                                               Threshold_CutOff_Avg[[j, marker]]
                                               },
                                           identifier = Identifier)
        ggdraw() +
          draw_plot(full_plot)
      })
    }
    names(plot_list) <- markers
    
    # Plot dimensions
    plot_width <- 5 * (length(plot_list) / plot_column_length)
    plot_height <- 5 * plot_column_length
    
    # Create output directory
    Histogram_dir <- paste0(format(Sys.time(), "%Y-%m-%d"),
                                     " Histogram_", Identifier,
                                     if (IsThereThresholdCategory) paste0("-", WhichToPlot))
    output_path <- file.path(Original_WD, output_dir_All, Histogram_dir)
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    
    # Pad plot list if needed (ex. if number of plots is not divisible by plot_column_length)
    remainder <- length(plot_list) %% plot_column_length
    if (remainder != 0) {
      NumberOfGridSlots <- ceiling(length(plot_list) / plot_column_length) * plot_column_length
      for (k in (length(plot_list) + 1):NumberOfGridSlots) {
        plot_list[[k]] <- ggplot() + theme_void()
      }
    }
    
    # Arrange and save grid
    grid_plot <- grid.arrange(grobs = plot_list,
                              nrow = ceiling(length(plot_list) / plot_column_length),
                              ncol = plot_column_length)
    
    ggsave(filename = file.path(output_path,
                                paste0("Histogram_", Identifier, "_Cluster", j,
                                       if (IsThereThresholdCategory) paste0("-", WhichToPlot),
                                       ".", OUTPUT.TYPE)),
           plot = grid_plot,
           width = plot_width, height = plot_height,
           units = "in", dpi = 300)
  }
}

if (PLOT_HISTOGRAMS) {
  # Plot ExprsMat ("RAW data") histograms for each marker by cluster
  plot_histograms_grid(Input_Dataframe = exprs_mat_Linear,
                       Extra_Columns = c("File_Num", "Cluster"),
                       Identifier = "ExprsMat")
  
  # Plot ValuesAboveThreshold ("RAW data" with Threshold_CutOff_Avg applied) histograms for each marker by cluster
  plot_histograms_grid(Input_Dataframe = exprs_mat_Linear_AboveThreshold,
                       Extra_Columns = c("File_Num", "Cluster"),
                       Identifier = "ValuesAboveThreshold")
  
  # Plot PercAboveThreshold ("RAW data" with Threshold_CutOff_Avg applied) histograms for each marker by cluster
  plot_histograms_grid(Input_Dataframe = PercAboveThreshold_Indexed_OutliersRemoved_REORDERED,
                       Identifier = "PercAboveThreshold",
                       Extra_Columns = c("File_Num", "Cluster", "File_Name", "Sample_Condition", "Reordered", "TimePoints", "Replicate"),
                       Include_Inset = FALSE)
}



#### Save cleaned up csv ######################################################
write.csv(x = PercAboveThreshold_Indexed_OutliersRemoved_REORDERED, 
          file = paste0("Expt", paste(WorkingSampleSet, collapse = ","), # Expt ID
                        paste0("_PercAboveThresh_", 
                               Set_Threshold*100, 
                               "thPerc",
                               "Of", sub("_.*", "", Threshold_Condition)), # Transformation info
                        "_", "RAWvalues", # Dataframe ID
                        if (IsThereThresholdCategory) { # IsThereThresholdCategory
                          paste0("-", WhichToPlot)
                        },
                        ".csv"))



#### Subset PercAboveThreshold data to exclude user-defined markers ###########
# Identify and remove columns based on user-defined Markers_OmitFromPlot
ColumnsToRemoveForPlot <- which(colnames(PercAboveThreshold_Indexed_OutliersRemoved_REORDERED) %in% Markers_OmitFromPlot)

ChannelsToKeepForHeatmap_PercAboveThresh <- PercAboveThreshold_Indexed_OutliersRemoved_REORDERED %>% 
  select(-c("File_Name", "Reordered", "File_Num", all_of(ColumnsToRemoveForPlot)))

# Calculate the average across replicates by Sample_Condition and Cluster
ChannelsToKeepForHeatmap_PercAboveThresh_RepsAvg <- ChannelsToKeepForHeatmap_PercAboveThresh %>%
  group_by(Sample_Condition, Cluster) %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

##### Could definitely write the average as a function, since it is used three times #####



#### Define normalization values for fold change calculation (i.e. average of replicates for Threshold_Condition) ####
# Define the timepoints needed for normalization
if (CtrlPerTimepoint == TRUE) {
  PostThreshold_NormTimepoints <- sort(unique(ChannelsToKeepForHeatmap_PercAboveThresh$TimePoints))
} else {
  PostThreshold_NormTimepoints <- Threshold_TimePoint
}

# Define the rows within ChannelsToKeepForHeatmap_PercAboveThresh to designate as normalization rows
Normalization_ROW <- c()
for (g in 1:length(PostThreshold_NormCondition)) {
  Normalization_ROW <- c(Normalization_ROW,
                         which(ChannelsToKeepForHeatmap_PercAboveThresh[1] == PostThreshold_NormCondition[g]))
}

# Defined the average of the normalization rows
Normalization_ROW_Avg <- ChannelsToKeepForHeatmap_PercAboveThresh[sort(Normalization_ROW), ] %>% 
  group_by(Sample_Condition, Cluster) %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame() %>%
  rename(NormTimepoint = TimePoints)



##### Calculate the fold change of the PercAboveThresh data relative to the Normalization_Row_Avg ####
ChannelsToKeepForHeatmap_FC <- ChannelsToKeepForHeatmap_PercAboveThresh %>%
  # Append a NormTimePoint column to ChannelsToKeepForHeatmap_PercAboveThresh
  mutate(NormTimepoint = if (CtrlPerTimepoint) TimePoints else PostThreshold_NormTimepoints) %>% 
  
  # Append Normalization_ROW_Avg to ChannelsToKeepForHeatmap_FC through NormTimePoint and add "_norm" to end of each appended column
  left_join(Normalization_ROW_Avg, 
            by = c("Cluster", "NormTimepoint" = "NormTimepoint"), 
            suffix = c("", "_norm")) %>%
  
  # Avoid division by zero
  mutate(across(ends_with("_norm"), 
                ~ ifelse(. == 0, NA, .))) %>% 
  
  # Calculate the Fold Change by dividing the ChannelsToKeepForHeatmap_FC columns by the "_norm" columns and add "FC_" to front of each appended column
  mutate(across(!c(Sample_Condition, TimePoints, Cluster, Replicate, NormTimepoint) & !ends_with("_norm"), 
                ~ . / get(paste0(cur_column(), "_norm")), 
                .names = "FC_{.col}")) %>%
  
  # Select only the columns that underwent Fold Change transformation and remove "FC_" from the column names
  select(Sample_Condition, Cluster, TimePoints, Replicate, starts_with("FC_")) %>% 
  rename_with(~ sub("^FC_", "", .), starts_with("FC_"))


# Calculate the average fold change across replicates by Sample_Condition and Cluster
ChannelsToKeepForHeatmap_FC_RepsAvg <- ChannelsToKeepForHeatmap_FC %>%
  group_by(Sample_Condition, Cluster) %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()



#### Calculate the log2(Fold Change) of the ChannelsToKeepForHeatmap_FC data ####
ChannelsToKeepForHeatmap_log2FC <- ChannelsToKeepForHeatmap_FC %>%
  mutate(across(5:ncol(ChannelsToKeepForHeatmap_FC), log2)) # Mutate across non-identifying columns

# Calculate the average log2(fold change) across replicates by Sample_Condition and Cluster
ChannelsToKeepForHeatmap_log2FC_RepsAvg <- ChannelsToKeepForHeatmap_log2FC %>%
  group_by(Sample_Condition, Cluster) %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()



#### Save ChannelsToKeepForHeatmap dataframes #################################
# Put ChannelsToKeepForHeatmap dataframes in a list
ChannelsToKeepForHeatmap_List <- list(RAWvalues = ChannelsToKeepForHeatmap_PercAboveThresh,
                                      RAWvalues_RepsAvg = ChannelsToKeepForHeatmap_PercAboveThresh_RepsAvg,
                                      FC = ChannelsToKeepForHeatmap_FC,
                                      FC_RepsAvg = ChannelsToKeepForHeatmap_FC_RepsAvg,
                                      log2FC = ChannelsToKeepForHeatmap_log2FC,
                                      log2FC_RepsAvg = ChannelsToKeepForHeatmap_log2FC_RepsAvg)

# Loop through each ChannelsToKeepForHeatmap dataframe and save
for (i in 1:length(ChannelsToKeepForHeatmap_List)) {
  ChannelsToKeepForHeatmap_FileName <- paste0("Expt", paste(WorkingSampleSet, collapse = ","), # Expt ID
                                              paste0("_PercAboveThresh_", 
                                              Set_Threshold*100, 
                                              "thPerc",
                                              "Of", sub("_.*", "", Threshold_Condition)), # Transformation info
                                       if (grepl("NC", paste(PostThreshold_NormCondition, collapse = ","))) {
                                         "_NormToNC"
                                         } else {
                                           "_NormTo0hr"
                                           }, # Normalization info
                                       "_", names(ChannelsToKeepForHeatmap_List[i]), # ChannelsToKeepForHeatmap dataframe ID
                                       ".csv")
  
  write.csv(x = ChannelsToKeepForHeatmap_List[[i]], 
            file = ChannelsToKeepForHeatmap_FileName)
  }


    
#### Put all dataframes into a list that can be looped through later ##########
ChannelsToKeepForHeatmap_ListOfDatasets <- list(ChannelsToKeepForHeatmap_PercAboveThresh,
                                                ChannelsToKeepForHeatmap_PercAboveThresh_RepsAvg,
                                                ChannelsToKeepForHeatmap_FC,
                                                ChannelsToKeepForHeatmap_FC_RepsAvg,
                                                ChannelsToKeepForHeatmap_log2FC,
                                                ChannelsToKeepForHeatmap_log2FC_RepsAvg)
    # Note to R learners that sometimes forget list indexing (like myself):
    #       To access a specific dataset: ChannelsToKeepForHeatmap_ListOfDatasets[[1]]

DataSetNames <- c("ChannelsToKeepForHeatmap_PercAboveThresh",
                  "ChannelsToKeepForHeatmap_PercAboveThresh_RepsAvg",
                  "ChannelsToKeepForHeatmap_FC",
                  "ChannelsToKeepForHeatmap_FC_RepsAvg",
                  "ChannelsToKeepForHeatmap_log2FC",
                  "ChannelsToKeepForHeatmap_log2FC_RepsAvg")



#### Other variables that don't necessarily need to be changed ################
TimepointsPerCond_All <- metadata_exprs_mat %>%
  mutate(Condition = sub("_.*", "", metadata_exprs_mat$Condition)) %>%
  group_by(Condition) %>% 
  summarize(TPs = length(unique(TimePoint)))

TimepointsPerCond <- max((TimepointsPerCond_All)[ , "TPs"])



#### For testing of the function ##############################################
# input_dataset <- ChannelsToKeepForHeatmap_ListOfDatasets[[6]]
# input_dataset_name <- substr(DataSetNames[6], nchar("ChannelsToKeepForHeatmap_"), nchar(DataSetNames[6]))
# input_dataset <- c()
# input_dataset_name <- c()



#### Define function to organize data and generate heatmaps ###################
HeatmapsFromDataset_Function <- function(input_dataset, input_dataset_name) {
  
  ### Prep data for Heatmap functions ###
  # Change any Inf values to NA, then remove PostThreshold_NormCondition values
  current_dataset <- data.frame(lapply(input_dataset, function(x) {
    x[is.infinite(x)] <- NA
    return(x)
    })) %>% 
  # Remove PostThreshold_NormCondition from the dataframe
    slice(-which(.data$Sample_Condition %in% PostThreshold_NormCondition))

  # Add 'NA' values to any time point that is missing data - this will create a row for each unique time point represented across all conditions
  if (length(ClusterOrder) == 0) {
    ClusterOrder <- c(1:NumOfClusters)
  }
  
  if (length(ClustersToRemove) > 0) {
    ClusterOrder <- ClusterOrder[!ClusterOrder %in% ClustersToRemove]
  }
  
  current_dataset_PlusNA <- current_dataset %>%
  # Remove time info from Sample_Condition for complete() function
    mutate(Sample_Condition = sub("_.*", "", Sample_Condition)) %>% 
  # 'Complete' the data frame based on Sample_Condition, Cluster, and unique TimePoints
    complete(Sample_Condition, Cluster, TimePoints = unique(TimePoints)) %>%
  # Add time info back to Sample_Condition and make cluster column a factor that is ordered by ClusterOrder
    mutate(Sample_Condition = paste0(Sample_Condition, "_", TimePoints, "hr"), 
           Cluster = factor(.data$Cluster, levels = ClusterOrder))
  
  # Reorder Cluster by its factor levels based on ClusterOrder (if provided) and the Marker columns based on MarkerOrder (if provided)
  if (length(SelectMarkersToPlot) > 0) {
    MarkerOrder <- gsub("-", ".", # adjust to match current_dataset_PlusNA colnames
                        MarkerOrder[MarkerOrder %in% c("Sample_Condition", "Cluster", "TimePoints", "Replicate", 
                                                       SelectMarkersToPlot)]) # Only include indexing and SelectMarkersToPlot markers
  } else {
    MarkerOrder <- gsub("-", ".", # adjust to match current_dataset_PlusNA colnames
                        MarkerOrder)
  }
  
  if (length(MarkerOrder) > 0) {
    current_dataset_PlusNA <- current_dataset_PlusNA[order(current_dataset_PlusNA$Cluster) , match(MarkerOrder, colnames(current_dataset_PlusNA))]
    
    current_dataset_PlusNA[ , colnames(current_dataset_PlusNA) %in% MarkerOrder]
  } else {
    current_dataset_PlusNA <- current_dataset_PlusNA[order(current_dataset_PlusNA$Cluster), ]
  }
  
  ### Define useful naming variables from the cleaned up dataset
  ConditionsNames <- unique(sub("(_[^_]+)$", "", current_dataset_PlusNA$Sample_Condition))
  MarkerNames <- unique(sub("(_[^_]+)$", "", colnames(select(current_dataset_PlusNA, 
                                                             -c("Sample_Condition", 
                                                                "TimePoints", 
                                                                "Replicate", 
                                                                "Cluster")))))  
  
  ### Shift data from each time point into a new and adjacent column 
  DataForHeatmap_temp <- current_dataset_PlusNA %>%
  # Remove the Replicate column
    select(-Replicate) %>%
  # Modify columns 
    mutate(Sample_Condition = sub("_.*", "", Sample_Condition), # Remove time info
           TimePoints = as.character(TimePoints)) %>% # Change type to character
  # Pivot data based on Sample_Condition and Cluster into new columns with names to reflect TimePoints
    pivot_wider(
      id_cols = c(Sample_Condition, Cluster),
      names_from = TimePoints,
      values_from = where(is.numeric)
    )
  
  ### Clean up matrices that will be inputs for a heatmap function: Group rows by Condition or by Cluster
  input_matrix_CondGrouped <- DataForHeatmap_temp %>% 
  # Order based on Sample_Condition levels
    mutate(Sample_Condition = factor(Sample_Condition, 
                                     levels = ConditionsNames)) %>% 
    arrange(Sample_Condition) %>% 
  # Remove columns used for identification
    select(-Sample_Condition, -Cluster) %>% 
  # Output a matrix - Heatmap function requirement
    as.matrix()
  
  input_matrix_ClustGrouped <- DataForHeatmap_temp %>% 
    # Order based on Cluster levels
    arrange(Cluster) %>% 
    # Remove columns used for identification
    select(-Sample_Condition, -Cluster) %>% 
    # Output a matrix - Heatmap function requirement
    as.matrix()
  
  ### Define the heatmap functions to organize by Condition or by Cluster ###
  # Write function to generate heatmap vertically split by condition or cluster across time with signaling markers ordered vertically
  HEATMAP_FUNCTION <- function(mat, grouping, markers, timepoints, conditions, label_cells) {
    stopifnot(ncol(mat) == length(markers) * length(timepoints))
    stopifnot(nrow(mat) == length(conditions) * length(unique(DataForHeatmap_temp$Cluster)))
    
    # Use "grouping" to check for incompatibilities in data organization and how to organize the heatmap
    if (grouping == "Cluster") {
      ROW_LABELS <- rep(conditions, times = length(unique(DataForHeatmap_temp$Cluster)))
      ROW_TITLES  <- unique(DataForHeatmap_temp$Cluster)
      ROW_SPLIT <- rep(unique(DataForHeatmap_temp$Cluster), each = length(conditions))
      ROW_TITLE_ROT <- 0
    } else if (grouping == "Condition") {
      ROW_LABELS <- rep(unique(DataForHeatmap_temp$Cluster), times = length(conditions))
      ROW_TITLES <- conditions
      ROW_SPLIT <- rep(conditions, each = length(unique(DataForHeatmap_temp$Cluster)))
      ROW_TITLE_ROT <- 90
    }
    
    # Create a variable to store generated heatmaps for each marker
    heatmaps <- list()
    
    # Loop through each marker to generate and store and individual heatmap for each marker in "heatmaps"
    for (i in seq_along(markers)) {
      # Define the dataset for the active marker
      marker <- markers[i]
      col_idx <- ((i - 1) * length(timepoints) + 1):(i * length(timepoints))
      sub_mat <- mat[, col_idx, drop = FALSE]
      
      # Define a function for assigning color to the cells of the heatmap
      custom_color_fun <- function(x, NAME) {
        
        if (grepl("_Perc", NAME) == TRUE) { # Percent above threshold data
          colorRamp2(c(0, max(x, na.rm = TRUE)), 
                     c("blue", "white", "red"))
          
        } else if (grepl("_FC", NAME) == TRUE) { # Fold change data
          if (max(x, na.rm = TRUE) >= 1.5) {
            colorRamp2(c(0, 0.6667, 1.4999, max(x, na.rm = TRUE)), 
                       c("blue", "white", "white", "red"))
          } else {
            colorRamp2(c(0, 0.6667, 1.4999), 
                       c("blue", "white", "white"))
          }
          
        } else if (grepl("_log2FC", NAME) == TRUE) { # log2(Fold change) data
          if (max(x, na.rm = TRUE) >= log2(1.5) &
              min(x, na.rm = TRUE) <= -log2(1.5)) {
            colorRamp2(c(-max(x, na.rm = TRUE), -0.584, 0.584, max(x, na.rm = TRUE)), 
                       c("blue", "white", "white", "red"))
          } else if (max(x, na.rm = TRUE) <= log2(1.5) &
                     min(x, na.rm = TRUE) <= -log2(1.5)) {
            colorRamp2(c(min(x, na.rm = TRUE), -0.584, 0.584, -min(x, na.rm = TRUE)), 
                       c("blue", "white", "white", "red"))
          } else if (max(x, na.rm = TRUE) >= log2(1.5) &
                     min(x, na.rm = TRUE) >= -log2(1.5)) {
            colorRamp2(c(-max(x, na.rm = TRUE), -0.584, 0.584, max(x, na.rm = TRUE)), 
                       c("blue", "white", "white", "red"))
          } else if (max(x, na.rm = TRUE) <= log2(1.5) &
                     min(x, na.rm = TRUE) >= -log2(1.5)) {
            colorRamp2(c(-log2(1.5), log2(1.5)), 
                       c("white", "white"))
          }
        }
      }
      
      # Generate heatmap for the active marker
      hm <- Heatmap(
        sub_mat,
        name = marker,
        
        cluster_rows = FALSE,
        row_title = ROW_TITLES,
        row_title_side = "left",
        row_title_gp = gpar(fontsize = 12, fontface = "bold"), # "plain", "bold", "italic", "oblique", and "bold.italic"
        row_title_rot = 90,
        show_row_names = TRUE,
        row_labels = ROW_LABELS,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8, fontface = "plain"), # "plain", "bold", "italic", "oblique", and "bold.italic"
        
        row_split = ROW_SPLIT,
        
        cluster_columns = FALSE,
        column_title = marker,
        column_title_rot = 45,
        show_column_names = TRUE,
        column_labels = timepoints,
        column_names_side = "top",
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 12, fontface = "bold"), # "plain", "bold", "italic", "oblique", and "bold.italic"
        
        col = custom_color_fun(sub_mat, input_dataset_name), # Scales the color legend
        show_heatmap_legend = FALSE, # Hide individual legends
        
        cell_fun = if (label_cells) {
          function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", sub_mat[i, j]), x, y, gp = gpar(fontsize = 12)) # Draw the cell value
            grid.rect(x, y, width, height, gp = gpar(col = "darkgray", lwd = 1.5, fill = NA)) # Draw cell border
          }
        } else {
          function(j, i, x, y, width, height, fill) {
            grid.rect(x, y, width, height, gp = gpar(col = "darkgray", lwd = 1.5, fill = NA)) # Draw cell border
          }
        }
        ,
        
        border = TRUE
      )
      
      heatmaps[[i]] <- hm
    }
    
    # Combine heatmaps for all markers
      Reduce(`+`, heatmaps) # vertical: `%v%` / horizontal: `+`
  }
  
  ### Generate a list of matrices to run throught the HEATMAP_FUNCTION
  input_matrix_list <- list(Condition = input_matrix_CondGrouped, 
                            Cluster = input_matrix_ClustGrouped)
  
  ### Use the HEATMAP_FUNCTION to generate and save the heatmaps
  input_matrix_heatmap <- c()
  for (s in 1:length(input_matrix_list)) {
    for (c in c(TRUE, FALSE)) {
    input_matrix_heatmap <- HEATMAP_FUNCTION(mat = input_matrix_list[[s]],
                                             grouping = names(input_matrix_list)[s],
                                             MarkerNames, # markers
                                             timepoints = sort(unique(current_dataset_PlusNA$TimePoints)),
                                             conditions = ConditionsNames,
                                             label_cells = c)
    
    print(input_dataset_name) # Print name as sanity check that the inputs are as expected
    
    if (c) {
      cells_labeled <- "_labeled"
    } else {
      cells_labeled <- ""
    }
    
    if (OUTPUT.TYPE == "png") {
      png(filename = paste0("Heatmap_By", names(input_matrix_list)[s], input_dataset_name, cells_labeled, ".png"),
          width = Heatmap_Width, 
          height = Heatmap_Height)
      draw(input_matrix_heatmap)
    } else if (OUTPUT.TYPE == "pdf") {
      pdf(file = paste0("Heatmap_By", names(input_matrix_list)[s], input_dataset_name, cells_labeled, ".pdf"),
          width = Heatmap_Width / 100, 
          height = Heatmap_Height / 100)
      draw(input_matrix_heatmap)
    } else {
      print("Please select output file type")
    }
    dev.off()
    }
  }
  print("End Function")
}


#### Use the above function to plot and save heatmaps #########################
MarkersInHeatmap <- colnames(ChannelsToKeepForHeatmap_ListOfDatasets[[1]])[!colnames(ChannelsToKeepForHeatmap_ListOfDatasets[[1]]) %in% 
                                                                             c("Sample_Condition", "Cluster", "TimePoints", "Replicate")]
ClustersInHeatmap <- unique(ChannelsToKeepForHeatmap_ListOfDatasets[[1]]$Cluster)
ConditionsInHeatmap <- unique(sub("_.*", "", ChannelsToKeepForHeatmap_ListOfDatasets[[1]]$Sample_Condition))

Heatmap_Width <- (30 * length(MarkersInHeatmap)) * TimepointsPerCond
Heatmap_Height <- (30 * length(ClustersInHeatmap)) * length(ConditionsInHeatmap[!ConditionsInHeatmap %in% "NoTreat"])


DatasetsToPlot <- c(
  # "ChannelsToKeepForHeatmap_FC",
  "ChannelsToKeepForHeatmap_FC_RepsAvg",
  # "ChannelsToKeepForHeatmap_log2FC",
  "ChannelsToKeepForHeatmap_log2FC_RepsAvg"
  )

for (t in which((DataSetNames %in% DatasetsToPlot) == TRUE)) {
  HeatmapsFromDataset_Function(ChannelsToKeepForHeatmap_ListOfDatasets[[t]], # input_dataset
                               substr(DataSetNames[t], nchar("ChannelsToKeepForHeatmap_"), nchar(DataSetNames[t]))) # input_dataset_name
}


#### Generate a generic legend template and save it ###########################
Generic_Legend_values <- c(-7:7) # 1:15

Generic_Legend <- Heatmap(
  Generic_Legend_values,
  name = "Legend",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  col = colorRamp2(c(-7, 0, 7), 
                   c("blue", "white", "red")),
  show_heatmap_legend = FALSE,
  border = TRUE
)

if (OUTPUT.TYPE == "png") {
  png(filename = paste0("Generic_Heatmap_Legend", ".png"),
      width = 50, 
      height = 500)
  draw(Generic_Legend)
} else if (OUTPUT.TYPE == "pdf") {
  pdf(file = paste0("Generic_Heatmap_Legend", ".pdf"),
      width = 0.5, # in inches
      height = 5) # in inches
draw(Generic_Legend)
} else {
  print("Please select output file type")
}
dev.off()


#### Return to original directory once everything is completed ################
setwd(Original_WD)


print("Complete! Have fun staring at your new plot!")

