
# Author: Jonathon Sewell
# Date: 2025-02-28
# Version: 9
# Update:
# (1) Fixed bugs additional bugs


# Purpose: 
# Generate dotplots of percent above a threshold (defined by the nth percent of a user-defined condition/time) 
# and colored by Fold Change relative to user-defined control




# # PURPOSE: Generate a Heatmap of the nth percent of all samples above a user-defined threshold by cluster based on a user-defined parameter
# # INPUT:  _ "filenums.csv"    _ "panel.csv"   _ "cluster_RX_exprs_means.csv" (X = clustering round)    
# #         _ "expression_matrix_analysis.csv"  _ "expression_matrix_other.csv"
# #         _ "metadata.csv" 
# #     ***NOTE: make sure "metadata.csv" includes the following columns: 
# #            _ File_Name   _ File_Num    _ Reordered       _ BarcodeSet  _ SampleCollection	
# #            _ Expt	      _ Condition	  _ Condition_Only  _ TimePoint   _ Replicate
# 
# # OUTPUT: _ Complete dataset thresholded by the user-defined parameters ("[...]_RAWvalues.csv")
# #         _ Fold Change data from the "RAWvalues" dataset relative to a user-defined condition(s) ("_FC.csv", _FC_RepsAvg.csv)
# #         _ log2(Fold change) data of the "FC" dataset ("_log2FC.csv", _log2FC_RepsAvg.csv)
# #         _ Heatmap of Fold Change data organized by cluster or by treatment condition (Heatmap_ByCluster_FC_RepsAvg.png, Heatmap_ByCondition_FC_RepsAvg.png)
# #         _ Heatmap of log2(FC) data organized by cluster or by treatment condition (Heatmap_ByCluster_log2FC_RepsAvg.png, Heatmap_ByCondition_log2FC_RepsAvg.png)
# 
# 
# # NOTES:  _ This script was originally modified from 06.3_Plot_Expression_HeatmapByCluster_BySample.R from Pipeline2
# #         _ To normalize data to a specified condition, the exprs_matrix data undergoes a Linear Transform from the input Arcsinh Transform data


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

#### Variables to Change #######################################################
### Set up variables used to input data ========================================
  ## Define how to initially subset the input dataframe based on larger-scale metadata
    WorkingSampleSet_Identifier <- "Expt" # i.e. BarcodeSet, SampleCollection, Expt
    WorkingSampleSet <- c(1, 2, 3, 4) # Select a sample set from the complete data set to work with (i.e. metadata$SampleCollection == 1)
    # Cheatsheet: NoTreat = "1" / ContDepr, Rescue, BNDF = "2, 3, 4" / K252a expt = "5" / Stock = "6" / Universal = "7"
      
  ## Define whether ThresholdCategory should be considered (only if 13.1_ExtractByMarkerExpression_HiLo_v3.R already used)
    IsThereThresholdCategory <- TRUE # TRUE / FALSE - TRUE if using output from 13.1_ExtractByMarkerExpression_HiLo_v3.R
    WhichToPlot <- c(11, 14, 41, 44) # hi / med / lo


### Set up variables for computational analysis ###
  ## Define whether to perform PseudoBulk analysis or cluster analysis
    DataBy <- "PseudoBulk" # "Clusters" / "PseudoBulk"
  
  ## Define how to set the threshold and on which timepoint to base this threshold
    Set_Threshold <- 0.95 # 1st Quartile = 0.25 / 2nd Quartile (Median) = 0.5 / 3rd Quartile = 0.75
    Threshold_TimePoint <- 0 # Must match value from metadata$TimePoint
  
  ## Define how to threshold the data (against time=0hr, or against another condition)
    NormalizeAgainst <- "0hr" # "0hr" / "Control Condition"
    # "Control Condition" if normalization is against a control of the same time point (ex. Expt at 1hr compared to Ctrl at 1hr); 
    # "0hr" if normalization is against Threshold_TimePoint (ex. Expt at 1hr compared to Ctrl at 0hr)

  ## Define Cluster information
    CLUSTER_ROUND <- 1 # Need to change this for the level of subsetting/subclustering!
    CLUSTERS_MERGED <- FALSE
    ClusterOrder <- c() # c() if ClusterOrder is to go un-defined
    ## Define order for all clusters, or order for clusters to be plotted 
    ## (if all are not included, only those in ClusterOrder will be plotted)
    # Progenitor: 1, 15
    # Unknown: 5, 10, 11, 14, 16, 19
    # Non-Neuronal: 3, 9, 17, 18, 20
    # Neuronal: 2, 4, 6, 7, 8, 12, 13
    # Neuron lineage: 1, 7, 2, 13, 4
    
    ## Other Examples:
    # c(1, 7, 2, 12, 13, 6, 8, 4,
    #   15, 3, 9, 17, 18, 20,
    #   5, 10, 11, 14, 16, 19
    #   )
    # rev(c(1, 7, 2, 13, 4)) # rev(c(3, 9, 15, 17, 18)) 
    # rev(c(3, 9, 15, 17, 18,
    #       1, 7, 2, 13, 4))
    
    
### Plotting information =======================================================
  ## Define output file type for all plots
    OUTPUT.TYPE <- "pdf" # choose png or pdf - for histogram and heatmaps
    
  ## Define whether to plot histograms of raw data
    PLOT_HISTOGRAMS <- FALSE # TRUE / FALSE
    
  ## Define which conditions to include on the plot
    TreatmentToAnalyze <- c("NoTreat", "BDNF")
    # Expt 1: "NoTreat"
    # Expt 2: "ContDepr"
    # Expt 3: "CompMed"
    # Expt 4: "BDNF"
    # Expt 5: "NoTreatTest" / BDNFTest" / "K252aTest" / "ContDeprTest"
    # Expt 6: "Stock"
    # Expt 7: "Universal"
    
  ## Define DotPlot variables
    # Define column number for both histogram and DotPlots
    plot_column_length <- 3
    # Scalar to adjust plot size
    SizePerPlot <- 5
    # Dimensions for single plots
    Plot_Width <- 3.5
    Plot_Height <- 6
    # Define whether to output individual plots, or a grid of plots
    PlotAsGrid <- TRUE # TRUE / FALSE
    # Include numerical values on plot
    ValuesOnPlot <- TRUE # TRUE / FALSE
    # Define whether the DotSize is standard across all markers, or varies by marker
    StandardizeDotSize <- "AcrossPlotting" # "AcrossDataset" (standardized across entire dataset) / "AcrossPlotting" (standardized across plotted data) / FALSE
    # Define whether the DotColor is standard across all markers, or varies by marker
    StandardizeDotColor <- FALSE # "AcrossDataset" (standardized across entire dataset) / "AcrossPlotting" (standardized across plotted data) / FALSE
    # Color variables
    ColorBy <- "FoldChange" # "FoldChange" / "log2FC"
    ColorToUse <- "A" # With viridis, A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo
    
  ## Define Marker information
    # Order of markers on output plot
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
    
    
### Define anything to be excluded from the analysis and/or plotting ===========
    
  # Are there any clusters that need to be removed from analysis (i.e. excluded from quantification, not just plotting)?
    ClustersToExlude <- c() 

  # Define how to subset further based on condition (Note: only need to include base condition, no timepoint information)
    ConditionsToExclude <- c() # c() if none are to be excluded
      # Expt 1: "NoTreat"
      # Expt 2: "ContDepr"
      # Expt 3: "CompMed"
      # Expt 4: "BDNF"
      # Expt 5: "NoTreatTest" / BDNFTest" / "K252aTest" / "ContDeprTest"
    
  # Define if there are any File_Name values to exclude (considered to be outliers or otherwise)
    FileNamesToExclude <- c() # Add any File_Name value(s) - ex. c("Sample001_normalized_BDNF_15min_1_Cells_To_Analyze_157Gated_quantile.fcs")

    
### Variables that are altered automatically based on user-defined variables
  ## If plotting by PseudoBulk, ClusterOrder is reset so all clusters are included in analysis
    if (DataBy == "PseudoBulk") {
      ClusterOrder <- c()
      }

  ## If plotting all markers from PseudoBulk data onto one plot, change PlotAsGrid
    if (DataBy == "PseudoBulk" & !isFALSE(StandardizeDotColor) & !IsThereThresholdCategory) {
      PlotAsGrid <- FALSE
    }
    
    
  ## Output what condition(s) to compare to for fold change calculations
    PostThreshold_NormCondition <- if (NormalizeAgainst == "Control Condition") {
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
    } else if (NormalizeAgainst == "0hr") {
      if (mean(WorkingSampleSet %in% c(1, 2, 3, 4)) > 0) {
        PostThreshold_NormCondition <- c("NoTreat_0hr")
      } else if (mean(WorkingSampleSet %in% c(5)) == 1) {
        PostThreshold_NormCondition <- c("NoTreat_0hr") # Make sure this matches metadata.csv
      } else {
        print("Make sure to include Expt 1")
        PostThreshold_NormCondition <- c("NoTreat_0hr")
      }
    } else {
      stop("Check NormalizeAgainst variable is defined")
    }
   
  ## If IsThereThresholdCategory == TRUE, ensure that Original_WD refers to directory generated by 13_ExtractByMarkerExpression.R (or a version of it)
    if (IsThereThresholdCategory) {
      MarkerThresholdedOn <- str_extract(basename(Original_WD), 
                                         "(?<=ExtractByMarkerExpression_)[^_]+_[^_]+")
      if (is.na(MarkerThresholdedOn)) {
        stop("Copy this script into output folder from 13_ExtractByMarkerExpression containing the necessary 'ThresholdCategory.csv' file")
      }
    } 



#### Variables that likely do not need to be changed ###########################
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

THRESHOLD_CATEGORY_filename <- "ThresholdCategory.csv"


#### Read in necessary files ###################################################
### Read in file numbers associated with the expression file
  FileNum <- read.table(FILE_NUM, header=TRUE, sep=",",
                        check.names=FALSE)

### Read in panel (will show which expression markers will be included)
  panel <- read.table(PANEL_INFILE, header=TRUE, sep=",",
                      check.names=FALSE)

### Read in metadata
  metadata <- read.table(METADATA_FILENAME, header=TRUE, sep=",",
                         check.names=FALSE)

### Read in cluster file
  clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                            check.names=FALSE)
  ## Define which clusters to remove from plot based on ClusterOrder values and ClustersToExlude values
  if (!is.null(ClusterOrder)) {
    ClustersToRemoveFromPlot <- sort(unique(clusters_in[[1]])[!unique(clusters_in[[1]]) %in% ClusterOrder])
    } else {
      ClustersToRemoveFromPlot <- ClustersToExlude
      ClusterOrder <- sort(unique(clusters_in[[1]])[!unique(clusters_in[[1]]) %in% ClustersToExlude])
      }

### Read in THRESHOLD_CATEGORY_filename file
  if (IsThereThresholdCategory) {
    ThresholdCategory <- read.table(THRESHOLD_CATEGORY_filename, header=TRUE, sep=",",
                                    check.names=FALSE)
    if (is.null(WhichToPlot)) {
      WhichToPlot <- sort(unique(ThresholdCategory[[1]]))
    }
  }

### Read in expression matrix and perform initial clean up of data
  exprs_mat_IN1 <- fread(paste0(Original_WD, "/", EXPRS_MAT_INFILE1), 
                         stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)
  
  exprs_mat_IN2 <- fread(paste0(Original_WD, "/", EXPRS_MAT_INFILE2), 
                         stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)
  exprs_mat_IN_ALL <- cbind(exprs_mat_IN1, exprs_mat_IN2)

  ## Combined relevant dataframes
  if (IsThereThresholdCategory) {
    AppendedFileNames <- c("File_Num", "Cluster", "ThreshCat")
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
  
  ## Define the FileNums to include based on WorkingSampleSet
  WorkingSampleSet_FileNums <- metadata$File_Num[metadata[ , WorkingSampleSet_Identifier] %in% WorkingSampleSet &
                                                   !sub("_.*", "", metadata$Condition) %in% ConditionsToExclude]
  ## Subset exprs_mat_IN on WorkingSampleSet_FileNums
  exprs_mat <- exprs_mat_IN %>% 
    dplyr::filter(File_Num %in% WorkingSampleSet_FileNums) #Subset based on WorkingSampleSet



#### Prep the input data so that it can undergo normalization #################
### Change the arcsinh transformed data to linear data
  Transform_Factor <- 5
  exprs_mat_Linear <- exprs_mat %>% 
    mutate(across(colnames(exprs_mat)[!colnames(exprs_mat) %in% AppendedFileNames],
                  ~ Transform_Factor * sinh(.x)
                  ))

  ## Remove any clusters that will not be included in the computational analysis
  if (length(ClustersToExlude) > 0) {
    exprs_mat_Linear <- exprs_mat_Linear[-which(exprs_mat_Linear$Cluster %in% ClustersToExlude), ]
  } else if (length(ClustersToExlude) == 0) {
    exprs_mat_Linear <- exprs_mat_Linear
  }



#### Create variables to be used later #########################################
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
if (length(SelectMarkersToPlot) > 0 & paste0(SelectMarkersToPlot, collapse = "") != "All") { # User-defined markers have been selected
  Markers_OmitFromPlot <- Markers_Plot[!Markers_Plot %in% SelectMarkersToPlot | 
                                         !Markers_Plot %in% Signaling_Markers]
  } else if (SelectMarkersToPlot == "All") { # All markers in Markers_Plot will be used
    Markers_OmitFromPlot <- NULL
    SelectMarkersToPlot <- Markers_Plot
    MarkerOrder <- c(MarkerOrder, Markers_Plot[!Markers_Plot %in% MarkerOrder])
    } else { # No user-defined markers have been selected, so only plot Signaling_Markers
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



#### Create folder for output files ############################################
output_dir_All <- paste0("12.1.1-", # Script ID,
                         format(Sys.time(), "%Y-%m-%d"), # Date info
                         "_DotPlot_", # Plot info
                         DataBy, "_", # PseudoBulk or by cluster
                         paste0(Set_Threshold*100, 
                                "thPerc",
                                "Of", sub("_.*", "", Threshold_Condition), "_"), # Threshold info
                         "Expt", paste(WorkingSampleSet, collapse = ","), # Expt ID
                         "_SizeByFracAbund_", # Size info
                         paste0("ColorBy", ColorBy, "_"), # Color info
                         ifelse(test = length(ClustersToExlude) > 1, 
                                yes = paste0("Excl-Cl", paste0(ClustersToExlude, collapse = ","), "_"),
                                no = "InclAllClusts_"), # Clusters included
                         ifelse(test = length(intersect(SelectMarkersToPlot, Markers_Plot)) == sum(!Markers_Plot %in% Markers_Analysis),
                                yes = "AllMarkers",
                                no = "SelectMarkers"),
                         ifelse(test = IsThereThresholdCategory,
                                yes = paste0("-", paste0(WhichToPlot, collapse = ",")),
                                no = "") # Hi/Med/Lo Threshold info
                         )
dir.create(output_dir_All)
setwd(paste0(Original_WD, "/", output_dir_All))



#### Generate primary data frame to loop off of to generate heatmaps ###########
# Define the File_Num(s) to set the threshold against
Threshold_FileNum <- metadata_exprs_mat[is.na(metadata_exprs_mat$TimePoint) == FALSE & 
                                          metadata_exprs_mat$TimePoint == Threshold_TimePoint, "File_Num"]

# Take the average of the portion of the File_Num defined by Threshold_TimePoint

# # For testing
# DataBy <- "PseudoBulk" # "PseudoBulk" / "Clusters"
# IsThereThresholdCategory <- TRUE # TRUE / FALSE

# Find the cutoff for each File_Num defined by Threshold_FileNum
if (IsThereThresholdCategory) {
  if (DataBy == "Clusters") {
    Threshold_CutOff_Avg <- exprs_mat_Linear %>% 
      # Filter out any ClustersToExlude and only include File_Nums defined by Threshold_FileNum
      dplyr::filter(!Cluster %in% ClustersToExlude, 
                    File_Num %in% Threshold_FileNum, # Only from defined Threshold_FileNum
                    ThreshCat %in% WhichToPlot) %>% # Only from defined ThreshCat
      # Determine which values in each Cluster and ThreshCat are above Set_Threshold for each replicate (i.e. File_Num)
      group_by(Cluster, File_Num, ThreshCat) %>%
      summarize(across(where(is.numeric), 
                       ~ quantile(.x, probs = Set_Threshold)), 
                .groups = "drop") %>% 
      # Take the average of each column by cluster and ThreshCat
      group_by(Cluster, ThreshCat) %>%
      summarize(across(where(is.numeric), 
                       mean), 
                .groups = "drop") %>%
      # Adjust dataframe for downstream analysis
      mutate(File_Num = 0) %>%
      relocate(File_Num) # Move File_Num to the first column
    
    } else if (DataBy == "PseudoBulk") {
      Threshold_CutOff_Avg <- exprs_mat_Linear %>% 
        # Filter out any ClustersToExlude and only include File_Nums defined by Threshold_FileNum
        dplyr::filter(!Cluster %in% ClustersToExlude, 
                      File_Num %in% Threshold_FileNum, # Only from defined Threshold_FileNum
                      ThreshCat %in% WhichToPlot) %>% # Only from defined ThreshCat
        # Determine which values in each ThreshCat are above Set_Threshold for each replicate (i.e. File_Num)
        group_by(File_Num, ThreshCat) %>%
        summarize(across(where(is.numeric), 
                         ~ quantile(.x, probs = Set_Threshold)), 
                  .groups = "drop") %>% 
        # Take the average of each column by ThreshCat
        group_by(ThreshCat) %>%
        summarize(across(where(is.numeric), 
                         mean), 
                  .groups = "drop") %>%
        # Adjust dataframe for downstream analysis
        mutate(File_Num = 0,
               Cluster = 0) %>%
        relocate(File_Num) # Move File_Num to the first column
      }
  
  } else {
    if (DataBy == "Clusters") {
      Threshold_CutOff_Avg <- exprs_mat_Linear %>% 
        # Filter out any ClustersToExlude and only include File_Nums defined by Threshold_FileNum
        dplyr::filter(!Cluster %in% ClustersToExlude, 
                      File_Num %in% Threshold_FileNum) %>% # Only from defined Threshold_FileNum
        # Determine which values in each Cluster are above Set_Threshold for each replicate (i.e. File_Num)
        group_by(Cluster, File_Num) %>%
        summarize(across(where(is.numeric), 
                         ~ quantile(.x, probs = Set_Threshold)), 
                  .groups = "drop") %>% 
        # Take the average of each column by cluster
        group_by(Cluster) %>%
        summarize(across(where(is.numeric), 
                         mean), 
                  .groups = "drop") %>%
        # Adjust dataframe for downstream analysis
        mutate(File_Num = 0) %>%
        relocate(File_Num) # Move File_Num to the first column
      
      } else if (DataBy == "PseudoBulk") {
        Threshold_CutOff_Avg <- exprs_mat_Linear %>% 
          # Filter out any ClustersToExlude and only include File_Nums
          dplyr::filter(!Cluster %in% ClustersToExlude, 
                        File_Num %in% Threshold_FileNum) %>% # Only from defined Threshold_FileNum
          # Determine which values are above Set_Threshold for each replicate (i.e. File_Num)
          group_by(File_Num) %>%
          summarize(across(where(is.numeric), 
                           ~ quantile(.x, probs = Set_Threshold)), 
                    .groups = "drop") %>% 
          # Take the average of each column
          summarize(across(where(is.numeric), 
                           mean), 
                    .groups = "drop") %>%
          # Adjust dataframe for downstream analysis
          mutate(File_Num = 0,
                 Cluster = 0) %>%
          relocate(File_Num) # Move File_Num to the first column
        }
    }



#### Find events above Threshold_CutOff_Avg in each File_Num across markers ####
if (DataBy == "Clusters") {
  if (IsThereThresholdCategory) {
    exprs_mat_Linear_AboveThreshold <- exprs_mat_Linear %>%
      # Only include ThreshCat defined by WhichToPlot
      dplyr::filter(ThreshCat %in% WhichToPlot) %>%
      # Find each value above the Threshold_CutOff_Avg (defined by ThreshCat) for each Cluster across all markers
      group_by(Cluster, ThreshCat) %>% 
      mutate(across(
        where(is.numeric),
        ~ ifelse(test = .x > as.numeric(Threshold_CutOff_Avg[Threshold_CutOff_Avg$Cluster == unique(Cluster) & 
                                                               Threshold_CutOff_Avg$ThreshCat == unique(ThreshCat), cur_column()]),
                 yes = .x, 
                 no = NA)
        ), 
        .groups = "drop")
    } else {
      exprs_mat_Linear_AboveThreshold <- exprs_mat_Linear %>%
        # Find each value above the Threshold_CutOff_Avg for each Cluster across all markers
        group_by(Cluster) %>% 
        mutate(across(
          where(is.numeric),
          ~ ifelse(test = .x > as.numeric(Threshold_CutOff_Avg[Threshold_CutOff_Avg$Cluster == unique(Cluster), cur_column()]),
                   yes = .x, 
                   no = NA)
          ), 
          .groups = "drop")
      }
  
  } else if (DataBy == "PseudoBulk") {
    if (IsThereThresholdCategory) {
      exprs_mat_Linear_AboveThreshold <- exprs_mat_Linear %>%
        # Only include ThreshCat defined by WhichToPlot
        dplyr::filter(ThreshCat %in% WhichToPlot) %>%
        # Find each value above the Threshold_CutOff_Avg (defined by ThreshCat) across all markers
        group_by(ThreshCat) %>% 
        mutate(across(
          where(is.numeric),
          ~ ifelse(test = .x > as.numeric(Threshold_CutOff_Avg[Threshold_CutOff_Avg$ThreshCat == unique(ThreshCat), cur_column()]),
                   yes = .x, 
                   no = NA)
          ), 
          .groups = "drop")
      } else {
        exprs_mat_Linear_AboveThreshold <- exprs_mat_Linear %>%
          # Mutate Cluster column to 0 to denote PseudoBulk and enable grouping to match Threshold_CutOff_Avg
          mutate(Cluster = 0) %>%
          group_by(Cluster) %>% 
          # Find each value above the Threshold_CutOff_Avg across all markers
          mutate(across(
            where(is.numeric),
            ~ ifelse(test = .x > as.numeric(Threshold_CutOff_Avg[ , cur_column()]),
                     yes = .x, 
                     no = NA)
            ), 
            .groups = "drop")
        }
    }



#### Calculate % of events above above Threshold_CutOff_Avg ####################
## Calculate percent above threshold for each File_Num and Cluster
if (DataBy == "Clusters") {
  if (IsThereThresholdCategory) {
    PercAboveThreshold <- exprs_mat_Linear_AboveThreshold %>% 
      # Find the percent of values above the Threshold_CutOff_Avg (defined by ThreshCat) for each File_Num across all markers
      # This will drop the Cluster column
      group_by(File_Num, Cluster, ThreshCat) %>% 
      summarize(across(colnames(exprs_mat_Linear)[!colnames(exprs_mat_Linear) %in% AppendedFileNames],
                       ~ mean(!is.na(.x),         # !is.an(.x) returns all values above threshold as TRUE
                              na.rm = TRUE) * 100 # mean of True/False gives the average (muliply by 100 for percent)
                       ), 
                .groups = "drop")
    } else {
      PercAboveThreshold <- exprs_mat_Linear_AboveThreshold %>% 
        # Find the percent of values above the Threshold_CutOff_Avg (defined by ThreshCat) for each File_Num across all markers
        # This will drop the Cluster column
        group_by(File_Num, Cluster) %>% 
        summarize(across(colnames(exprs_mat_Linear)[!colnames(exprs_mat_Linear) %in% AppendedFileNames],
                         ~ mean(!is.na(.x),         # !is.an(.x) returns all values above threshold as TRUE
                                na.rm = TRUE) * 100 # mean of True/False gives the average (muliply by 100 for percent)
                         ), 
                  .groups = "drop")
      }
  
  } else if (DataBy == "PseudoBulk") {
    if (IsThereThresholdCategory) {
      PercAboveThreshold <- exprs_mat_Linear_AboveThreshold %>% 
        # Find the percent of values above the Threshold_CutOff_Avg (defined by ThreshCat) for each File_Num across all markers
        # This will drop the Cluster column
        group_by(File_Num, ThreshCat) %>% 
        summarize(across(colnames(exprs_mat_Linear)[!colnames(exprs_mat_Linear) %in% AppendedFileNames],
                         ~ mean(!is.na(.x),         # !is.an(.x) returns all values above threshold as TRUE
                                na.rm = TRUE) * 100 # mean of True/False gives the average (muliply by 100 for percent)
                         ), 
                  .groups = "drop") %>% 
        # Add back Cluster column with 0 values to denote PseudoBulk analysis
        mutate(Cluster = 0, 
               .after = File_Num)
      } else {
        PercAboveThreshold <- exprs_mat_Linear_AboveThreshold %>% 
          # Find the percent of values above the Threshold_CutOff_Avg (defined by ThreshCat) for each File_Num across all markers
          # This will drop the Cluster column
          group_by(File_Num) %>% 
          summarize(across(colnames(exprs_mat_Linear)[!colnames(exprs_mat_Linear) %in% AppendedFileNames],
                           ~ mean(!is.na(.x),         # !is.an(.x) returns all values above threshold as TRUE
                                  na.rm = TRUE) * 100 # mean of True/False gives the average (muliply by 100 for percent)
          ), 
          .groups = "drop") %>% 
          # Add back Cluster column with 0 values to denote PseudoBulk analysis
          mutate(Cluster = 0, 
                 .after = File_Num)
        }
    }

## Organize dataframe for easier downstream analysis
# Add 'NA' values to any Cluster that is missing a FileNum - this will create a row for each Cluster represented across all File_Nums
PercAboveThreshold_PlusNA <- PercAboveThreshold %>%
  complete(Cluster, File_Num)

# Add identifying information from metadata_exprs_mat
PercAboveThreshold_Indexed <- merge(metadata_exprs_mat[ , c("File_Name", "Condition", "Reordered", "File_Num", "TimePoint", "Replicate")],
                                    PercAboveThreshold_PlusNA,
                                    by = "File_Num")

# Remove any samples that are considered outliers
if (!is.null(FileNamesToExclude)) {
  PercAboveThreshold_Indexed_OutliersRemoved <- PercAboveThreshold_Indexed %>%
    dplyr::filter(!File_Name %in% FileNamesToExclude)
} else {
  PercAboveThreshold_Indexed_OutliersRemoved <- PercAboveThreshold_Indexed
}

# Reorder data frame
PercAboveThreshold_Indexed_OutliersRemoved_REORDERED <- PercAboveThreshold_Indexed_OutliersRemoved %>% 
  arrange(Cluster, Reordered)



#### Plot histogram by marker with respective Threshold_CutOff_Avg value #######
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
                                     if (IsThereThresholdCategory) paste0("-", paste0(WhichToPlot, collapse = ",")))
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
                                       if (IsThereThresholdCategory) paste0("-", paste0(WhichToPlot, collapse = ",")),
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
                       Extra_Columns = c("File_Num", "Cluster", "File_Name", "Condition", "Reordered", "TimePoint", "Replicate"),
                       Include_Inset = FALSE)
}



#### Save cleaned up csv #######################################################
write.csv(x = PercAboveThreshold_Indexed_OutliersRemoved_REORDERED, 
          file = paste0("Expt", paste(WorkingSampleSet, collapse = ","), # Expt ID
                        paste0("_PercAboveThresh_", 
                               Set_Threshold*100, 
                               "thPerc",
                               "Of", sub("_.*", "", Threshold_Condition)), # Transformation info
                        "_", "RAWvalues", # Dataframe ID
                        if (IsThereThresholdCategory) { # IsThereThresholdCategory
                          paste0("-", paste0(WhichToPlot, collapse = ","))
                        },
                        ".csv"))



#### Subset PercAboveThreshold data to exclude user-defined markers ############
# Identify and remove columns based on user-defined Markers_OmitFromPlot
DataToPlot_PercAboveThresh <- PercAboveThreshold_Indexed_OutliersRemoved_REORDERED %>% 
  select(-c("File_Name", "Reordered", "File_Num", all_of(Markers_OmitFromPlot)))

# Define function to average input dataframes
df_Avg_function <- function(df) {
  if (IsThereThresholdCategory) {
    if (DataBy == "PseudoBulk") {
      df %>%
        # Drops Replicate column (averaged out) and Cluster column
        group_by(Condition, TimePoint, ThreshCat) %>%
        summarize(across(all_of(SelectMarkersToPlot), 
                         ~ mean(.x, na.rm = TRUE)), 
                  .groups = "drop"
                  ) %>%
        # Add back Cluster column
        mutate(Cluster = 0, 
               .after = Condition) %>% 
        as.data.frame()
      } else if (DataBy == "Clusters") {
        df %>%
          # Drops Replicate column (averaged out)
          group_by(Condition, Cluster, TimePoint, ThreshCat) %>%
          summarize(across(all_of(SelectMarkersToPlot), 
                           ~ mean(.x, na.rm = TRUE)), 
                    .groups = "drop"
                    ) %>%
          as.data.frame()
        }
    } else {
      if (DataBy == "PseudoBulk") {
        df %>%
          # Drops Replicate column (averaged out) and Cluster column
          group_by(Condition, TimePoint) %>%
          summarize(across(all_of(SelectMarkersToPlot), 
                           ~ mean(.x, na.rm = TRUE)), 
                    .groups = "drop"
                    ) %>%
          # Add back Cluster column
          mutate(Cluster = 0, 
                 .after = Condition) %>%
          as.data.frame()
        } else if (DataBy == "Clusters") {
          df %>%
            # Drops Replicate column (averaged out)
            group_by(Condition, TimePoint, Cluster) %>%
            summarize(across(all_of(SelectMarkersToPlot), 
                             ~ mean(.x, na.rm = TRUE)), 
                      .groups = "drop"
                      ) %>%
            as.data.frame()
          }
      }
  }

# Calculate the average across replicates by Condition and Cluster
DataToPlot_PercAboveThresh_RepsAvg <- df_Avg_function(DataToPlot_PercAboveThresh)



#### Define normalization values based on Threshold_Condition ##################
# Define the timepoints needed for normalization
if (NormalizeAgainst == "Control Condition") {
  PostThreshold_NormTimepoints <- sort(unique(DataToPlot_PercAboveThresh$TimePoint))
} else {
  PostThreshold_NormTimepoints <- Threshold_TimePoint
}

# Define the rows within DataToPlot_PercAboveThresh to designate as normalization rows
if (DataBy == "PseudoBulk") {
  Normalization_ROW_Avg <- DataToPlot_PercAboveThresh %>% 
    dplyr::filter(TimePoint %in% PostThreshold_NormTimepoints) %>% 
    # Drops Replicate and Cluster columns
    group_by(across(all_of(setdiff(colnames(DataToPlot_PercAboveThresh), c(SelectMarkersToPlot, "Cluster", "Replicate"))))) %>% 
    summarize(across(all_of(SelectMarkersToPlot),
                     ~ mean(.x, na.rm = TRUE)), 
              .groups = "drop") %>% 
    mutate(Cluster = 0,
           .after = TimePoint) %>%
    rename(NormTimepoint = TimePoint)
  } else if (DataBy == "Clusters") {
    Normalization_ROW_Avg <- DataToPlot_PercAboveThresh %>% 
      dplyr::filter(TimePoint %in% PostThreshold_NormTimepoints) %>% 
      # Drops Replicate column
      group_by(across(all_of(setdiff(colnames(DataToPlot_PercAboveThresh), c(SelectMarkersToPlot, "Replicate"))))) %>% 
      summarize(across(all_of(SelectMarkersToPlot),
                       ~ mean(.x, na.rm = TRUE)), 
                .groups = "drop") %>%
      rename(NormTimepoint = TimePoint)
    }


##### Normalize PercAboveThresh data relative to the Normalization_Row_Avg #####
## Calculate the fold change against the respective Normalization_Row_Avg values
DataToPlot_FC <- DataToPlot_PercAboveThresh %>%
  # Append a NormTimePoint column to DataToPlot_PercAboveThresh
  mutate(NormTimepoint = if (NormalizeAgainst == "Control Condition") TimePoint else PostThreshold_NormTimepoints) %>% 
  # Append Normalization_ROW_Avg to DataToPlot_FC through Cluster, ThreshCat (if present), and NormTimePoint
  # Then add "_norm" to end of each appended column
  left_join(Normalization_ROW_Avg, 
            by = c(colnames(DataToPlot_PercAboveThresh)[!colnames(DataToPlot_PercAboveThresh) %in% 
                                                                        c("Condition", "TimePoint", "Replicate", SelectMarkersToPlot)],
                   "NormTimepoint" = "NormTimepoint"), 
            suffix = c("", "_norm")) %>%
  # Avoid division by zero
  mutate(across(ends_with("_norm"), 
                ~ ifelse(. == 0, NA, .))) %>% 
  # Calculate the Fold Change by dividing the DataToPlot_FC columns by the "_norm" columns and add "FC_" to front of each appended column
  mutate(across(all_of(SelectMarkersToPlot),
                ~ . / get(paste0(cur_column(), "_norm")),
                .names = "FC_{.col}")) %>%
  # Select only the columns that underwent Fold Change transformation
  select(colnames(DataToPlot_PercAboveThresh)[!colnames(DataToPlot_PercAboveThresh) %in% SelectMarkersToPlot],
         starts_with("FC_")) %>% 
  # Remove "FC_" from the column names
  rename_with(~ sub("^FC_", "", .), starts_with("FC_"))


# Calculate the average fold change across replicates by Condition and Cluster
DataToPlot_FC_RepsAvg <- df_Avg_function(DataToPlot_FC)


## Calculate the log2(Fold Change) of the DataToPlot_FC data
DataToPlot_log2FC <- DataToPlot_FC %>%
  mutate(across(all_of(SelectMarkersToPlot), 
                ~ log2(.x)
                )) %>% # Mutate across non-identifying columns
  mutate(across(all_of(SelectMarkersToPlot), 
                ~ na_if(., Inf) # Remove any Inf (i.e. log2(0))
                )) %>%
  mutate(across(all_of(SelectMarkersToPlot), 
                ~ na_if(., -Inf) # Remove any -Inf (i.e. log2(0))
                ))

# Calculate the average log2(fold change) across replicates by Condition and Cluster
DataToPlot_log2FC_RepsAvg <- df_Avg_function(DataToPlot_log2FC)



#### Save DataToPlot dataframes ################################################
# Put DataToPlot dataframes in a list
DataToPlot_ListOfDatasets <- list(RAWvalues = DataToPlot_PercAboveThresh,
                                  RAWvalues_RepsAvg = DataToPlot_PercAboveThresh_RepsAvg,
                                  FC = DataToPlot_FC,
                                  FC_RepsAvg = DataToPlot_FC_RepsAvg,
                                  log2FC = DataToPlot_log2FC,
                                  log2FC_RepsAvg = DataToPlot_log2FC_RepsAvg)
DataSetNames <- names(DataToPlot_ListOfDatasets)

# Loop through each DataToPlot dataframe and save
for (i in 1:length(DataToPlot_ListOfDatasets)) {
  DataToPlot_FileName <- paste0("Expt", paste(WorkingSampleSet, collapse = ","), # Expt ID
                                              paste0("_PercAboveThresh_", 
                                              Set_Threshold*100, 
                                              "thPerc",
                                              "Of", sub("_.*", "", Threshold_Condition)), # Transformation info
                                       if (grepl("NC", paste(PostThreshold_NormCondition, collapse = ","))) {
                                         "_NormToNC"
                                         } else {
                                           "_NormTo0hr"
                                           }, # Normalization info
                                       "_", names(DataToPlot_ListOfDatasets[i]), # DataToPlot dataframe ID
                                       ".csv")
  
  write.csv(x = DataToPlot_ListOfDatasets[[i]], 
            file = DataToPlot_FileName)
  }



#### Organize and plot data by DotPlot #########################################
## Define function for plotting DotPlots (size = % Abundance, color = normalized values)
Plotting_Function <- function(DotSize_name = DataSetNames[2], DotColor_name, SelectCluster = NULL) {
  # # # For testing Plotting_Function
  # DotSize_name <- DataSetNames[2]
  # DotColor_name <- DataSetNames[4] # 4 (FC) / 6 (log2FC)
  # SelectCluster <- 0 # ClusterOrder / HowManyLoops / 0 (PseudoBulk) / NULL (Clusters, No ThresholdCategory) / 4 (Clusters, ThresholdCategory)
  # i <- "pS6_175Lu" # "pS6_175Lu" / "pATM_159Tb"
  
  InputDataframe_DotSize <- DataToPlot_ListOfDatasets[[DotSize_name]]
  InputDataframe_DotColor <- DataToPlot_ListOfDatasets[[DotColor_name]]
  Treatments <- TreatmentToAnalyze # unique(gsub("_.*", "", InputDataframe_DotSize[["Condition"]]))
  if (is.null(SelectCluster)) {
    Clusters <- unique(InputDataframe_DotSize[["Cluster"]])
  } else {
    Clusters <- SelectCluster
  }
  Marker <- colnames(InputDataframe_DotSize)[colnames(InputDataframe_DotSize) %in% SelectMarkersToPlot] # "pPLK1_113In"
  
  # Combine InputDataframe_DotSize and InputDataframe_DotColor_norm data frames
  DataToAnalyze <- rbind(InputDataframe_DotSize %>% mutate(Dot = "Size", .before = everything()),
                         InputDataframe_DotColor %>% mutate(Dot = "Color", .before = everything())
                         )
  
  # Filter by Treatments and Cluster, then change clusters/timepoints/treatment to factor
  DataToAnalyze_adjusted <- DataToAnalyze[gsub("_.*", "", DataToAnalyze[["Condition"]]) %in% Treatments, ] %>% 
    dplyr::filter(!Cluster %in% ClustersToRemoveFromPlot & 
                    Cluster %in% Clusters) %>% 
    mutate(Cluster = factor(Cluster,
                            levels = if (DataBy == "PseudoBulk") { 0 } else { ClusterOrder }),
           TimePoint = factor(as.character(TimePoint), 
                               levels = as.character(sort(unique(TimePoint)))),
           Treatment = factor(gsub("_.*", "", .data[["Condition"]]),
                              levels = TreatmentToAnalyze),
           .before = everything()
           )
  
  # Define the MaxSize for dots - will be used to standardize dot size if selected
  if (StandardizeDotSize == "AcrossDataset") {
    MaxSize <- InputDataframe_DotSize[which(InputDataframe_DotSize[["Dot"]] == "Size"), ] %>% 
      select(-c("Treatment", "Dot", "Condition", "Cluster", "TimePoint", names(InputDataframe_DotSize)[grepl("ThreshCat", names(InputDataframe_DotSize))])) %>% 
      abs() %>% 
      max(na.rm = TRUE)
  } else if (StandardizeDotSize == "AcrossPlotting") {
    MaxSize <- DataToAnalyze_adjusted[which(DataToAnalyze_adjusted[["Dot"]] == "Size"), ] %>% 
      select(-c("Treatment", "Dot", "Condition", "Cluster", "TimePoint", names(DataToAnalyze_adjusted)[grepl("ThreshCat", names(DataToAnalyze_adjusted))])) %>% 
      abs() %>% 
      max(na.rm = TRUE)
  }
  
  # Define the MaxColor for dots - will be used to standardize dot size if selected
  if (StandardizeDotSize == "AcrossDataset") {
    MaxColor <- InputDataframe_DotColor[which(InputDataframe_DotColor[["Dot"]] == "Color"), ] %>% 
      select(-c("Treatment", "Dot", "Condition", "Cluster", "TimePoint", names(InputDataframe_DotColor)[grepl("ThreshCat", names(InputDataframe_DotColor))])) %>% 
      abs() %>% 
      max(na.rm = TRUE)
  } else if (StandardizeDotSize == "AcrossPlotting") {
    MaxColor <- DataToAnalyze_adjusted[which(DataToAnalyze_adjusted[["Dot"]] == "Color"), ] %>% 
      select(-c("Treatment", "Dot", "Condition", "Cluster", "TimePoint", names(DataToAnalyze_adjusted)[grepl("ThreshCat", names(DataToAnalyze_adjusted))])) %>% 
      abs() %>% 
      max(na.rm = TRUE)
    }
  
  
  ## Loop through clusters and plot each marker over time
  # # Define list into which plots are saved
  # Plot_List <- list()
  
  # Define what will be plotted
  if (DataBy == "PseudoBulk" & !isFALSE(StandardizeDotColor) & !IsThereThresholdCategory) {
    WhatToPlot <- "PseudoBulk" 
    # x <- WhatToPlot
  } else {
    WhatToPlot <- Marker
    # x <- WhatToPlot
  }

  
  Plot_List <- sapply(WhatToPlot, function (x) {
    
    ## Adjust input value as necessary
    if (x == "PseudoBulk") {
      DataToAnalyze_CurrentCluster_ToPlot <- DataToAnalyze_adjusted %>% 
        select(all_of(c(colnames(DataToAnalyze_adjusted)[!colnames(DataToAnalyze_adjusted) %in% SelectMarkersToPlot])), 
               all_of(Marker)) %>% 
        pivot_longer(!c(Treatment, Dot, Condition, Cluster, TimePoint, ThreshCat),
                     names_to = "Marker",
                     values_to = "Value") %>% 
        pivot_wider(names_from = Dot, 
                    values_from = "Value")
    } else {
      DataToAnalyze_CurrentCluster_ToPlot <- DataToAnalyze_adjusted %>% 
        select(all_of(c(colnames(DataToAnalyze_adjusted)[!colnames(DataToAnalyze_adjusted) %in% SelectMarkersToPlot])), 
               all_of(x)) %>% 
        pivot_wider(names_from = Dot, 
                    values_from = x)
      }
    
    ## Add min/max color and size values if StandardizeDotColor == FALSE or StandardizeDotSize == FALSE
    if (isFALSE(StandardizeDotColor)) {
      MaxColor <- max(DataToAnalyze_CurrentCluster_ToPlot[["Color"]], na.rm = TRUE)
    }
    
    if (isFALSE(StandardizeDotSize)) {
      MaxSize <- max(DataToAnalyze_CurrentCluster_ToPlot[["Size"]], na.rm = TRUE)
    }
    
    ## Define color variables
    if (grepl("log2FC", DotColor_name)) {
      # Set the minimum value
      MinColor <- -1 * MaxColor
      
      # In the case that the max/min values are outside of the +/-log2(1.5) range
      if (MaxColor < log2(1.5)) { # max(abs(DataToAnalyze_CurrentCluster_ToPlot$Color), na.rm = TRUE)
        Color_Gradient_colors <- c("white", # viridis(3, option = ColorToUse)[1]
                                   "white") # viridis(3, option = ColorToUse)[2]
        Color_Gradient_values <- c(0,
                                   1)
        
      } else {
        Color_Gradient_colors <- c("blue", # viridis(3, option = ColorToUse)[1]
                                   "white", "white",
                                   "red") # viridis(3, option = ColorToUse)[2]
        Color_Gradient_values <- c(0, 
                                   0.5 - (log2(1.5) / MaxColor), # max(abs(DataToAnalyze_CurrentCluster_ToPlot$Color), na.rm = TRUE)
                                   0.5 + (log2(1.5) / MaxColor), # max(abs(DataToAnalyze_CurrentCluster_ToPlot$Color), na.rm = TRUE)
                                   1)
      }
      
    } else if (grepl("FC", DotColor_name)) {
      # Set the minimum value
      MinColor <- 0
      
      if (MaxColor > 1.5) { # max(DataToAnalyze_CurrentCluster$Color, na.rm = TRUE)
        Color_Gradient_colors <- c("blue", # viridis(3, option = ColorToUse)[1]
                                   "white", "white", 
                                   "red") # viridis(3, option = ColorToUse)[2]
        Color_Gradient_values <- c(0, 
                                   1/1.5, 1.5, 
                                   MaxColor) / MaxColor # max(DataToAnalyze_CurrentCluster_ToPlot$Color, na.rm = TRUE)
      } else {
        Color_Gradient_colors <- c("blue", # viridis(3, option = ColorToUse)[1]
                                   "white", "white")
        Color_Gradient_values <- c(0, 
                                   1/1.5, 1.5) / 1.5
      }
    }
    
    ## Generate and save Plots
    if (IsThereThresholdCategory) {
      ## Plot data by ThreshCat
      Output_Plot <- DataToAnalyze_CurrentCluster_ToPlot %>% 
        ggplot(aes(x = TimePoint, 
                   y = factor(ThreshCat,
                              levels = WhichToPlot)
        )) +
        
        geom_point(aes(size = Size,
                       color = Color)
        ) +
        
        scale_size_continuous(name = "Percent\nAbundance",
                              breaks = if (max(DataToAnalyze_CurrentCluster_ToPlot$Size, na.rm = TRUE) < 10) {
                                seq(from = 0, to = 10, by = 1)
                              } else {
                                seq(from = 0, to = 100, by = 10)
                              },
                              limits = c(0, MaxSize)     # fix the data domain to 0–70 across all plots
        ) +
        
        scale_color_gradientn(limits = c(MinColor, MaxColor),
                              colors = Color_Gradient_colors,
                              values = Color_Gradient_values
        ) +
        
        labs(title = paste0(x),
             # subtitle = paste0("Percent Abundance colored by ",
             #                   ifelse(test = grepl("_FC", DotColor_name),
             #                          yes = "Fold Change",
             #                          no = "log2FC")),
             x = "TimePoint",
             y = paste0("Threshold Category: ", MarkerThresholdedOn, " Level"),
             color =  ifelse(test = grepl("log2FC", DotColor_name),
                             yes = "log2FC\n(Marker/0hr)",
                             no = "Fold Change\n(Marker/0hr)")
        ) +
        
        theme(panel.grid = element_blank(),
              plot.title = element_text(hjust = 0),
              plot.subtitle = element_text(hjust = 0),
              panel.background = element_rect("lightgray") # element_rect("lightgray") / element_blank()
        )
      
    } else if (x == "PseudoBulk") {
      Output_Plot <- DataToAnalyze_CurrentCluster_ToPlot %>% 
        ggplot(aes(x = TimePoint, 
                   y = factor(Marker,
                              levels = rev(unique(Marker)))
        )) +
        
        geom_point(aes(size = Size,
                       color = Color)
        ) +
        
        scale_size_continuous(name = "Percent\nAbundance",
                              breaks = if (max(DataToAnalyze_CurrentCluster_ToPlot$Size, na.rm = TRUE) < 10) {
                                seq(from = 0, to = 10, by = 1)
                              } else {
                                seq(from = 0, to = 100, by = 10)
                              },
                              limits = c(0, MaxSize)
        ) +
        
        scale_color_gradientn(limits = c(MinColor, MaxColor),
                              colors = Color_Gradient_colors,
                              values = Color_Gradient_values
        ) +
        
        labs(title = paste0(x),
             subtitle = paste0("TREATMENTS: ", paste0(TreatmentToAnalyze, collapse = ","), 
                               " / Threshold: ", 
                               (Set_Threshold*100), "% [", Threshold_TimePoint, 
                               "hr] / NORMALIZATION: ", 
                               ColorBy, " against ", NormalizeAgainst
                               ),
             x = "TimePoint",
             y = "Cluster",
             color = "Fold Change\n(Marker/0hr)"
        ) +
        
        theme(panel.grid = element_blank(),
              plot.title = element_text(hjust = 0),
              plot.subtitle = element_text(hjust = 0, size = 6),
              panel.background = element_rect("lightgray") # element_rect("lightgray") / element_blank()
        )
      } else { 
      ## Plot data by Cluster
      Output_Plot <- DataToAnalyze_CurrentCluster_ToPlot %>% 
        ggplot(aes(x = TimePoint, 
                   y = factor(Cluster,
                              levels = rev(ClusterOrder))
        )) +
        
        geom_point(aes(size = Size,
                       color = Color)
        ) +
        
        scale_size_continuous(name = "Percent\nAbundance",
                              breaks = if (max(DataToAnalyze_CurrentCluster_ToPlot$Size, na.rm = TRUE) < 10) {
                                seq(from = 0, to = 10, by = 1)
                              } else {
                                seq(from = 0, to = 100, by = 10)
                              },
                              limits = c(0, MaxSize)
        ) +
        
        scale_color_gradientn(limits = c(MinColor, MaxColor),
                              colors = Color_Gradient_colors,
                              values = Color_Gradient_values
        ) +
        
        labs(title = paste0(x),
             # subtitle = paste0("Percent Abundance colored by ",
             #                   ifelse(test = grepl("_FC", DotColor_name),
             #                          yes = "Fold Change",
             #                          no = "log2FC")),
             x = "TimePoint",
             y = "Cluster",
             color = "Fold Change\n(Marker/0hr)"
        ) +
        
        theme(panel.grid = element_blank(),
              plot.title = element_text(hjust = 0),
              plot.subtitle = element_text(hjust = 0),
              panel.background = element_rect("lightgray") # element_rect("lightgray") / element_blank()
        )
    }
    
    if (ValuesOnPlot) {
      Output_Plot <- Output_Plot +
        geom_text(aes(label = round(Color, digits = 2)), # Add Color values
                  nudge_y = 0.2,
                  size = 3,
                  color = "blue2") +
        geom_text(aes(label = round(Size, digits = 2)), # Add Size values
                  nudge_y = -0.2,
                  size = 3,
                  color = "orange2") + 
        labs(caption = "Color in Blue / Size in Orange") +
        theme(plot.caption = element_text(hjust = 0))
    } else {
      Output_Plot <- Output_Plot +
        # geom_text(aes(label = round(Color, digits = 2)), # Add Color values
        #           nudge_y = 0.2,
        #           size = 3,
        #           color = "blue2") +
        # geom_text(aes(label = round(Size, digits = 2)), # Add Size values
        #           nudge_y = -0.2,
        #           size = 3,
        #           color = "orange2") + 
        labs(caption = paste0("Max Color = ", 
                              round(max(DataToAnalyze_CurrentCluster_ToPlot[["Color"]], na.rm = TRUE), digits = 2),
                              " / Max Size = ",
                              round(max(DataToAnalyze_CurrentCluster_ToPlot[["Size"]], na.rm = TRUE), digits = 2))
             ) +
        theme(plot.caption = element_text(hjust = 0))
    }
    
    # ### Save plot into list
    # Plot_List[[x]] <- Output_Plot
    }, 
    simplify = FALSE)
  
  print("Done with loop through markers")
  
  # Output the plot list
  return(Plot_List) 
  
  }



# # For testing
# Plotting_Function(DotSize_name = DataSetNames[2],
#                   DotColor_name = DataSetNames[4],
#                   SelectCluster = NULL) # if plotting by ThreshCat, select specfic cluster



## Generate list of plots
if (ColorBy == "FoldChange") {
  SizeBy_Value <- 2
  ColorBy_value <- 4
} else if (ColorBy == "log2FC") {
  SizeBy_Value <- 2
  ColorBy_value <- 6
}

## Define how many times to loop through function to save plots
if (DataBy == "PseudoBulk") {
    HowManyLoops <- 0
    } else if (IsThereThresholdCategory) {
      HowManyLoops <- ClusterOrder
      } else if (DataBy == "Clusters") {
          HowManyLoops <- "Clusters"
          } else {
            HowManyLoops <- NULL
            }
  


if(PlotAsGrid) { # Generate grid, if selected, then save
  lapply(HowManyLoops, function(clusterID) {
    if (clusterID == "Clusters") {
      clusterID <- ClusterOrder
      }
    
    ## Run Plotting_Function
    OutputPlots_List <- Plotting_Function(DotSize_name = DataSetNames[SizeBy_Value],
                                          DotColor_name = DataSetNames[ColorBy_value],
                                          SelectCluster = clusterID)
    
    
    height_scalar <- if (IsThereThresholdCategory) {
      length(unique(DataToPlot_ListOfDatasets[[DataSetNames[SizeBy_Value]]][["ThreshCat"]]))
    } else if (DataBy == "PseudoBulk") {
      length(intersect(colnames(DataToPlot_ListOfDatasets[[DataSetNames[SizeBy_Value]]]), 
                       Markers_Plot))
    } else {
      length(unique(DataToPlot_ListOfDatasets[[DataSetNames[SizeBy_Value]]][["Cluster"]]))
    }
    
    
    ## Arrange plots
    # Define plot dimensions for the outputted plot
    plot_height <- height_scalar * (length(OutputPlots_List) / plot_column_length) 
    plot_width <- SizePerPlot * plot_column_length
    
    # Pad the plot list if needed (ex. if number of plots is not divisible by plot_column_length)
    remainder <- length(OutputPlots_List) %% plot_column_length
    if (remainder != 0) {
      NumberOfGridSlots <- ceiling(length(OutputPlots_List) / plot_column_length) * plot_column_length
      for (k in (length(OutputPlots_List) + 1):NumberOfGridSlots) {
        OutputPlots_List[[k]] <- ggplot() + theme_void()
        }
      }
    
    
    # Arrange and add title
    PlotsInGrid <- grid.arrange(grobs = OutputPlots_List[match(MarkerOrder[MarkerOrder %in% names(OutputPlots_List)], 
                                                               names(OutputPlots_List))],
                                nrow = ceiling(length(OutputPlots_List) / plot_column_length),
                                ncol = plot_column_length)
    
    grid_plot <- grid.arrange(textGrob(paste0("TREATMENTS: ", paste0(TreatmentToAnalyze, collapse = ","), 
                                              " / Threshold: ", 
                                              (Set_Threshold*100), "% [", Threshold_TimePoint, 
                                              "hr] / NORMALIZATION: ", 
                                              ColorBy, " against ", NormalizeAgainst, 
                                              " / BULK or CLUSTER: ",
                                              if (DataBy == "PseudoBulk") {
                                                "PseudoBulk"
                                              } else if (DataBy == "Clusters" & !IsThereThresholdCategory) {
                                                ifelse(test = length(ClustersToRemoveFromPlot) > 1, 
                                                       yes = paste("Clusters", paste0(ClusterOrder, collapse = ", ")),
                                                       no = "AllClusters")
                                              } else if (DataBy == "Clusters" & IsThereThresholdCategory) {
                                                paste("Cluster", clusterID)
                                              }), 
                                       gp = gpar(fontsize = 18, fontface = "bold")), 
                              PlotsInGrid, 
                              ncol = 1, heights = c(0.05, 0.9))
    
    ## Save arranged plots
    ggsave(filename = file.path(Original_WD,
                                output_dir_All,
                                paste0("DotPlot_", # Plot info
                                       paste0(TreatmentToAnalyze, collapse = ","), # Treatment info
                                       "_SizeByFracAbund_", # Size info
                                       paste0("ColorBy", ColorBy, "_"), # Color info
                                       ifelse(test = sum(names(OutputPlots_List) %in% SelectMarkersToPlot) == sum(!Markers_Plot %in% Markers_Analysis),
                                              yes = "AllMarkers",
                                              no = paste0(sum(names(OutputPlots_List) %in% SelectMarkersToPlot), "SelectMarkers")),
                                       if (DataBy == "PseudoBulk") {
                                         "PseudoBulk"
                                       } else if (DataBy == "Clusters" & !IsThereThresholdCategory) {
                                         ifelse(test = length(ClustersToRemoveFromPlot) > 1, 
                                                yes = paste("Clusters", paste0(ClusterOrder, collapse = ", ")),
                                                no = "_AllClusters")
                                       } else if (DataBy == "Clusters" & IsThereThresholdCategory) {
                                         paste0("_Cluster", clusterID)
                                       }, # Clusters included
                                       ifelse(test = StandardizeDotColor != FALSE,
                                              yes = paste0("StdColor", StandardizeDotSize),
                                              no = ""), 
                                       ifelse(test = StandardizeDotSize != FALSE,
                                              yes = paste0("_", "StdSize", StandardizeDotSize),
                                              no = ""),
                                       "-Grid", # Indicate grid
                                       ".", OUTPUT.TYPE) # file type info
                                ),
           plot = grid_plot,
           width = plot_width, height = plot_height,
           units = "in", dpi = 300
           )
    })
} else {
  Plot_Height <- round(((Plot_Height / length(c(ClusterOrder, ClustersToRemoveFromPlot))) * length(ClusterOrder) * 1.2), 
                       digits = 1)
  if (DataBy == "PseudoBulk" & !isFALSE(StandardizeDotColor)) {
    Plot_Width <- Plot_Width * 1.5
    Plot_Height <- Plot_Width * 2
  }
  
  lapply(HowManyLoops, function(clusterID) {
    
    ## Run Plotting_Function
    OutputPlots_List <- Plotting_Function(DotSize_name = DataSetNames[SizeBy_Value],
                                          DotColor_name = DataSetNames[ColorBy_value],
                                          SelectCluster = clusterID)
    
    ## Loop through each plot and save in output folder
    for (j in names(OutputPlots_List)[nchar(names(OutputPlots_List)) > 0]) {
      ggsave(filename = file.path(Original_WD,
                                  output_dir_All,
                                  paste0("DotPlot_", # Plot info
                                         paste0(TreatmentToAnalyze, collapse = ","), # Treatment info
                                         "_SizeByFracAbund_", # Size info
                                         paste0("ColorBy", ColorBy, "_"), # Color info
                                         ifelse(test = length(ClustersToRemoveFromPlot) > 1, 
                                                yes = paste0("Clusters", paste0(ClusterOrder, collapse = ","), "_"),
                                                no = "AllClusters_"), # Clusters included
                                         ifelse(test = StandardizeDotColor != FALSE,
                                                yes = paste0("StdColor", StandardizeDotSize, "_"),
                                                no = ""), 
                                         ifelse(test = StandardizeDotSize != FALSE,
                                                yes = paste0("StdSize", StandardizeDotSize, "_"),
                                                no = ""),
                                         j, # Marker info
                                         ".", OUTPUT.TYPE) # file type info
                                  ),
      plot = OutputPlots_List[[j]],
      width = Plot_Width, # Width in inches
      height = Plot_Height, # Height in inches
      # round(((Plot_Height / length(c(ClusterOrder, ClustersToRemoveFromPlot))) * length(ClusterOrder) * 1.2), digits = 1)
      
      dpi = 300 # Resolution (use 300 for print quality)
      )
    }
    
    
    
  })
}



#### Return to original directory once everything is completed ################
setwd(Original_WD)


print("Complete! Have fun staring at your new plot!")

