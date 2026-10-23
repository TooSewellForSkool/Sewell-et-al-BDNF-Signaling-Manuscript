
# Title: LineGraphAnalysis
# Version: 10.1
# Author: Jonathon Sewell
# Date: 2026-02-25

# Update:
# (1) Updated to be used with 12.1.1 script

# Purpose: Select conditions and markers to plot individually

library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(cowplot)
library(viridis)
library(gridExtra)
library(grid)

Original_WD <- getwd()

### User-defined Variables ====================================================
## Define marker-related variables
# Compare selected markers within a given cluster
PlotBy_marker <- FALSE # TRUE / FALSE - compare select markers across 

# Plot All markers or only selected individual markers
MarkerSelect <- "Individual" # "Individual" / "All"

if (PlotBy_marker) {
  MarkerSelect <- "Individual"
  }

# Select individual marker(s) to be analyzed, if MarkerSelect = "Individual"
IndividualMarkers <- c(# "TrkB_s_158Gd",
                       # "TrkB_i_174Yb",
                       # "pTrkB_170Er",
                       # "pPLCG-2_144Nd",
                       "pERK_171Yb",
                       "pAKt_168Er",
                       # "pSrc_154Sm",
                       "pCreb_176Lu",
                       # "p-p38_142Nd",
                       "pS6_175Lu",
                       # "pSTAT3_153Eu",
                       # "pSTAT5_147Sm",
                       "pGSK3b_139La",
                       # "pPLK1_113In",
                       # "H3K9ac_115In",
                       # "pNFKB_166Er",
                       # "p75NTR_161Dy",
                       "p-cJUN_162Dy"#,
                       # "p-p53_151Eu",
                       # "pATM_159Tb",
                       # "yH2AX_165Ho",
                       # "CC3_173Yb"
                       )

# Select whether there is hi/med/lo information
IsThereThresholdCategory <- FALSE # TRUE / FALSE
# # hi / med / lo / 
# # "p75NTR-hi,TrkB_i-lo" / "p75NTR-hi,TrkB_i-med"/
# # "p75NTR-med,TrkB_i-lo" / "p75NTR-med,TrkB_i-med" / "p75NTR-med,TrkB_i-hi"/ 
# # "p75NTR-lo,TrkB_i-med" / "p75NTR-lo,TrkB_i-hi"


## Define cluster-related variables
# Select the cluster(s) to be analyzed
AnalyzeAllClusters <- FALSE # TRUE / FALSE - TRUE if all clusters are being analyzed; FALSE if only ClustersToAnalyze are analyzed, or for PseudoBulk analysis
ClustersToAnalyze <- c(4) # c() if using PseudoBulk analysis or if plotting by "Treatment"
    # Progenitor: 1, 15
    # Unknown: 5, 10, 11, 14, 19
    # Non-Neuronal: 3, 9, 17, 18, 20
    # Neuronal: 2, 4, 6, 7, 8, 12, 13
    # Neuron lineage: 1, 7, 2, 13, 4

if (AnalyzeAllClusters) { # Cluster component for naming purposes
  ClusterName <- "AnalyzeAllClusters"
  } else if (!AnalyzeAllClusters & length(ClustersToAnalyze) == 0) {
    ClusterName <- "PseudoBulk"
    } else if (!AnalyzeAllClusters & length(ClustersToAnalyze) > 0) {
      ClusterName <- paste(ClustersToAnalyze, collapse = ",")
      }


## Select the treatment(s) to be analyzed 
TreatmentToAnalyze <- c("NoTreat", "BDNF", "CompMed", "ContDepr")
# "NoTreat" / "ContDepr" / "CompMed"  / "BDNF" / "NoTreatTest" / "ContDeprTest" / "BDNFTest" / "K252aTest"
TreatmentName <- paste(TreatmentToAnalyze, collapse = ",")


## Define how to plot data
if (PlotBy_marker) { # automatically selects "Marker" if PlotBy_marker <- TRUE
  if (length(IndividualMarkers) > 1) {
    PlotBy <- "Marker"
    } else { # Checks whether there are multiple markers to compare
      stop("Select more IndividualMarkers to compare")
      }
  } else if(IsThereThresholdCategory) { # automatically selects "ThresholdCategory" if there is hi/med/lo information
    PlotBy <- "ThresholdCategory"
    } else if (length(ClustersToAnalyze) > 0 & length(TreatmentToAnalyze) > 2) { 
      PlotBy <- "Treatment"
      } else if (length(ClustersToAnalyze) == 0) { # automatically selects "Treatment" if no clusters are selected
        PlotBy <- "Treatment"
        } else if (length(ClustersToAnalyze) > 0) {
          PlotBy <- "Cluster"
          }


## Define how to organize the time variable
# Define if an inset is necessary to visualize early timepoints (Categorical = no inset / Continuous = add inset)
TimePoint_VariableType <- "Categorical" # "Categorical" / "Continuous"

# Select the upper time point that defines an "early time point" for an inset plot
EarlyTimePoint_Cutoff <- 5 # hrs


## Define the outputted plot info
PlotAsGrid <- TRUE # TRUE / FALSE
plot_column_length <- 1 # Scalar for histogram plotting

if (PlotBy == "Cluster") {
  ColorToUse <- "Rainbow" # Rainbow to color as in ID UMAP
  } else {
    ColorToUse <- "D" # With viridis, A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo
    }

Plot_Width <- 7
Plot_Height <- 6
Output_File_Type <- "pdf" # "pdf", "jpeg", "tiff", "png", "bmp", "svg"


### Read in and prep data for analysis ========================================
input_file_transformation <- "_FC" # "_FC" / "hr_RAWvalues"
input_file_type <- "\\.csv$"


input_filenames <- list.files(Original_WD, pattern = input_file_type) %>% # finds all files in Original_WD of the input_file_type
  (\(x) x[grepl(input_file_transformation, x)])() # pass the list of file names to an anonymous function that filters for input_file_transformation


# Function to detect file type to be read
read_input_file <- function(path) {
  ext <- tools::file_ext(path)
  switch(ext,
         csv = read_csv(path,show_col_types = FALSE),
         tsv = read_tsv(path),
         txt = read_table(path),
         stop("Unsupported file type: ", ext)
  )
}

# Alternative to using the read_input_file function, though it assumes "," delimited and may be slower than read_csv()
# read.table(file.path(Original_WD, input_filenames[i]),
#                                         header = TRUE,
#                                         sep = ",") %>% 
#                                mutate("X" = sub("_.*", "", input_filenames[i])) %>% 
#                                rename(InputFile = X)


# Read in files and combine files of a similar transformation (i.e. _RepsAvg with _RepsAvg)
# This method is flexible if you have multiple "_FC" and "_FC_RepsAvg" files to combine
# If only one is read in, it still works the same
input_dataframe_all <- c()
input_dataframe_all_avg <- c()
for (i in 1:length(input_filenames)) {
  
  if (!grepl("_RepsAvg", input_filenames[i])) {
    input_dataframe_all <- rbind(input_dataframe_all, 
                                 read_input_file(file.path(Original_WD, input_filenames[i])) %>% 
                                   mutate('...1' = sub("_.*", "", input_filenames[i]),
                                          ThresholdCategory = sub(".*FC-([^\\.]+)\\.csv", "\\1", input_filenames[i])) %>% 
                                   rename(InputFile = '...1')) # Read in file and remove first column (indices)
  } else if (grepl("_RepsAvg", input_filenames[i])) {
    input_dataframe_all_avg <- rbind(input_dataframe_all_avg, 
                                     read_input_file(file.path(Original_WD, input_filenames[i])) %>% 
                                       mutate('...1' = sub("_.*", "", input_filenames[i]),
                                              ThresholdCategory = sub(".*RepsAvg-([^\\.]+)\\.csv", "\\1", input_filenames[i])) %>% 
                                       rename(InputFile = '...1')) # Read in file and remove first column (indices)
  }
  
}

if(IsThereThresholdCategory) { # if there is hi/med/lo information, change ThresholdCategory to a factor
  factor_order <- c("hi", 
                    unique(input_dataframe_all$ThresholdCategory)[-which(unique(input_dataframe_all$ThresholdCategory) %in% 
                                                                           c("hi", "lo"))], 
                    "lo")
  
  input_dataframe_all$ThresholdCategory <- factor(input_dataframe_all$ThresholdCategory, 
                                                  levels = factor_order)
  input_dataframe_all_avg$ThresholdCategory <- factor(input_dataframe_all_avg$ThresholdCategory, 
                                                      levels = factor_order)
}


input_dataframe_list <- list("input_dataframe" = input_dataframe_all, 
                             "input_dataframe_avg" = input_dataframe_all_avg)

RemoveOverlappingConditions <- function(DATAFRAME, DATAFRAME_NAME) {
  
  ColumnOfInterest <- "InputFile"
  
  # Step 1: Define the reference group (e.g., "Expt1")
  reference_group <- unique(DATAFRAME[[ColumnOfInterest]])[1]
  
  # Step 2: Get overlapping Condition values between reference and all others
  reference_conditions <- dplyr::filter(DATAFRAME, DATAFRAME[[ColumnOfInterest]] == reference_group) %>% pull(Condition)
  
  # Step 3: Identify all other InputFile groups
  other_groups <- setdiff(unique(DATAFRAME[[ColumnOfInterest]]), reference_group)
  
  # Step 4: Filter out rows from other groups with overlapping Condition
  filtered_dataframe <- dplyr::filter(DATAFRAME, !(DATAFRAME[[ColumnOfInterest]] %in% other_groups &
                                                     DATAFRAME[["Condition"]] %in% reference_conditions))
  
  return(filtered_dataframe)
}

# Rename the dataframes in input_dataframe_list, remove any overlapping conditions and add "Treatment" column
for (j in 1:length(input_dataframe_list)) {
  assign(names(input_dataframe_list[j]), 
         RemoveOverlappingConditions(input_dataframe_list[[j]], 
                                     names(input_dataframe_list[j])) %>% 
           mutate(Treatment = sub("_[^_]+$", "", .data$Condition), .after = "Condition"))
}


if (MarkerSelect == "All") {
  SelectMarker <- colnames(input_dataframe)[!colnames(input_dataframe) %in% c("InputFile", "Condition", "Treatment", "TimePoint", "Replicate", "Cluster", "ThresholdCategory")]
} else if (MarkerSelect == "Individual") {
  SelectMarker <- IndividualMarkers
}


### Make a function to plot line graphs =======================================
# Create output directory
output_dir <- paste0(format(Sys.time(), "%Y-%m-%d"), # date info
                     " Linegraphs_", # plot info
                     ifelse (test = IsThereThresholdCategory & PlotBy == "ThresholdCategory", # plot by ThresholdCategory, or other
                             yes = paste0("By", PlotBy,"-CompareCat_"), 
                             no = paste0("By", PlotBy, "_")),
                     paste0(TreatmentName, "_"), # Treatments plotted
                     ifelse(test = ClusterName == "AnalyzeAllClusters" | ClusterName == "PseudoBulk",
                            yes = paste0(ClusterName),
                            no = paste0("Cl", ClusterName)) # Clusters included
                     )
output_path <- file.path(Original_WD,  output_dir)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)







Max_Values_ByMarker <- sapply(SelectMarker, function (column) {
  max(input_dataframe[[column]], na.rm = TRUE)
})

Min_Values_ByMarker <- sapply(SelectMarker, function (column) {
  ceiling(sd(input_dataframe[[column]], na.rm = TRUE)) * -1
  })



Plotting_Function <- function(InputDataframe, Treatments, Clusters, Marker) {
  
  # # For testing Plotting_Function
  # InputDataframe <- input_dataframe
  # Treatments <- TreatmentToAnalyze
  # Clusters <- ClustersToAnalyze # 2
  # Marker <- SelectMarker # "pPLK1_113In" 
  
  #### Prep data frame for plotting, given the input parameters ####
  df <- InputDataframe
  
  # Filter by Treatments
  if (length(Treatments) > 0) {
    DataToAnalyze_AllClusters <- df[df$Treatment %in% Treatments, ]
  } else if (length(Treatments) == 0) {
    DataToAnalyze_AllClusters <- df
  }
  
  # Filter by Clusters, then change Clusters to a factor or just character
  if (!AnalyzeAllClusters & length(Clusters) > 0) { # User-selected clusters
    DataToAnalyze <- DataToAnalyze_AllClusters[DataToAnalyze_AllClusters$Cluster %in% Clusters, ] %>% 
      mutate(Cluster = factor(as.character(Cluster),
                              levels = Clusters))
    } else if (AnalyzeAllClusters | length(Clusters) == 0) { # PseudoBulk OR AnalyzeAllClusters
      DataToAnalyze <- DataToAnalyze_AllClusters %>% 
        mutate(Cluster = as.character(Cluster))
      }
  
  # If user-defined, change TimePoint to a factor
  if (TimePoint_VariableType =="Categorical") {
    DataToAnalyze <- DataToAnalyze %>% mutate(TimePoint = factor(as.character(DataToAnalyze$TimePoint), 
                                                                  levels = as.character(unique(DataToAnalyze$TimePoint))))
  } else if (TimePoint_VariableType =="Continuous") {
    DataToAnalyze <- DataToAnalyze
  }
  
  # Change Treatment to a factor
  DataToAnalyze <- DataToAnalyze %>% 
    mutate(Treatment = factor(Treatment,
                              levels = TreatmentToAnalyze))
  
  
  ## Define list into which plots are saved
  Plot_List <- list()
  
  if(PlotBy_marker) { # If comparing markers
    #### Loop through each cluster and save plots into a list ####
    for (c in unique(DataToAnalyze$Cluster)) {
      ### Set up data frames for plotting ###
      ## Perform 0-1 normalization so markers can be compared
      # Set up and apply a function to normalize on 0-1 scale
      range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}
      
      # Normalize data so that max/min are defined using all clusters
      DataToPlot_01 <- DataToAnalyze %>% 
        mutate(across(all_of(colnames(DataToAnalyze)[!colnames(DataToAnalyze) %in% 
                                                       c("InputFile", "Condition", "Treatment", "Cluster", "TimePoint", "Replicate", "ThresholdCategory")]
                             ),
                      ~ range01(.x)
                      ))
      
      
      ## Define the working data frame
      DataToPlot <- DataToPlot_01[ , c("Condition", "Treatment", "TimePoint", "Cluster", "Replicate", all_of(Marker), "ThresholdCategory")] %>% 
        dplyr::filter(Cluster == c) # c
      
      
      ## Organize data to plot markers on same graph
      DataToPlot_pivoted <- DataToPlot %>% 
        pivot_longer(cols = all_of(Marker))
      
      
      ## Determine mean and standard deviation for each timepoint in each treatment
      DataToPlot_pivoted_summary <- DataToPlot_pivoted %>%
        group_by(Condition, Treatment, Cluster, TimePoint, ThresholdCategory, name) %>% 
        summarize(sd = sd(value),
                  value = mean(value),
                  .groups = NULL)
      
      
      ## Copy 0hr entry for Treatments without a 0hr TimePoint - helps with visualization
      # Identify treatments with and without 0-hour entries
      treatments_with_0hr <- DataToPlot_pivoted_summary %>%
        dplyr::filter(TimePoint == 0) %>%
        pull(Treatment) %>%
        unique()
      treatments_without_0hr <- setdiff(unique(DataToPlot_pivoted_summary$Treatment), 
                                        treatments_with_0hr)
      
      # Create a function to generate a 0hr input for treatments_without_0hr (same input as treatments_with_0hr)
      append_0hr <- lapply(treatments_without_0hr, function(trt) {
        DataToPlot_pivoted_summary %>% 
          dplyr::filter(TimePoint == 0) %>% # Create synthetic 0-hour entries
          mutate(Condition = trt, # Change naming to match current treatments_without_0hr
                 Treatment = trt)
      }) %>%
        bind_rows()
      
      # Append copied 0hr entries to original data and change Treatment to factor
      DataToPlot_pivoted_summary <- bind_rows(DataToPlot_pivoted_summary, append_0hr) %>%  
        mutate(Treatment = factor(Treatment,
                                  levels = TreatmentToAnalyze))
      
      
      ### Set up plots ###
      ## Define colors to use
      if (ColorToUse == "Rainbow") {
        numclusters <- length(unique(input_dataframe$Cluster)) # Calculate number of total clusters
        colors <- rainbow(numclusters, alpha = 1) # Generate color palette
        colors <- colors[Clusters]
      } else {
        colors <- viridis(n = length(unique(DataToPlot_pivoted_summary[["name"]])),
                          option = ColorToUse)
      }
      
      
      ## Plot all timepoints
      All_Timepoints <- ggplot(DataToPlot_pivoted, 
             aes(x = TimePoint, y = value, 
                 color = name)) +
        geom_point(size = 1, 
                   alpha = 0.5
                   ) +
        
        # # If you prefer error bars to ribbon for plotting SD:
        # geom_errorbar(data = DataToPlot_pivoted_summary,
        #               aes(x = TimePoint, y = value, 
        #                   ymin = value - sd, ymax = value + sd),
        #               linewidth = 0.2
        #               ) +
        
        # If you prefer ribbon to error bars for plotting SD:
        geom_ribbon(data = DataToPlot_pivoted_summary,
                    aes(x = TimePoint, y = value, 
                        ymin = value - sd, ymax = value + sd, # Ribbon represents sd
                        fill = name, 
                        group = name
                        ),
                    alpha = 0.2,
                    color = NA,
                    inherit.aes = FALSE
                    ) + 
        
        geom_point(data = DataToPlot_pivoted_summary[which(DataToPlot_pivoted_summary$Treatment %in% 
                                                             treatments_without_0hr), ], # Only plot treatments_without_0hr
                   aes(y = value),
                   size = 2
                   ) + 
        
        geom_line(data = DataToPlot_pivoted_summary[which(DataToPlot_pivoted_summary$Treatment %in% 
                                                    treatments_without_0hr), ], # Only plot treatments_without_0hr
                  aes(y = value, group = name),
                  linewidth = 0.5
                  ) +    
        geom_point(data = DataToPlot_pivoted_summary[which(DataToPlot_pivoted_summary$Treatment %in% 
                                                     treatments_with_0hr), ], # Plot treatments_with_0hr as top layer
                   aes(y = value),
                   size = 2
                   ) + 
        
        labs(color = paste0("Signaling Marker\nCluster ", c),
             fill = paste0("Signaling Marker\nCluster ", c),
             # title = if (input_file_transformation == "_FC") {
             #   if (length(Clusters) == 0) {
             #     paste0("Fold Change (", c, ") over time: ",  paste0(ClusterName, collapse = ","))
             #   } else {
             #     paste0("Fold Change (", c, ") over time: Cluster ",  paste0(ClusterName, collapse = ","))
             #   }
             # } else if (input_file_transformation == "_log2FC") {
             #   if (length(Clusters) == 0) {
             #     paste0("log2FC(", c, ") over time ",  paste0(ClusterName, collapse = ","))
             #   } else {
             #     paste0("log2FC(", c, ") over time: Cluster ", paste0(ClusterName, collapse = ","))
             #   }
             # },
             x = "Time (hr)",
             y = if (input_file_transformation == "_FC") {
               paste0("Fold Change (Marker / 0hr) - 0-1 Normalized")
             } else if (input_file_transformation == "_log2FC") {
               paste0("log2FC(Marker / 0hr) - 0-1 Normalized")
             }
             ) +
        
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        
        ylim(0, 1) +
        
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      Output_Plot <- ggdraw() +
        draw_plot(All_Timepoints)
      
      
      ### Save plot into list
      Plot_List[[c]] <- Output_Plot
      
      }
    } else { # If not comparing markers
      #### Loop through each marker in Markers and save plots into a list ####
      for (b in Marker) { # Indicates cluster (HowManyLoops > 1) or ClustersToAnalyze (HowManyLoops = 1)
        ### Set up data frames for plotting ###
        # Define the working data frame
        DataToPlot <- DataToAnalyze[ , c("Condition", "Treatment", "TimePoint", "Cluster", "Replicate", b, "ThresholdCategory")] %>% 
          dplyr::filter(Cluster %in% if (length(ClustersToAnalyze) > 1) { ClustersToAnalyze } else { unique(DataToAnalyze$Cluster) }
                        )
        
        # Calculate summary statistics
        DataToPlot_summary <- DataToPlot %>%
          group_by(Condition, Treatment, Cluster, TimePoint, ThresholdCategory) %>% 
          summarize(Mean = mean(.data[[b]]),
                    SD = sd(.data[[b]]))
        
        ### Copy 0hr entry for Treatments without a 0hr TimePoint - helps with visualization ###
        # Identify treatments with and without 0-hour entries
        treatments_with_0hr <- DataToPlot_summary %>%
          dplyr::filter(TimePoint == 0) %>%
          pull(Treatment) %>%
          unique()
        treatments_without_0hr <- setdiff(unique(DataToPlot_summary$Treatment), 
                                          treatments_with_0hr)
        
        # Create a function to generate a 0hr input for treatments_without_0hr (same input as treatments_with_0hr)
        append_0hr <- lapply(treatments_without_0hr, function(trt) {
          DataToPlot_summary %>% 
            dplyr::filter(TimePoint == 0) %>% # Create synthetic 0-hour entries
            mutate(Condition = trt, # Change naming to match current treatments_without_0hr
                   Treatment = trt)
          }) %>%
          bind_rows()
        
        # Append copied 0hr entries to original data and change Treatment to factor
        DataToPlot_summary <- bind_rows(DataToPlot_summary, append_0hr) %>%  
          mutate(Treatment = factor(Treatment,
                                    levels = TreatmentToAnalyze))
        
        ### Set up plots ###
        # Define colors to use
          if (PlotBy == "Treatment" | 
              PlotBy == "ThresholdCategory") {
            colors <- c("lightblue4", 
                        viridis(n = length(unique(DataToPlot_summary[[PlotBy]])),
                                option = ColorToUse)[-3]
                        )
            } else if (ColorToUse == "Rainbow") {
              numclusters <- length(unique(input_dataframe$Cluster)) # Calculate number of total clusters
              colors <- rainbow(numclusters, alpha = 1) # Generate color palette
              colors <- colors[Clusters]
              } else {
                colors <- viridis(n = length(unique(DataToPlot_summary[[PlotBy]])),
                                  option = ColorToUse)
                }
          
          # Plot all TimePoint
          All_Timepoints <- ggplot(DataToPlot, 
                                   aes(x = TimePoint, y = .data[[b]], 
                                       color = !!sym(PlotBy))) +
            geom_point(size = 1, 
                       alpha = 0.5
            ) +
            
            # # If you prefer error bars to ribbon for plotting SD:
            # geom_errorbar(data = DataToPlot_summary,
            #               aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD),
            #               linewidth = 0.2
            #               ) +
            
            # If you prefer ribbon to error bars for plotting SD:
            geom_ribbon(data = DataToPlot_summary,
                        aes(x = TimePoint, y = Mean, 
                            ymin = Mean - SD, ymax = Mean + SD, # Ribbon represents SD
                            fill = !!sym(PlotBy), 
                            group = !!sym(PlotBy)
                        ),
                        alpha = 0.2,
                        color = NA,
                        inherit.aes = FALSE
                        ) + 
            
            geom_point(data = DataToPlot_summary[which(DataToPlot_summary$Treatment %in% 
                                                         treatments_without_0hr), ], # Only plot treatments_without_0hr
                       aes(y = Mean),
                       size = 2
                       ) + 
            geom_line(data = DataToPlot_summary[which(DataToPlot_summary$Treatment %in% 
                                                        treatments_without_0hr), ], # Only plot treatments_without_0hr
                      aes(y = Mean, group = !!sym(PlotBy)),
                      linewidth = 0.5
                      ) +    
            geom_point(data = DataToPlot_summary[which(DataToPlot_summary$Treatment %in% 
                                                         treatments_with_0hr), ], # Plot treatments_with_0hr as top layer
                       aes(y = Mean),
                       size = 5,
                       alpha = 0.75,
                       shape = 18
                       ) + 
          
            labs(title = b,
                 color = PlotBy,
                 fill = PlotBy,
                 # title = if (input_file_transformation == "_FC") {
                 #   if (length(Clusters) == 0) {
                 #     paste0("Fold Change (", b, ") over time: ",  paste0(ClusterName, collapse = ","))
                 #   } else {
                 #     paste0("Fold Change (", b, ") over time: Cluster ",  paste0(ClusterName, collapse = ","))
                 #   }
                 # } else if (input_file_transformation == "_log2FC") {
                 #   if (length(Clusters) == 0) {
                 #     paste0("log2FC(", b, ") over time ",  paste0(ClusterName, collapse = ","))
                 #   } else {
                 #     paste0("log2FC(", b, ") over time: Cluster ", paste0(ClusterName, collapse = ","))
                 #   }
                 # },
                 x = "Time (hr)",
                 y = if (input_file_transformation == "_FC") {
                   paste0("Fold Change (", b, " / 0hr)")
                   } else if (input_file_transformation == "_log2FC") {
                     paste0("log2FC(", b, " / 0hr)")
                     }
                 ) +
            
            scale_color_manual(values = colors) +
            scale_fill_manual(values = colors) +
            
            ylim(Min_Values_ByMarker[[b]], Max_Values_ByMarker[[b]]+abs(Min_Values_ByMarker[[b]])) +
            
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          
          
          # If TimePoint is a continuous variable, plot early TimePoint
          if (TimePoint_VariableType == "Continuous") {
            Early_Timepoints <- ggplot(DataToPlot[which(DataToPlot[, "TimePoint"] < EarlyTimePoint_Cutoff), ], 
                                       aes(x = TimePoint, y = .data[[b]], color = !!sym(PlotBy))) +
              geom_point(size = 1, 
                         alpha = 0.5) +
              
              # # If you prefer error bars to ribbon for plotting SD:
              # geom_errorbar(data = DataToPlot_summary[which(DataToPlot_summary[, "TimePoint"] < EarlyTimePoint_Cutoff), ],
              #               aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD),
              #               linewidth = 0.2) + 
              
              # # If you prefer ribbon to error bars for plotting SD:
              geom_ribbon(data = DataToPlot_summary[which(DataToPlot_summary[, "TimePoint"] < EarlyTimePoint_Cutoff), ],
                          aes(x = TimePoint, y = Mean, 
                              ymin = Mean - SD, ymax = Mean + SD, # Ribbon represents SD
                              fill = !!sym(PlotBy), 
                              group = !!sym(PlotBy)
                              ),
                          alpha = 0.2,
                          color = NA,
                          inherit.aes = FALSE
                          ) + 
              
              geom_point(data = DataToPlot_summary[which(DataToPlot_summary[, "TimePoint"] < EarlyTimePoint_Cutoff) &
                                                     which(DataToPlot_summary$Treatment %in% treatments_without_0hr), ], # Only plot treatments_without_0hr
                         aes(y = Mean),
                         size = 2
                         ) + 
              geom_line(data = DataToPlot_summary[which(DataToPlot_summary[, "TimePoint"] < EarlyTimePoint_Cutoff) &
                                                    which(DataToPlot_summary$Treatment %in% treatments_without_0hr), ], # Only plot treatments_without_0hr
                        aes(y = Mean, group = !!sym(PlotBy)),
                        linewidth = 0.5
                        ) +
              geom_point(data = DataToPlot_summary[which(DataToPlot_summary[, "TimePoint"] < EarlyTimePoint_Cutoff) &
                                                     which(DataToPlot_summary$Treatment %in% treatments_with_0hr), ], # Plot treatments_with_0hr as top layer
                         aes(y = Mean),
                         size = 2
                         ) + 
              
              labs(# color = PlotBy,
                # title = paste0("", b, " over time"),
                x = NULL, #"Time"
                y = NULL # b
                ) +
              
              scale_color_manual(values = colors) +
              scale_fill_manual(values = colors) +
              
              theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                    legend.position = "none",
                    plot.background = element_rect(color = "black", fill = "white", size = 1)
                    )
            }
          
          if (TimePoint_VariableType == "Categorical") {
            Output_Plot <- ggdraw() +
              draw_plot(All_Timepoints)
            } else if (TimePoint_VariableType == "Continuous") {
              Output_Plot <- ggdraw() +
                draw_plot(All_Timepoints) +
                draw_plot(Early_Timepoints, 
                          x = 0.47, y = 0.54, width = 0.4, height = 0.4) 
                    #     x = 0.07,   y = 0.54, width = 0.4, height = 0.4   for y-axis > 10
                    #     x = 0.1,    y = 0.54, width = 0.4, height = 0.4    for y-axis < 10
                    #     x = 0.1,    y = 0.44, width = 0.5, height = 0.5
                    # Right-shifted: 
                    #     x = 0.375,  y = 0.44, width = 0.5, height = 0.5
                    #     x = 0.47,   y = 0.54, width = 0.4, height = 0.4
              }
          
          ### Save plot into list
          Plot_List[[b]] <- Output_Plot
          }
      }
      print("Done with loop through markers")
        
        
      
      return(Plot_List) # output the plot list
    }
  
# Plotting_Function(InputDataframe = input_dataframe,
#                   Treatments = TreatmentToAnalyze,
#                   Clusters = ClustersToAnalyze, # ClustersToAnalyze
#                   Marker = SelectMarker)






#### Generate and save plots ####
if(PlotAsGrid) { # Generate grid, if selected, then save

  #### Define how many times to apply the Plotting_Function (based on if data is PseudoBulk or clustered)
  if (AnalyzeAllClusters & !PlotBy_marker) { # If all clusters are to be plotted separately
    HowManyLoops <- sort(unique(input_dataframe$Cluster))
    } else { # If clusters are being compared, or if PseudoBulk
      HowManyLoops <- 1
      }
  
  #### Apply Plotting_Function to generate a list of plots, arrange as grid, then save
  lapply(HowManyLoops, function(clusterID) {
    
    if (length(HowManyLoops) == 1) {
      clusterID <- ClustersToAnalyze
      }
    
    ## Generate list of plots
    OutputPlots_List <- Plotting_Function(InputDataframe = input_dataframe,
                                          Treatments = TreatmentToAnalyze,
                                        Clusters = clusterID,
                                        Marker = SelectMarker
                                        )
    
    ## Arrange plots
    # Define plot dimensions for the outputted plot
    plot_height <- 5 * (length(OutputPlots_List) / plot_column_length)
    plot_width <- 5 * plot_column_length
    
    # Pad the plot list if needed (ex. if number of plots is not divisible by plot_column_length)
    remainder <- length(OutputPlots_List) %% plot_column_length
    if (remainder != 0) {
      NumberOfGridSlots <- ceiling(length(OutputPlots_List) / plot_column_length) * plot_column_length
      for (k in (length(OutputPlots_List) + 1):NumberOfGridSlots) {
        OutputPlots_List[[k]] <- ggplot() + theme_void()
      }
    }
    
    # Arrange and add title
    PlotsInGrid <- grid.arrange(grobs = OutputPlots_List,
                              nrow = ceiling(length(OutputPlots_List) / plot_column_length),
                              ncol = plot_column_length)
    
    grid_plot <- grid.arrange(textGrob(ifelse(test = ClusterName != "PseudoBulk" & ClusterName != "AnalyzeAllClusters",
                                              yes = paste0("Cluster", paste0(clusterID, collapse = ","), " Plots"),
                                              no = paste0(ClusterName, " Plots")), 
                          gp = gpar(fontsize = 18, fontface = "bold")), 
                          PlotsInGrid, 
                 ncol = 1, heights = c(0.05, 0.9))
    
    ## Save arranged plots
    ggsave(filename = file.path(output_path,
                                paste0("Linegraph_", # plot info
                                       ifelse(test = ClusterName == "PseudoBulk",
                                              yes = paste0(ClusterName, "_"),
                                              no = "Clusters_"), # What values are depicted?
                                       paste0("By", PlotBy), # PlotBy info
                                       if (IsThereThresholdCategory) paste0("-CompareCat"), # Threshold info
                                       "-Grid", # Indicate grid
                                       ifelse(test = ClusterName != "PseudoBulk" & ClusterName != "AnalyzeAllClusters",
                                              yes = paste0("_Cluster", paste0(clusterID, collapse = ",")),
                                              no = ""), # Clusters included
                                       ".", Output_File_Type) # file type info
                                ),
           plot = grid_plot,
           width = plot_width, height = plot_height,
           units = "in", dpi = 300
           )
    })
  
  } else { # Or save each file
    #### Define how many times to apply the Plotting_Function (based on if data is PseudoBulk or clustered)
    if (AnalyzeAllClusters) { # If all clusters are to be plotted separately
      HowManyLoops <- sort(unique(input_dataframe$Cluster))
      } else { # If clusters are being compared, or if PseudoBulk
        HowManyLoops <- 1
        }
    
    #### Apply Plotting_Function to generate a list of plots, arrange as grid, then save
    lapply(HowManyLoops, function(clusterID) {
      
      if (length(HowManyLoops) == 1) {
        clusterID <- ClustersToAnalyze
        }
      
      ## Generate list of plots
      OutputPlots_List <- Plotting_Function(InputDataframe = input_dataframe,
                                                         Treatments = TreatmentToAnalyze,
                                                         Clusters = clusterID,
                                                         Marker = SelectMarker
                                                         )
      
      
      ## Create output folder for each clusterID and set as working directory
      CLUSTER.FOLDER.NAME <- paste0("Cluster ", paste0(clusterID, collapse = ","))
      CLUSTER.FOLDER <- paste0(output_path, "/", CLUSTER.FOLDER.NAME)
      dir.create(CLUSTER.FOLDER)
      
      setwd(CLUSTER.FOLDER)
    
    
      ## Loop through each plot and save in output folder
      for (i in names(OutputPlots_List)) {
        ggsave(filename = paste0("Linegraph_", # plot info
                             paste0("By", PlotBy, "_"), # PlotBy info
                             paste0("Cluster", paste0(clusterID, collapse = ","), "-"), # Cluster info
                             i, # Marker info
                             ".", Output_File_Type), # file type info
               plot = OutputPlots_List[[i]],
               width = Plot_Width, # Width in inches
               height = Plot_Height, # Height in inches
               dpi = 300 # Resolution (use 300 for print quality)
               )
        }
      
      ## Reset to previous working directory
      setwd(output_path)
      
      })
  }


#### Reset to original directory ####
setwd(Original_WD)

paste("Have fun staring at your plots!!!")