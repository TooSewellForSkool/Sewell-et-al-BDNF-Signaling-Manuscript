
# Author: Jonathon Sewell
# Date: 2026-02-24
# Version: 7
# Update:
# (1) Addressed bugs



### Packages to call (or install) ###
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(cowplot)
library(ggsignif)
library(purrr)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(grid)


Original_WD <- getwd() # setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) to set source file location as working directory



#### Read in user-defined files and combine files of a similar transformation (i.e. _RepsAvg with _RepsAvg)
# Input file information
input_file_transformation <- "_FC" # Don't select "_RepsAvg" - this will prevent plotting of individual points
input_file_type <- "csv"

input_filename <- list.files(Original_WD, pattern = paste0("\\.", input_file_type, "$")) %>% # finds all files in Original_WD of the input_file_type
  (\(x) x[grepl(paste0(input_file_transformation, ".", input_file_type), x)])() # filter for input_file_transformation

# Plot All markers or only selected individual markers
MarkerSelect <- "Individual" # "Individual" / "All"

# Select individual marker(s) to be analyzed, if MarkerSelect = "Individual"
IndividualMarkers <- c(#"TrkB_s_158Gd"#,
                       # "TrkB_i_174Yb",
                       # "pTrkB_170Er",
                       "pPLCG-2_144Nd",
                       # "pERK_171Yb",
                       # "pAKt_168Er",
                       "pSrc_154Sm",
                       # "pCreb_176Lu",
                       "p-p38_142Nd",
                       # "pS6_175Lu",
                       "pSTAT3_153Eu",
                       "pSTAT5_147Sm",
                       # "pGSK3b_139La"#,
                       "pPLK1_113In",
                       "H3K9ac_115In",
                       "pNFKB_166Er",
                       # "p75NTR_161Dy",
                       "p-cJUN_162Dy",
                       "p-p53_151Eu",
                       "pATM_159Tb",
                       "yH2AX_165Ho",
                       "CC3_173Yb"
                       )

# Define the clusters to include
ClustersToAnalyze <- c()  # c(1, 7, 2, 12, 13, 6, 8, 4) / c() if using PseudoBulk data / "All" if using all clusters
# Progenitor: 1, 15
# Unknown: 5, 10, 11, 14, 19
# Non-Neuronal: 3, 9, 17, 18, 20
# Neuronal: 2, 4, 6, 7, 8, 12, 13
# Neuron lineage: 1, 7, 2, 13, 4


if (is.null(ClustersToAnalyze)) {
  ClusterName <- "PseudoBulk"
} else if (unique(ClustersToAnalyze == "All")) { # Cluster component for naming purposes
  ClusterName <- "AnalyzeAllClusters"
} else if (length(ClustersToAnalyze) > 0) {
  ClusterName <- paste(ClustersToAnalyze, collapse = ",")
}

# Define the conditions to compare (adjust names as needed to match the 'Treatment' column of metadata.csv)
CONDITIONS_TO_COMPARE <- c("ContDepr", "BDNFTest", "BDNF+K252aTest")
CONDITIONS_TO_COMPARE_Naming <- factor(x = c("ContDepr", "BDNFTest", "K252aTest"),
                                       levels = c("ContDepr", "BDNFTest", "K252aTest"))

# Define Plotting variables
STEP_INCREASE <- 0.15 # fraction of total height of plot to increase separation between error bars
ColorToUse <- "C" # With viridis, A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo
SameColorsAsUMAP <- FALSE # TRUE / FALSE - color all bars one color that matches the UMAP color of the cluster being plotted

Variable_yAxis_Scale <- TRUE # TRUE = y-axis scaled for each plot / FALSE = y-axis scaled by marker

PlotAsGrid <- TRUE # TRUE / FALSE
plot_column_length <- 5 # Scalar for histogram plotting

Plot_Dimensions <- c(4 , 5) # Width x Height (in.)
OUTPUT.TYPE <- "pdf" # choose png or pdf - for histogram and heatmaps



#### Read in input file and make adjustments ###################################
# Read in input file
input_dataframe <- read_csv(file.path(Original_WD, input_filename), 
                            show_col_types = FALSE) %>% 
  mutate('...1' = sub("_.*", "", input_filename)) %>% # First column now identifies the Expt the input data is from
  rename(InputFile = '...1') # Changes first column name


# Redefine 'Cluster', if clusters are not defined (NA)
if (length(unique(input_dataframe$Cluster)) > 1) {
  input_dataframe <- input_dataframe
  } else { # For Pseudobulk dataframe, this changes value from 'NA' to 0, which cleans things up a bit (maybe)
    input_dataframe <- input_dataframe %>% mutate(Cluster = 0)
  }

# Define unique clusters within the input file


if (is.null(ClustersToAnalyze)) {
  CLUSTERS_TO_PLOT <- 0
} else if (unique(ClustersToAnalyze == "All")) {
  CLUSTERS_TO_PLOT <- unique(input_dataframe$Cluster)
} else if (length(ClustersToAnalyze) > 0) {
  CLUSTERS_TO_PLOT <- ClustersToAnalyze
} else {
  CLUSTERS_TO_PLOT <- NULL
}


# if (length(unique(input_dataframe$Cluster)) > 1) {
#   CLUSTERS_TO_COMPARE <- unique(input_dataframe$Cluster)
# } else {
#   CLUSTERS_TO_COMPARE <- NULL
# }

# Define markers (column titles) within the input file
# MARKERS_TO_COMPARE <- setdiff(colnames(input_dataframe), 
#                               c("InputFile", "Sample_Condition", "TimePoints", "Replicate", "Cluster"))

if (MarkerSelect == "All") {
  MARKERS_TO_COMPARE <- setdiff(colnames(input_dataframe), 
                          c("InputFile", "Sample_Condition", "TimePoints", "Replicate", "Cluster"))
  } else if (MarkerSelect == "Individual") {
    MARKERS_TO_COMPARE <- IndividualMarkers
    }


# Define average of input_dataframe
input_dataframe_avg <- input_dataframe %>% 
  group_by(Sample_Condition, TimePoints) %>% 
  summarize(across(all_of(MARKERS_TO_COMPARE), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

#### Prep data frame to plot and define the plotting function ##################
# Define plotting function
BARPLOT_FUNCTION <- function(MARKER, CLUSTER) {
  # For testing function
  # CLUSTER <- CLUSTERS_TO_PLOT # 1
  # MARKER <- MARKERS_TO_COMPARE # "pERK_171Yb" # MARKERS_TO_COMPARE["TrkB_i_174Yb"]
  # c <- "pERK_171Yb" # 0 # 19
  
  # Define the working dataframe by subsetting on function's input variables
  dataframeToPlot <- input_dataframe %>% 
    dplyr::filter(grepl(paste(CONDITIONS_TO_COMPARE_Naming, collapse = "|"), Sample_Condition),
                  Cluster %in% CLUSTER) %>% 
    select(all_of(c("InputFile", "Sample_Condition", "Cluster", "TimePoints", MARKER))) %>% 
    mutate(Treatment = factor(sub("_.*", "", Sample_Condition), 
                              levels = CONDITIONS_TO_COMPARE),
           .before = "Cluster")
  
  # Define colors to use
  if (SameColorsAsUMAP) {
    unique_cluster <- sort(unique(input_dataframe$Cluster))
    cluster_colors <- scales::hue_pal()(length(unique_cluster))
    colors <- rep(cluster_colors[sort(CLUSTER)], times = length(unique(dataframeToPlot$Treatment)))
    } else {
      colors <- viridis(n = length(unique(dataframeToPlot$Treatment)),
                        option = ColorToUse)
      }
  
  # Define max/min value for MARKER to scale y-axis
  MaxValue_MARKER <- dataframeToPlot %>% 
    dplyr::filter(Cluster %in% CLUSTER) %>% 
    select(all_of(MARKER)) %>% 
    max(na.rm = TRUE)
  
  # MinValue_MARKER <- max(ceiling(sd(dataframeToPlot[colnames(dataframeToPlot) %in% MARKER], 
  #                                   na.rm = TRUE) / 0.05) * 0.05 # Round to nearest 0.05
  # ) * -1 # Make negative to set lowest value
  
  
  ## Define list into which plots are saved
  Plot_List <- list()
  
  
  if (ClusterName == "PseudoBulk") {
    VariableToLoop <- MARKER
    
  } else {
    VariableToLoop <- CLUSTER
    MinValue_MARKER <- max(ceiling(sd(dataframeToPlot[[MARKER]], 
                                      na.rm = TRUE) / 0.05) * 0.05 # Round to nearest 0.05
    ) * -1 # Make negative to set lowest value
  }
  
  
  
  #### Loop through each cluster and save plots into a list ####
  for (c in VariableToLoop) {
    # Define variables based on whether this is PseudoBulk data or clustered
    if (ClusterName == "PseudoBulk") {
      dataframeToPlot_filtered <- dataframeToPlot  %>% 
        select(all_of(c("InputFile", "Sample_Condition", "Treatment", "Cluster", "TimePoints", c)))
      
      MinValue_MARKER <- max(ceiling(sd(dataframeToPlot[[c]],
                                        na.rm = TRUE) / 0.05) * 0.05 # Round to nearest 0.05
      ) * -1 # Make negative to set lowest value
      
      PlotVariable <- c
      
      PlotTitle <- paste0("1hr Treatment: Pseudobulk - ", c)
    
      } else {
        
        dataframeToPlot_filtered <- dataframeToPlot %>% 
          dplyr::filter(Cluster == c) # c
        
        PlotVariable <- MARKER
        
        PlotTitle <- paste0("1hr Treatment: Cluster ", c)
        
        }
    
    ### Generate a bargraph that groups TimePoints by Treatment ###
    BARPLOT <- ggplot(data = dataframeToPlot_filtered,
                      aes(x = Treatment, y = .data[[PlotVariable]])) +
      
      geom_bar(data = dataframeToPlot_filtered %>% 
                 group_by(TimePoints, Treatment) %>% 
                 summarize(StdDev = sd(.data[[PlotVariable]]),
                           !!PlotVariable := mean(.data[[PlotVariable]]),
                           .groups = "drop"),
               aes(fill = factor(Treatment)),
               color = "black",
               alpha = 0.5,
               stat = "identity", 
               position = position_dodge(width = 0.9)
               ) +
      
      geom_point(aes(fill = factor(Treatment)), 
                 color = "black",
                 alpha = 1,
                 position = position_dodge(width = 0.9), 
                 size = 2, 
                 shape = 21, 
                 stroke = 0.5
      ) +
      
      geom_errorbar(data = dataframeToPlot_filtered %>% 
                      group_by(TimePoints, Treatment) %>% 
                      summarize(StdDev = sd(.data[[PlotVariable]]),
                                !!PlotVariable := mean(.data[[PlotVariable]]),
                                .groups = "drop"),
                    aes(ymin = .data[[PlotVariable]] - StdDev,
                        ymax = .data[[PlotVariable]] + StdDev, 
                        group = factor(TimePoints)
                    ),
                    position = position_dodge(width = 0.9),
                    width = 0.2,
                    color = "black"
      ) +
      
      geom_signif(data = dataframeToPlot_filtered %>%
                    select(all_of(PlotVariable), Treatment, TimePoints),
                  color = "black",
                  comparisons = combn(CONDITIONS_TO_COMPARE, 2, simplify = FALSE),
                  test = "t.test",
                  test.args = list(paired = FALSE),
                  y_position = max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) + (0.1 * max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE)),
                  step_increase = if (Variable_yAxis_Scale) {
                    STEP_INCREASE
                    } else if (max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) > (MaxValue_MARKER*0.75)) {
                    STEP_INCREASE
                    } else if (max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) < (MaxValue_MARKER*0.75) &
                               max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) > (MaxValue_MARKER*0.5)) { 
                      STEP_INCREASE * 2
                    } else if (max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) < (MaxValue_MARKER*0.5) &
                               max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) > (MaxValue_MARKER*0.25)) { 
                      STEP_INCREASE * 4
                    } else {
                      STEP_INCREASE * 8
                    }
                    , # STEP_INCREASE
                  map_signif_level = function(p) sprintf("p = %.2g", p)
                  ) +
      
      scale_fill_manual(values = colors) +
      
      labs(title = PlotTitle, 
           x = "Treatment",
           y = paste0(PlotVariable, sub("_", ": ",input_file_transformation), " relative to ",
                      input_dataframe_avg$TimePoints[which(round(input_dataframe_avg[[PlotVariable]], digits = 2) == 1)], "hr"),
           fill = "Treatment"
           )
    
    
    if (Variable_yAxis_Scale) {
      BARPLOT <- BARPLOT + 
        scale_y_continuous(breaks =  if (max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) < 2) {
          seq(0, ceiling(max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE)* 1.5), 
              by = 0.2)
        } else if (max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) < 5) {
          seq(0, ceiling(max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE)* 1.5), 
              by = 0.5)
        } else if (max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) < 10) {
          seq(0, ceiling(max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE)* 1.5), 
              by = 1)
        } else if (max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) < 15) {
          seq(0, ceiling(max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE)* 1.5), 
              by = 2)
        } else if (max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE) < 20) {
          seq(0, ceiling(max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE)* 1.5), 
              by = 5)
        } else {
          seq(0, ceiling(max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE)* 1.5))
        },
        limits = c(0, round(max(dataframeToPlot_filtered[[PlotVariable]], na.rm = TRUE)* 1.5, digits = 1))
        ) +
        
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    } else {
      BARPLOT <- BARPLOT +
        ylim(MinValue_MARKER, 
             MaxValue_MARKER + (abs(MinValue_MARKER) * 3)) + 
        
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    
    Output_Plot <- ggdraw() +
      draw_plot(BARPLOT)
    
    
    ### Save plot into list
    Plot_List[[c]] <- Output_Plot
    
    
    }
  
  
  print("Done with loop through markers")
  
  
  return(Filter(f = Negate(f = is.null), Plot_List[VariableToLoop]) # Remove any NULL elements
         ) # output the plot list
  }



#### Create folder for output files based on the current CONDITIONS_TO_COMPARE ####
output_dir_CONDITION <- paste0("16_", 
                               format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), 
                               "_COMPARE_", paste(CONDITIONS_TO_COMPARE, collapse = ","))
dir.create(output_dir_CONDITION)
setwd(paste0(Original_WD, "/", output_dir_CONDITION))




#### Define how to apply BARPLOT_FUNCTION across PlotVarablbeToLoop
if (ClusterName == "PseudoBulk") {
  PlotVarablbeToLoop <- ClusterName
  
} else {
  PlotVarablbeToLoop <- MARKERS_TO_COMPARE
}


#### Generate and save plots ####
if(PlotAsGrid) { # Generate grid, if selected, then save
  
  #### Apply Plotting_Function to generate a list of plots, arrange as grid, then save
  lapply(PlotVarablbeToLoop, function(SelectMarker) {
    
    if (SelectMarker == "PseudoBulk"){
      SelectMarker <- MARKERS_TO_COMPARE
      }
    
    ## Generate list of plots
    OutputPlots_List <- BARPLOT_FUNCTION(CLUSTER = CLUSTERS_TO_PLOT,
                                         MARKER = SelectMarker
                                         ) 
    
    ## Arrange plots
    # Define plot dimensions for the outputted plot
    plot_height <- 5 * ceiling(length(OutputPlots_List) / plot_column_length)
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
                                              yes = paste0("Cluster", paste0(CLUSTERS_TO_PLOT, collapse = ","), " Plots"),
                                              no = paste0(ClusterName, " Plots")), 
                                       gp = gpar(fontsize = 18, fontface = "bold")), 
                              PlotsInGrid, 
                              ncol = 1, heights = c(0.05, 0.9))
    
    ## Save arranged plots
    ggsave(filename = file.path(paste0(Original_WD, "/", output_dir_CONDITION),
                                paste0("BarPlot_", # plot info
                                       ifelse(test = ClusterName == "PseudoBulk",
                                              yes = paste0(ClusterName, "_"),
                                              no = "Clusters_"), # What values are depicted?
                                       ifelse(test = ClusterName != "PseudoBulk" & ClusterName != "AnalyzeAllClusters",
                                              yes = paste0("_Cluster", paste0(CLUSTERS_TO_PLOT, collapse = ",")),
                                              no = ClusterName), # Clusters included
                                       paste0("_", paste0(sub("(*_*)([^_]*)$", "", SelectMarker), collapse = "-")),
                                       "-Grid", # Indicate grid
                                       ".", OUTPUT.TYPE) # file type info
                                ),
    plot = grid_plot,
    width = plot_width, height = plot_height,
    units = "in", dpi = 300
    )
  })
  
} else { # Or save each file
  if (length(CLUSTERS_TO_PLOT) > 1) {
    for (a in CLUSTERS_TO_PLOT) {
      output_dir_CLUSTER <- paste0("COMPARE_", paste(CONDITIONS_TO_COMPARE, collapse = ","),
                                   "_Cluster", a)
      dir.create(output_dir_CLUSTER)
      setwd(paste0(Original_WD, "/", output_dir_CONDITION, "/", output_dir_CLUSTER))
      
      for (b in 1:length(MARKERS_TO_COMPARE)) {
        BARPLOT_FUNCTION(CLUSTER = a,
                         MARKER = MARKERS_TO_COMPARE[b])
        
        
        
        ggsave(filename = paste0(MARKERS_TO_COMPARE[b],
                                 "_", unique(input_dataframe$InputFile),
                                 "_Barplot",
                                 if (length(CLUSTERS_TO_PLOT) > 1) {
                                   paste0("_Cluster", a)
                                 } else {
                                   "_Pseudobulk"
                                 },
                                 ".", OUTPUT.TYPE),
               plot = BARPLOT,
               width = Plot_Dimensions[1],
               height = Plot_Dimensions[2],
               units = c("in")
               )
        
        
      }
      setwd(paste0(Original_WD, "/", output_dir_CONDITION))
    }
    
  } else {
    for (b in 1:length(MARKERS_TO_COMPARE)) {
      BARPLOT_FUNCTION(CLUSTER = 0,
                       MARKER = MARKERS_TO_COMPARE[b])
      
      ggsave(filename = paste0(MARKERS_TO_COMPARE[b],
                               "_", unique(input_dataframe$InputFile),
                               "_Barplot",
                               "_Pseudobulk",
                               ".", OUTPUT.TYPE),
             plot = BARPLOT,
             width = Plot_Dimensions[1],
             height = Plot_Dimensions[2],
             units = c("in")
             )
      }
    }
  }



#### Reset working directory ###################################################
setwd(Original_WD)

print("Complete! Have fun staring at your new plot!")



