
# Author: Jonathon Sewell
# Date: 2026-03-04
# Version: 9

# Updates:
# (1) Updated line graph to include stats


# Purpose: Choose between different plots (Biaxial, Histogram, Percentage, BoxAndWhisker)



library(ggfortify)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(cowplot)
library(grid)
library(gridExtra)
library(ggridges)


INPUT.FOLDER <- getwd() # setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) to set source file location as working directory


#### User-defined input parameters #############################################
## Select which type(s) of plot(s) to generate
PlotBiaxials <- FALSE # TRUE / FALSE
PlotPercentages <- FALSE # TRUE / FALSE
PlotHistograms <- FALSE # TRUE / FALSE
PlotBoxAndWhisker <- FALSE # TRUE / FALSE
PlotLineGraph <- TRUE # TRUE / FALSE


## Define values related to the selected marker(s) to plot
# Marker(s) to threshold on:
SELECTED.MARKER <- c(x = "TrkB_i_174Yb",
                     y = "TrkB_s_158Gd"
)
# p75NTR_161Dy   / TrkB_i_174Yb
# Sox2_160Gd     / BLBP_172Yb
# Tuj1_89Y       / CXC3CR1_141Pr
# Vimentin_143Nd / TrkB_s_158Gd
# pCreb_176Lu    / pS6_175Lu  / pERK_171Yb


## Define how Threshold Categorizations (ex Lo, Med, Hi) will be defined
# User-defined threshold values:
EXPRESSION.DELINEATION.VALUES <- c(x = list(c(1.5)),
                                   y = list(c(0.15, 1.4, 3.7)) # from QuantileValues <- c(1/3, 2/3)
)
# p75NTR_161Dy = c(0.45, 1.45, 3.85) /   TrkB_i_174Yb = c(0.15, 1.4, 3.7) / TrkB_s_158Gd = c(0.5, 0.9, 2.8) or c(0.7, 3) or 
# Sox2_160Gd = 2.45   /   BLBP_172Yb = 1.95
# Tuj1_89Y = 3.15     /   CXC3CR1_141Pr = 2.35
# Vimentin_143Nd = 0.5
# These values can be defined by looking at expression histograms/violin plots, or by using alternative justification

# Threshold values defined by quantiles:
DelineateByQuantile <- FALSE # TRUE / FALSE
QuantileValues <- c(0.5) # c(1/3, 2/3) / c(0.5)

# Operator(s): ">=" or "<=" to indicate "greater than" or "less than" the corresponding EXPRESSION.DELINEATION.VALUES
THRESHOLD.OPERATOR <- c(x = list(c(">=", ">=", ">=")), # c(">=", ">=") / c(">=", ">=", ">=")
                        y = list(c(">=", ">=", ">=")) # c(">=", ">=") / c(">=", ">=", ">=")
                        )

# Have you justified your EXPRESSION.DELINEATION.VALUES values?
DELINEATION.VALUES.justified <- TRUE # TRUE / FALSE


## Define cluster-related variables
CLUSTER.ROUND <- 1 # Need to change this for the level of subsetting/subclustering!
CLUSTERS.MERGED <- FALSE # TRUE / FALSE

ClusterSelection <- "Individual" # "Individual" / "PseudoBulk" / "All"
ClusterToPlot <- c(1, 2, 4, 7, 13, 15) # c(1, 7, 2, 13, 4)  - automatically adjusts if ClusterSelection <- "PseudoBulk" or "All"


## Define conditions to plot
ConditionsToPlot <- c("NoTreat", "BDNF") 
# c("NoTreat", "ContDepr", "CompMed", "BDNF", "BDNFTest", "BDNF+K252aTest")


## Define variables related to plotting parameters
# Define what to plot
PercentsOnPlots <- FALSE # TRUE / FALSE (for Biaxial and Histogram - automatic for Percentage plots)

MarkerToPlot <- "TrkB_i_174Yb" # NULL / "pAKt_168Er" / "TrkB_s_158Gd" / "Tuj1_89Y" / "CHAT_152Sm"
COLOR_BY_biaxial <- "TimePoint" # SignalingMarker / ThresholdCategory_x / ThresholdCategory_y / ThresholdCategory_xy / TimePoint / Cluster / Cluster2
ColorToUse <- "Default" # viridis color option / "CellID" (color as clusters) / "Default" (use default color palette)
# viridis color options: A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo

PLOT_BY_histogram <- setNames(paste0("ThreshCat_", SELECTED.MARKER[which(!is.na(SELECTED.MARKER))]), 
                              names(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))]))
RidgePlot_By <- "TimePoint"
      # NULL / "TimePoint" / "Treatment" / "Cluster" / "ThreshCat_x" / "ThreshCat_y" / "ThreshCat_xy" / paste0("ThreshCat_", SELECTED.MARKER[["x"]])

FacetBy <- NULL
      # NULL / "TimePoint" / "Treatment" / "Cluster" / "ThreshCat_xy" / paste0("ThreshCat_", SELECTED.MARKER[["x"]])

PLOT_BY_linegraph <- "Cluster"
      # # NULL / "Treatment" / "Cluster" / "ThreshCat_xy" / paste0("ThreshCat_", SELECTED.MARKER[["x"]])
ClusterToPlot_linegraph <- c(1, 2, 4, 7, 13, 15)
Normalize_LineGraph <- "FoldChange" # "FoldChange" / "log2FC" / FALSE
LineGraph_Yaxis_Normalization <- 2 # "Full_Dataset" / "Conditions_Dataset" / "Working_Dataset" / 2 (i.e. hardcode)


# Define output parameters
PlotAsGrid <- TRUE # TRUE / FALSE
if (length(MarkerToPlot) == 1) {
  PlotAsGrid <- FALSE
  }
plot_column_length <- 5 # Scalar for histogram plotting
size_per_plot <- 5

point_size <- 0.1
OUTPUT_HEIGHT <- 4
OUTPUT_WIDTH <- 7
Plot_Dimensions <- c(4 , 5) # Width x Height (in.)

BIAXIAL_OUTPUT_FORMAT <- "png" # "png" or "pdf"
HISTOGRM_OUTPUT_FORMAT <- "pdf" # "png" or "pdf"
PERCENT_OUTPUT_FORMAT <- "pdf" # "png" or "pdf"
BOXandWHISKER_OUTPUT_FORMAT <- "pdf" # "png" or "pdf"
LineGraph_OUTPUT_FORMAT <- "pdf" # "png" or "pdf"

## Make changes based on user-defined variables
if (ClusterSelection != "Individual") {
  ClusterToPlot <- ClusterSelection
}

if (ClusterSelection == "PseudoBulk") {
  PlotAsGrid <- FALSE
}

if (ClusterSelection == "Individual" & is.null(ClusterToPlot_linegraph)) {
  ClusterToPlot_linegraph <- ClusterToPlot
} 

if (COLOR_BY_biaxial == "Cluster") {
  ColorToUse <- "CellID"
  } else if (is.null(PLOT_BY_linegraph) | isFALSE(PlotLineGraph)) {
    ColorToUse <- ColorToUse
    } else if (PlotLineGraph & PLOT_BY_linegraph == "Cluster") {
      ColorToUse <- "CellID"
      } else {
        ColorToUse <- ColorToUse
        }


#### Defined parameters ########################################################
## Input basenames
CLUSTERS.BASENAME <- "/cluster_RX_assigns.csv"
LAYOUT.BASENAME <- "/umap_xy_RX.csv"

## Input filenames
CLUSTERS.FILENAME <- sub("RX", paste0("R", CLUSTER.ROUND), CLUSTERS.BASENAME) #Cluster file name adjusted for round of clustering
if (CLUSTERS.MERGED) { #Adjust cluster file name if clusters were merged
  CLUSTERS.FILENAME <- sub(".csv", "_merged.csv", CLUSTERS.FILENAME)
}
if (CLUSTER.ROUND > 1) {
  CLUSTERS.FILENAME.Primary <- sub("RX", paste0("R", CLUSTER.ROUND-1), CLUSTERS.BASENAME) #Cluster file name adjusted for round of clustering
  if (CLUSTERS.MERGED) { #Adjust cluster file name if clusters were merged
    CLUSTERS.FILENAME.Primary <- sub(".csv", "_merged.csv", CLUSTERS.FILENAME)
  }
}

LAYOUT.FILENAME <- sub("RX", paste0("R", CLUSTER.ROUND), LAYOUT.BASENAME) # UMAP_layout.csv
if (CLUSTER.ROUND > 1) {
  LAYOUT.FILENAME.Primary <- sub("RX", paste0("R", CLUSTER.ROUND-1), LAYOUT.BASENAME)
}

PANEL.FILENAME <- "/panel.csv"
METADATA.FILENAME <- "/metadata.csv"
FILENUMS.FILENAME <- "/filenums.csv"
EXPRS.ANALYSIS.FILENAME <- "/expression_matrix_analysis.csv"
EXPRS.OTHER.FILENAME <- "/expression_matrix_other.csv"

print("Input parameters loaded")



#### Read in necessary files ###################################################
## Read in files
layout.in <- read.csv(paste0(INPUT.FOLDER, LAYOUT.FILENAME))
if (CLUSTER.ROUND > 1) {
  layout.in.primary <- read.csv(paste0(INPUT.FOLDER, LAYOUT.FILENAME.Primary))
}

clusters.in <- read.csv(paste0(INPUT.FOLDER, CLUSTERS.FILENAME))
if (CLUSTER.ROUND > 1) {
  clusters.in.primary <- read.csv(paste0(INPUT.FOLDER, CLUSTERS.FILENAME.Primary))
}
# Adjust ClusterToPlot as needed
if (ClusterSelection != "Individual") {
  ClusterToPlot <- sort(unique(clusters.in[[1]]))
}

panel.in <- read.csv(paste0(INPUT.FOLDER, PANEL.FILENAME))
metadata.in <- read.csv(paste0(INPUT.FOLDER, METADATA.FILENAME))
filenums.in <- read.csv(paste0(INPUT.FOLDER, FILENUMS.FILENAME))
exprs.analysis.in <- fread(paste0(INPUT.FOLDER, EXPRS.ANALYSIS.FILENAME), stringsAsFactors = FALSE, data.table = FALSE)
exprs.other.in <- fread(paste0(INPUT.FOLDER, EXPRS.OTHER.FILENAME), stringsAsFactors = FALSE, data.table = FALSE)

# Define markers to plot, if not already defined by user
if (is.null(MarkerToPlot)) {
  MarkerToPlot <- panel.in$Fixed_Param[panel.in$Plot == 1]
}



print("Necessary files read in")



#### Threshold data based on thresholding parameters ###########################
## Combine exprs_mat data frames to prep thresholding
# thresholding.ct.df <- merge(layout.clusters.id, layout.concat.transform.id)
thresholding.ct.df <- cbind(merge(filenums.in %>% rename(File_Num = "V1"), 
                                  metadata.in[, c("File_Num", "TimePoint", "Treatment", "Replicate")], by = "File_Num"),
                            if (CLUSTER.ROUND > 1) {
                              cbind(clusters.in.primary %>% rename(Cluster = "R1"),
                                    clusters.in %>% rename(Cluster2 = "R2"))
                            } else {
                              clusters.in %>% rename(Cluster = "R1")
                            },
                            exprs.analysis.in, 
                            exprs.other.in)

## Check that you have justified your EXPRESSION.DELINEATION.VALUES - if not, use histogram
if (DELINEATION.VALUES.justified == FALSE) {
  stop("Use histogram (or other method) to justify your selections — If justified, change DELINEATION.VALUES.justified to TRUE")
}

## Modify threshold parameters based on user defined variables and thresholding.ct.df
if (DelineateByQuantile) {
  # Update EXPRESSION.DELINEATION.VALUES to values based on QuantileValues
  EXPRESSION.DELINEATION.VALUES <- sapply(names(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))]), function (n) {
    EXPRESSION.DELINEATION.VALUES[[n]] <- list(quantile(x = thresholding.ct.df[[SELECTED.MARKER[[n]]]],
                                                        probs = QuantileValues) %>% unname())
  })
}


## Generate a ThresholdCategory column for each of SELECTED.MARKER to categorize each event
thresholding.ct.df.category <- cbind(thresholding.ct.df,
                                     sapply(paste0("ThreshCat_", SELECTED.MARKER[which(!is.na(SELECTED.MARKER))]), function(n) {
                                       c(rep(NA, times = nrow(thresholding.ct.df)))
                                       }) %>% as.data.frame()
                                     )


## Apply user-defined thresholding parameters to define each event in the ThresholdCategory column
for (i in names(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))])) {
  # Define the current 'marker'
  marker <- SELECTED.MARKER[[i]]
  
  # Establish values to define ThreshCat values
  if (DelineateByQuantile) {
    value <- c(min(thresholding.ct.df.category[[marker]]), # minimum value
               quantile(x = thresholding.ct.df.category[[marker]],
                        probs = c(1/3, 2/3)) %>% unname(), # delination value(s)
               max(thresholding.ct.df.category[[marker]]) # maximum value
    )
  } else {
    value <- c(min(thresholding.ct.df.category[[marker]]), # minimum value
               as.vector(EXPRESSION.DELINEATION.VALUES[[i]]), # delination value(s)
               max(thresholding.ct.df.category[[marker]]) # maximum value
    )
  }
  # operator <- as.vector(THRESHOLD.OPERATOR[grepl(i, names(THRESHOLD.OPERATOR))])
  
  # Add column to thresholding.ct.df.category to define ThreshCat for each 'marker'
  for (j in seq(length(value)-1)) {
    thresholding.ct.df.category[dplyr::between(x = thresholding.ct.df.category[[marker]], 
                                               left = value[j], # >=
                                               right = value[j+1] # <=
    ), 
    paste0("ThreshCat_", SELECTED.MARKER[[i]])] <- j
  }
}



#### Plot Biaxials ########################################
## Create and set Output folder
# Update plots to generate if necessary
if (length(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))]) != 2) {
  PlotBiaxials <- FALSE
  }

# Name the plots that are to be generated
PlotsToGenerate <- c()
if (PlotBiaxials) {
  PlotsToGenerate <- c(PlotsToGenerate, "Biaxial")
}
if (PlotHistograms) {
  PlotsToGenerate <- c(PlotsToGenerate, "Histogram")
}
if (PlotPercentages) {
  PlotsToGenerate <- c(PlotsToGenerate, "Percentage")
}
if (PlotBoxAndWhisker) {
  PlotsToGenerate <- c(PlotsToGenerate, "BoxAndWhisker")
}
if (PlotLineGraph) {
  PlotsToGenerate <- c(PlotsToGenerate, "LineGraph")
}

## Output folder
ALL.THRESHOLD.NAME <- sapply(names(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))]), function (n) {
  paste0(SELECTED.MARKER[[n]], "-", 
         paste0(ifelse(test = THRESHOLD.OPERATOR[[n]][seq(length(EXPRESSION.DELINEATION.VALUES[[n]]))] == ">=",
                       yes = "Gr",
                       no = "Ls"),
                EXPRESSION.DELINEATION.VALUES[[n]] %>% round(digits = 2),
                collapse = "-")
         )
  }) %>% 
  unname() %>% 
  paste0(collapse = "_")

# Create and set folder to to which plots will be saved
output_dir <- paste0(format(Sys.time(), "%Y-%m-%d"), # Date info
                     " Plot_", 
                     ifelse(test = sum(grepl("LineGraph", PlotsToGenerate)) == 1,
                            yes = paste0(paste0(PlotsToGenerate, collapse = ","), "-", Normalize_LineGraph),
                            no = paste0(PlotsToGenerate, collapse = ",")),
                     ifelse(test = paste0(PlotsToGenerate, collapse = "") == "LineGraph",
                            yes = "",
                            no = paste0("_", ALL.THRESHOLD.NAME)),
                     "_", ClusterSelection)
dir.create(output_dir)

setwd(paste0(INPUT.FOLDER, "/", output_dir))


## Define variables to use in plotting
# Define conditions and their order
Conditions <- unique(thresholding.ct.df.category$Treatment)[match(unique(thresholding.ct.df.category$Treatment),
                                                                  c("NoTreat", "ContDepr", "CompMed", "BDNF",
                                                                    "BDNFTest", "BDNF+K252aTest"))]


## Process dataframe to plot and store COLOR_BY_value in a list
# Filter on ConditionsToPlot and ClusterToPlot, then generate co-exprs column
DataToPlot_filtered <- thresholding.ct.df.category %>% 
  # Filter appropriately
  dplyr::filter(Treatment %in% ConditionsToPlot) %>% 
  # Create column to denote co-expression levels
  unite(col = "ThreshCat_xy", contains("ThreshCat_"), sep = "", remove = FALSE) %>% 
  mutate(ThreshCat_xy = factor(ThreshCat_xy,
                               levels = sort(unique(ThreshCat_xy))))


# If plotting PseudoBulk data, mutate Cluster to 0
if (ClusterSelection == "PseudoBulk") {
  DataToPlot_filtered <- DataToPlot_filtered %>% 
    mutate(Cluster = "PseudoBulk")
}


# Define values that can be used to color plots, which will be selected from user-defined COLOR_BY_biaxial
COLOR_BY_value <- list(SignalingMarker = DataToPlot_filtered[colnames(DataToPlot_filtered) %in% MarkerToPlot],
                       ThresholdCategory_x = factor(as.character(DataToPlot_filtered[[paste0("ThreshCat_", 
                                                                                             SELECTED.MARKER[["x"]])]]), 
                                                    levels = sort(unique(DataToPlot_filtered[[paste0("ThreshCat_", 
                                                                                                     SELECTED.MARKER[["x"]])]]))),
                       ThresholdCategory_y = factor(as.character(DataToPlot_filtered[[paste0("ThreshCat_", 
                                                                                             SELECTED.MARKER[["y"]])]]), 
                                                    levels = sort(unique(DataToPlot_filtered[[paste0("ThreshCat_", 
                                                                                                     SELECTED.MARKER[["y"]])]]))),
                       ThresholdCategory_xy = factor(as.character(DataToPlot_filtered[["ThreshCat_xy"]]), 
                                                     levels = sort(unique(DataToPlot_filtered[["ThreshCat_xy"]]))),
                       Treatment = factor(as.character(DataToPlot_filtered$Treatment), 
                                          levels = Conditions[Conditions %in% unique(DataToPlot_filtered$Treatment)]),
                       TimePoint = factor(as.character(DataToPlot_filtered$TimePoint), 
                                          levels = sort(unique(DataToPlot_filtered$TimePoint))),
                       Cluster = factor(as.character(DataToPlot_filtered$Cluster), 
                                        levels = sort(unique(DataToPlot_filtered$Cluster))),
                       Cluster2 = if (CLUSTER.ROUND > 1) {
                         factor(as.character(DataToPlot_filtered$Cluster2), 
                                levels = sort(unique(DataToPlot_filtered$Cluster2)))
                       } else {
                         NULL
                       }
)
Biaxial_COLOR_BY <- COLOR_BY_value[COLOR_BY_biaxial]
# Histogram_PLOT_BY <- COLOR_BY_value[PLOT_BY_histogram]

# Add color info to DataToPlot_filtered dataframe
DataToPlot <- DataToPlot_filtered %>% 
  mutate(BiaxialColor = Biaxial_COLOR_BY[[1]])

# Define min/max values for each axis
MaxValue_y <- ceiling(max(DataToPlot[[SELECTED.MARKER[["y"]]]], na.rm = TRUE))
MinValue_y <- min(DataToPlot[[SELECTED.MARKER["y"]]], na.rm = TRUE) + (MaxValue_y * -0.025)

MaxValue_x <- ceiling(max(DataToPlot[[SELECTED.MARKER[["x"]]]], na.rm = TRUE))
MinValue_x <- min(DataToPlot[[SELECTED.MARKER["x"]]], na.rm = TRUE) + (MaxValue_x * -0.025)


#### Define function for different plots to generate ###########################
## Function to plot biaxials ===================================================
Biaxial_Function <- function(CurrentCluster) {
  ## For testing
  # CurrentCluster <- 4 # "PseudoBulk" / 4
  
  ## Set up necessary dataframes for plotting
  DataToPlot_CurrentCluster <- DataToPlot %>% 
    dplyr::filter(Cluster == CurrentCluster)
  
  Percentages_Cluster <- DataToPlot_CurrentCluster %>%
    count(Cluster, ThreshCat_xy, name = "n") %>%
    complete(Cluster,ThreshCat_xy,
             fill = list(n = 0)
             ) %>%
    group_by(Cluster) %>%
    mutate(Percent = 100 * n / sum(n)) %>%
    ungroup()
  
  
  # Define colors for each cluster to keep them consistent
  if (ColorToUse == "CellID") { 
    color_mapping <- setNames(scales::hue_pal()(length(unique(thresholding.ct.df.category[["Cluster"]]))),
                              sort(unique(thresholding.ct.df.category[["Cluster"]])))
  } else if (ColorToUse == "Default") {
    color_mapping <- setNames(hcl(h = seq(15, 375, 
                                          length = length(unique(DataToPlot[["BiaxialColor"]])) + 1), 
                                  l = 65, 
                                  c = 100)[1:length(unique(DataToPlot[["BiaxialColor"]]))],
                              as.character(sort(unique(DataToPlot[["BiaxialColor"]]))))
  } else {
    color_mapping <- viridis(n = length(unique(DataToPlot[["BiaxialColor"]])),
                             option = ColorToUse)
  }
  
  ## Save plot into BIAXIAL
  BIAXIAL <- DataToPlot_CurrentCluster %>% 
    ggplot(aes(x = .data[[SELECTED.MARKER[["x"]]]],
               y = .data[[SELECTED.MARKER[["y"]]]])) +
    geom_point(aes(color = .data[["BiaxialColor"]]), 
               size = point_size, 
               alpha = 0.5) + 
    geom_density_2d(bins = 10,
                    color = "black",
                    linewidth = 0.25,
                    alpha = 0.75) +
    # Make plot adjustments
    scale_color_manual(values = color_mapping[sort(unique(dplyr::filter(DataToPlot_CurrentCluster, 
                                                                        DataToPlot_CurrentCluster$Cluster == CurrentCluster)[["BiaxialColor"]]))]) +
    scale_x_continuous(limits = c(MinValue_x, MaxValue_x), 
                       breaks = seq(0, MaxValue_x, by = 1)) +
    scale_y_continuous(limits = c(MinValue_y, MaxValue_y),
                       breaks = seq(0, MaxValue_y, by = 1)) +
    # Add lines to denote EXPRESSION.DELINEATION.VALUES
    geom_vline(xintercept = as.vector(EXPRESSION.DELINEATION.VALUES[["x"]]),
               linetype = 'longdash', color = "black") +
    geom_hline(yintercept = as.vector(EXPRESSION.DELINEATION.VALUES[["y"]]),
               linetype = 'longdash', color = "black") +
    # Adjust labels
    labs(color = if (names(Biaxial_COLOR_BY) == "SignalingMarker") {
      MarkerToPlot
    } else if (grepl("Threshold", names(Biaxial_COLOR_BY))) {
      paste0("Cluster ", CurrentCluster , "\nThreshold\nCategory")
    } else {
      paste0("Cluster ", CurrentCluster, "\n", names(Biaxial_COLOR_BY))
    }
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  ## Label as necessary
  if (PercentsOnPlots) {
    label_x_coordinates <- c("11" = EXPRESSION.DELINEATION.VALUES[["x"]]*0.9, "12" = EXPRESSION.DELINEATION.VALUES[["x"]]*0.9,
                             "13" = EXPRESSION.DELINEATION.VALUES[["x"]]*0.9,
                             "21" = MaxValue_x*0.5, "22" = MaxValue_x*0.5, "23" = MaxValue_x*0.5,
                             "31" = MaxValue_x*0.9, "32" = MaxValue_x*0.9, "33" = MaxValue_x*0.9
                             )
    
    label_y_coordinates <- c("11" = 0, "12" = MaxValue_y*0.5, "13" = MaxValue_y*0.95,
                             "21" = 0, "22" = MaxValue_y*0.5, "23" = MaxValue_y*0.95,
                             "31" = 0, "32" = MaxValue_y*0.5, "33" = MaxValue_y*0.95
                             )
    
    BIAXIAL <- BIAXIAL +
      annotate(geom = "text",
               x = unname(label_x_coordinates[which(names(label_x_coordinates) %in% Percentages_Cluster$ThreshCat_xy)]),
               y = unname(label_y_coordinates[which(names(label_y_coordinates) %in% Percentages_Cluster$ThreshCat_xy)]),
               label = paste0(round(Percentages_Cluster$Percent, digits = 1), "%")
               )
    }
  
  ## Return plot
  Output_BiaxialPlot <- ggdraw() +
    draw_plot(BIAXIAL)
  
  return(Output_BiaxialPlot)
  }

## Function to plot histograms =================================================
Histogram_Function <- function(CurrentCluster, VaribaleToPlot) {
  ## For testing
  # CurrentCluster <- 4 # "PseudoBulk" / 4
  # VaribaleToPlot <- "ThreshCat_TrkB_i_174Yb" # "ThreshCat_TrkB_i_174Yb" / "ThreshCat_p75NTR_161Dy" / "ThreshCat_TrkB_s_158Gd" (from PLOT_BY_histogram)
  
  ## Define varaibles
  i <- SELECTED.MARKER[grepl(sub("ThreshCat_", "", VaribaleToPlot), SELECTED.MARKER)]
  MAX.SELECTED.MARKER <- if (names(i) == "x") { MaxValue_x } else if (names(i) == "y") { MaxValue_y }
  
  ## Set up necessary dataframes for plotting
  DataToPlot_CurrentCluster <- DataToPlot %>% 
    mutate(across(all_of(starts_with("ThreshCat_")),
                         ~ factor(.x, levels = sort(unique(.x)))
                  )) %>% 
    dplyr::filter(Cluster %in% CurrentCluster)
  
  # Percentages_Cluster <- DataToPlot_CurrentCluster %>% 
  #   group_by(across(all_of(c("Cluster", VaribaleToPlot))), .drop = FALSE) %>% 
  #   summarize(Percent = (n() / nrow(DataToPlot_CurrentCluster) * 100),
  #             .groups = "drop")
  
  Percentages_Cluster <- DataToPlot_CurrentCluster %>%
    count(Cluster, !!rlang::sym(VaribaleToPlot), name = "n") %>%
    complete(
      Cluster,
      !!rlang::sym(VaribaleToPlot),
      fill = list(n = 0)
    ) %>%
    group_by(Cluster) %>%  # or just RidgePlot_By for global denominator
    mutate(Percent = 100 * n / sum(n)) %>%
    ungroup()
  
  
  ## Save plot into HISTOGRAM
  if (!is.null(RidgePlot_By)) {
    # Define color values
    if (grepl("ThrechCat", RidgePlot_By)) { 
      ColorNum <- length(unique(DataToPlot[[RidgePlot_By]])) # Define number of total values
      colors <- scales::hue_pal()(length(ColorNum)) # Generate color palette
      ColorValues <- rev(colors[CurrentCluster]) # Select the only colors from the plotted clusters
      } else if (RidgePlot_By == "Cluster") {
      ColorNum <- length(unique(thresholding.ct.df.category[[RidgePlot_By]])) # Define number of total values
      colors <- rainbow(ColorNum, alpha = 1) # Generate color palette
      ColorValues <- rev(colors[CurrentCluster]) # Select the only colors from the plotted clusters
    } else {
      ColorNum <- length(unique(DataToPlot_CurrentCluster[[RidgePlot_By]])) # Define number of total values
      ColorValues <- rev(viridis(n = ColorNum, option = "C"))
      }
    
    # Generate histrogram plot
    HISTOGRAM <- DataToPlot_CurrentCluster %>% 
      mutate(!!RidgePlot_By := factor(.data[[RidgePlot_By]], levels = rev(sort(unique(.data[[RidgePlot_By]]))))) %>% 
      ggplot(aes(x = .data[[i]], y = .data[[RidgePlot_By]])) +
      geom_density_ridges(aes(fill = .data[[RidgePlot_By]]),
                          # quantile_lines = TRUE,                 # vertical quantile lines
                          # quantiles = c(0.25, 0.5, 0.75),        # Q1, median, Q3
                          # vline_linetype = c("dashed","solid","dashed"),
                          # vline_color = "black",
                          # vline_width = 0.75,
                          alpha = 0.75
                          ) +
      
      # facet_wrap(facets = vars(ThreshCat_TrkB_i_174Yb)) + 
      
      scale_fill_manual(values = ColorValues) +
      geom_vline(xintercept = EXPRESSION.DELINEATION.VALUES[[names(i)]],
                 linetype = 'longdash') +
      scale_x_continuous(limits = c(0, MAX.SELECTED.MARKER), breaks = seq(0, MAX.SELECTED.MARKER, by = 1)) + 
      labs(title = paste(i, ": Cluster ", CurrentCluster))
    
    if (!is.null(FacetBy)) {
      HISTOGRAM <- HISTOGRAM + facet_wrap(facets = vars(!!sym(FacetBy)), 
                                          labeller = labeller(.default = function(x) paste0(FacetBy, ": ", x))
                                          )
      }
    
    } else if (is.null(RidgePlot_By)) {
      HISTOGRAM <- DataToPlot_CurrentCluster %>% 
        ggplot(aes(x = .data[[i]])) +
        geom_histogram(aes(y = after_stat(density)),
                       fill = "lightgray",
                       bins = 175) +
        geom_density(lwd = 1.2,
                     position = "stack",
                     linetype = 2,
                     colour = "magenta"
                     ) +
        geom_vline(xintercept = EXPRESSION.DELINEATION.VALUES[[names(i)]],
                   linetype = 'longdash') + 
        labs(title = paste(i, ": Cluster ", CurrentCluster)) +
        # annotate(geom = "text",
        #          x = c(0, EXPRESSION.DELINEATION.VALUES[[names(i)]]),
        #          y = c(0.05 * 1:length(Percentages_Cluster[[VaribaleToPlot]])),
        #          # quantile(density(as.numeric(DataToPlot_CurrentCluster[[VaribaleToPlot]]))[["y"]], probs = 0.75, na.rm = TRUE) -
        #          # ((quantile(density(as.numeric(DataToPlot_CurrentCluster[[VaribaleToPlot]]))[["y"]], probs = 0.75, na.rm = TRUE)*0.1) *
        #          # 0:(length(Percentages_Cluster[[VaribaleToPlot]])-1)),
        #          label = paste0(Percentages_Cluster[[VaribaleToPlot]], "=\n" , round(Percentages_Cluster$Percent, digits = 1), "%")
        # ) +
        scale_x_continuous(limits = c(0, MAX.SELECTED.MARKER), 
                           breaks = seq(0, MAX.SELECTED.MARKER, by = 1))
      
      if (!is.null(FacetBy)) {
        HISTOGRAM <- HISTOGRAM + facet_wrap(facets = vars(!!sym(FacetBy)), 
                                            labeller = labeller(.default = function(x) paste0(FacetBy, ": ", x))
                                            )
        }
      
      } else {
        stop("Check that RidgePlot_By is filled appropriately")
        }
  
  ## Return plot
    HistPlot_Output <- ggdraw() +
      draw_plot(HISTOGRAM)
    
    return(HistPlot_Output)
    }
  
## Function to plot percentages in lolipop plots ===============================
Percent_Function <- function(CurrentCluster, VaribaleToPlot) {
  ## For testing
  # CurrentCluster <- 4 "PseudoBulk" # "PseudoBulk" / 4
  # VaribaleToPlot <- "ThreshCat_TrkB_i_174Yb" # "ThreshCat_TrkB_i_174Yb" / "ThreshCat_p75NTR_161Dy"  / "ThreshCat_TrkB_s_158Gd" (from PLOT_BY_histogram)
  
  ## Set up necessary dataframes for plotting
  DataToPlot_CurrentCluster <- DataToPlot %>% 
    dplyr::filter(Cluster == CurrentCluster) %>% 
    select(all_of(c("TimePoint", "Treatment", "Replicate", "Cluster")),
           all_of(starts_with("ThreshCat_")))
  
  ## Save plot into PercPlot
  if (!is.null(RidgePlot_By)) {
    Percentages_Cluster <- DataToPlot_CurrentCluster %>%
      select(all_of(c("TimePoint", "Treatment", "Replicate", "Cluster")),
             all_of(starts_with("ThreshCat_"))) %>% 
      group_by(across(everything())) %>%
      count(name = "n") %>%
      # complete(Cluster, !!rlang::sym(RidgePlot_By), !!rlang::sym(VaribaleToPlot),
      #          fill = list(n = 0)
      #          ) %>%
      group_by(Cluster, !!rlang::sym(RidgePlot_By)) %>%  # or just RidgePlot_By for global denominator
      mutate(Percent = 100 * n / sum(n),
             Cluster = factor(Cluster,  levels = sort(unique(Cluster))),
             !!VaribaleToPlot := factor(.data[[VaribaleToPlot]],  levels = sort(unique(.data[[VaribaleToPlot]]))),
             !!RidgePlot_By := factor(.data[[RidgePlot_By]],  levels = sort(unique(.data[[RidgePlot_By]])))
             ) %>%
      ungroup()
    
    Percentages_Cluster_summary <- Percentages_Cluster %>% 
      group_by(across(all_of(c("Cluster", "TimePoint", VaribaleToPlot)))) %>% 
      summarize(SD_Perc = sd(Percent),
                Percent = mean(Percent),
                .groups = "drop")
    
    
    PercPlot <- ggplot(data = Percentages_Cluster,
                       aes(x = .data[[RidgePlot_By]], y = Percent,
                           fill = .data[[VaribaleToPlot]])) +
      geom_point(aes(color = .data[[VaribaleToPlot]]),
                 size = 1,
                 # position = position_dodge(width = 0.9),
                 alpha = 0.75
                 ) + 
      # # If you prefer error bars to ribbon for plotting SD:
      # geom_errorbar(data = Percentages_Cluster_summary,
      #               aes(y = .data[["Percent"]],
      #                   ymin = .data[["Percent"]] - .data[["SD_Perc"]],
      #                   ymax = .data[["Percent"]] + .data[["SD_Perc"]]),
      #               width = 0.3,
      #               linewidth = 0.2
      #               ) +
      
      # If you prefer ribbon to error bars for plotting SD:
      geom_ribbon(data = Percentages_Cluster_summary,
                  aes(x = .data[[RidgePlot_By]], y = Percent, 
                      ymin = Percent - SD_Perc, ymax = Percent + SD_Perc, # Ribbon represents SD
                      fill = .data[[VaribaleToPlot]], 
                      group = !!sym(VaribaleToPlot)
                      ),
                  alpha = 0.2,
                  color = NA,
                  inherit.aes = FALSE
                  ) +
      geom_point(data = Percentages_Cluster_summary,
                 aes(color = .data[[VaribaleToPlot]]),
                 size = 3,
                 # position = position_dodge(width = 0.9),
                 alpha = 0.75
                 ) + 
      geom_line(data = Percentages_Cluster_summary,
                aes(y = Percent, group = !!sym(VaribaleToPlot), color = .data[[VaribaleToPlot]]),
                linewidth = 0.5
                ) +
      labs(title = paste("Cluster ", CurrentCluster)) # + 
      # theme(legend.position = "none") 
    
    } else if (is.null(RidgePlot_By)) { 
      Percentages_Cluster <- DataToPlot_CurrentCluster %>%
        count(Cluster, .data[[VaribaleToPlot]], name = "n") %>%
        complete(
          Cluster,
          .data[[VaribaleToPlot]],
          fill = list(n = 0)
          ) %>%
        group_by(Cluster) %>%  # or just RidgePlot_By for global denominator
        mutate(Percent = 100 * n / sum(n)) %>%
        ungroup()
      
      PercPlot <- ggplot(data = Percentages_Cluster,
                         aes(x = .data[[VaribaleToPlot]], y = Percent)) +
        geom_segment(aes(x = .data[[VaribaleToPlot]], xend = .data[[VaribaleToPlot]],
                         y = 0, yend = Percent)) +
        geom_point(aes(color = .data[[VaribaleToPlot]] %>% as.factor()),
                   size = 5) +
        geom_text(aes(label = round(Percent, digits = 2)),
                  nudge_y = 0,
                  nudge_x = 0.325,
                  size = 4) +
        labs(title = paste("Cluster ", CurrentCluster),
             color = VaribaleToPlot)
      
      }
  
  ## Return plot
  PercPlot_Output <- ggdraw() +
    draw_plot(PercPlot)

  return(PercPlot_Output)
  }



## Function to plot box-and-whisker plots ======================================
BoxAndWhisher_Function <- function(CurrentCluster, VaribaleToPlot) {
  ## For testing
  # CurrentCluster <- 4 "PseudoBulk" # "PseudoBulk" / 4
  # VaribaleToPlot <- "ThreshCat_TrkB_s_158Gd" # "ThreshCat_TrkB_i_174Yb" / "ThreshCat_p75NTR_161Dy" (from PLOT_BY_histogram)
  
  ## Define varaibles
  i <- SELECTED.MARKER[grepl(sub("ThreshCat_", "", VaribaleToPlot), SELECTED.MARKER)]
  MAX.SELECTED.MARKER <- if (names(i) == "x") { MaxValue_x } else if (names(i) == "y") { MaxValue_y }
  
  ## Set up necessary dataframes for plotting
  DataToPlot_CurrentCluster <- DataToPlot %>% 
    mutate(!!VaribaleToPlot := factor(.data[[VaribaleToPlot]],
                                      levels = sort(unique(.data[[VaribaleToPlot]])))) %>% 
    dplyr::filter(Cluster == CurrentCluster) %>% 
    select(all_of(c("TimePoint", "Treatment", "Replicate", "Cluster", unname(i), VaribaleToPlot))) %>% 
    mutate(across(all_of(c("TimePoint", "Treatment", "Replicate", "Cluster")),
                  ~ factor(.x, levels = sort(unique(.x)))
                  ))
  
  
  BoxAndWhiskerPlot <- ggplot(data = DataToPlot_CurrentCluster,
         aes(x = !!sym(unname(i)),
             y = Replicate)) + 
    geom_boxplot() + 
    facet_wrap(facets = vars(TimePoint),
               nrow = 6, ncol = 1) +
    scale_x_continuous(limits = c(0, MAX.SELECTED.MARKER), breaks = seq(0, MAX.SELECTED.MARKER, by = 1)) + 
    labs(title = paste(i, ": Cluster ", CurrentCluster))
  
  
  ## Return plot
  BoxAndWhiskerPlot_Output <- ggdraw() +
    draw_plot(BoxAndWhiskerPlot)
  
  return(BoxAndWhiskerPlot_Output)
  }


RunStatistics <- TRUE   # TRUE/FALSE
List_p_values <- TRUE   # TRUE/FALSE

StatisticalComparison_Type  <- "TimePoint"   # "Condition", "TimePoint", "Cluster"
StatisticalComparison_Value <- 0             # e.g. 0 for 0hr
Statistical_Test <- "t.test"                 # "t.test" or "wilcox"
Alpha <- 0.05

### Stats Function ============================================================
MakeStats_ConditionByTime <- function(df, marker, ref_treatment = "NoTreat", ref_timepoint = 0, test = "wilcox") {
  
  # df <- DataToPlot_Working # DataToPlot
  # marker <- "TrkB_s_158Gd"
  # ref_treatment <- "NoTreat"
  # ref_timepoint <- 0
  # test <- "t.test" Statistical_Test # "t.test"
  
  
  
  ref_vals <- df %>%
    dplyr::filter(Treatment == ref_treatment, as.character(TimePoint) == as.character(ref_timepoint)) %>%
    pull(.data[[marker]])
  
  if (length(ref_vals) < 2) {
    warning("Reference group has fewer than 2 observations - cannot run statistics")
    return(NULL)
  }
  
  df %>%
    # dplyr::filter(as.character(TimePoint) != as.character(ref_timepoint)) %>%
    # mutate(TimePoint_chr = as.character(TimePoint)) %>%
    group_by(across(all_of(c("TimePoint", "Cluster")))) %>%
    summarise(
      p_value = tryCatch({
        sub_vals <- pick(everything()) %>% pull(.data[[marker]])
        if (length(sub_vals) < 2) return(NA_real_)
        if (test == "t.test") {
          t.test(sub_vals, ref_vals, paired = FALSE)$p.value
        } else {
          wilcox.test(sub_vals, ref_vals, paired = FALSE)$p.value
        }
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      Significance = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      ), 
      REF = rep(ref_timepoint, time = .data[[StatisticalComparison_Type]] %>% length()),
      TEST = .data[[StatisticalComparison_Type]]
    )
}


## Function to plot box-and-whisker plots ======================================
LineGraph_Function <- function(CurrentCluster, CurrentMarker) {
  # CurrentCluster <- ClusterToPlot_linegraph # "PseudoBulk" / 4 / # HowManyLoops = ClusterToPlot / ClusterToPlot_linegraph
  # CurrentMarker <- "TrkB_i_174Yb" # "pERK_171Yb" / "TrkB_s_158Gd" / "TrkB_i_174Yb" (from MarkerToPlot)
  
  LineGraph_Output <- list()
  
  ## Set up necessary dataframes for plotting
  # Mutate necessary columns
  #, select those that will be used, and filter based on input parameters
  DataToPlot_Initial <- DataToPlot %>% 
    mutate(!!PLOT_BY_linegraph := factor(as.character(.data[[PLOT_BY_linegraph]]), 
                                         levels = sort(unique(.data[[PLOT_BY_linegraph]]))),
           TimePoint = factor(TimePoint, 
                              levels = sort(unique(TimePoint)))
           )  %>% 
    select(all_of(c("File_Num", "TimePoint", "Treatment", "Replicate", "Cluster", FacetBy, PLOT_BY_linegraph, CurrentMarker)))
  
  # Change PLOT_BY_linegraph if it is the same as the CurrentMarker (makes downstream calculations more flexible)
  if (PLOT_BY_linegraph == CurrentMarker) {
    PLOT_BY_linegraph <- NULL
    }
  
  # Take the average of all events per group
  DataToPlot_Initial_AvgPerRep <- DataToPlot_Initial %>% 
    group_by(across(all_of(all_of(c("TimePoint", "Treatment", "Replicate", "Cluster", FacetBy, PLOT_BY_linegraph))))) %>% 
    summarize(!!CurrentMarker := mean(.data[[CurrentMarker]]),
              .groups = "drop")
  
  # Normalize based on Normalize_LineGraph
  if (isFALSE(Normalize_LineGraph)) {
    DataToPlot_Working <- DataToPlot_Initial_AvgPerRep
  } else {
    # Take the means of TimePoint = 0 to be used for Fold Change normalization
    Means_0hr <- DataToPlot_Initial_AvgPerRep[DataToPlot_Initial_AvgPerRep[["TimePoint"]] == 0, ] %>% 
      group_by(across(all_of(all_of(c("TimePoint", "Cluster", FacetBy, PLOT_BY_linegraph))))) %>% 
      summarize(!!CurrentMarker := mean(.data[[CurrentMarker]]),
                .groups = "drop") %>% 
      rename_with(~ paste0(.x, "_Norm"),
                  .cols = -all_of(c("Cluster", FacetBy, PLOT_BY_linegraph)))
    
    # Take the Fold Change relative to TimePoint = 0
    DataToPlot_Working <- DataToPlot_Initial_AvgPerRep %>%
      # Append a TimePoint_Norm column to DataToPlot_PercAboveThresh
      mutate(TimePoint_Norm = factor(0, 
                                     levels = sort(unique(TimePoint)))
      ) %>% 
      # Append Means_0hr to DataToPlot_FC through Cluster, ThreshCat (if present), and TimePoint_Norm
      # Then add "_norm" to end of each appended column
      left_join(Means_0hr, 
                # by = c(colnames(DataToPlot_Initial_AvgPerRep)[!colnames(DataToPlot_Initial_AvgPerRep) %in% 
                #                                                        c("Treatment", "TimePoint", "Replicate", PLOT_BY_linegraph, CurrentMarker)]#,
                #        # "NormTimepoint" = "NormTimepoint"
                # ), 
                # suffix = c("", "_norm")
                ) %>%
      # Avoid division by zero
      mutate(across(ends_with("_Norm"), 
                    ~ ifelse(. == 0, NA, .))) %>% 
      # Calculate the Fold Change by dividing the DataToPlot_FC columns by the "_norm" columns and add "FC_" to front of each appended column
      mutate(across(all_of(CurrentMarker),
                    ~ . / get(paste0(cur_column(), "_Norm")),
                    .names = "FC_{.col}")) %>%
      # Select only the columns that underwent Fold Change transformation
      select(colnames(DataToPlot_Initial_AvgPerRep)[!colnames(DataToPlot_Initial_AvgPerRep) %in% CurrentMarker],
             starts_with("FC_")) %>% 
      # Remove "FC_" from the column names
      rename_with(~ sub("^FC_", "", .), starts_with("FC_")) %>%
      # Change any "0 / Norm = 0" values to NA
      mutate(across(all_of(CurrentMarker), 
                    ~ ifelse(. == 0, NA, .)))
    
    if (Normalize_LineGraph == "log2FC") {
      DataToPlot_Working <- DataToPlot_Working %>% 
        mutate(!!CurrentMarker := log2(.data[[CurrentMarker]]))
      }
    } 
  
  # Take the average for each replicate
  DataToPlot_Working_Summary <- DataToPlot_Working %>% 
    group_by(across(all_of(all_of(c("TimePoint", "Treatment", "Cluster", FacetBy, PLOT_BY_linegraph))))) %>% 
    summarize(StdDev = sd(.data[[CurrentMarker]]),
              !!CurrentMarker := mean(.data[[CurrentMarker]]),
              .groups = "drop")
  
  # Set min/max yaxis values
  if (LineGraph_Yaxis_Normalization == "Full_Dataset") {
    Max_Yaxis_Value <- max(thresholding.ct.df.category[[CurrentMarker]], na.rm = TRUE)
    Min_Yaxis_Value <- min(thresholding.ct.df.category[[CurrentMarker]], na.rm = TRUE)
  } else if (LineGraph_Yaxis_Normalization == "Conditions_Dataset") {
    Max_Yaxis_Value <- max(DataToPlot[[CurrentMarker]], na.rm = TRUE)
    Min_Yaxis_Value <- min(DataToPlot[[CurrentMarker]], na.rm = TRUE)
  } else if (LineGraph_Yaxis_Normalization == "Working_Dataset") {
    Max_Yaxis_Value <- max(DataToPlot_Working[[CurrentMarker]], na.rm = TRUE)
    Min_Yaxis_Value <- min(DataToPlot_Working[[CurrentMarker]], na.rm = TRUE)
  } else if (is.numeric(LineGraph_Yaxis_Normalization)) {
    Max_Yaxis_Value <- LineGraph_Yaxis_Normalization
    Min_Yaxis_Value <- ifelse(test = Normalize_LineGraph == "log2FC",
                              yes = -1 * LineGraph_Yaxis_Normalization,
                              no = 0)
  } else {
    Max_Yaxis_Value <- max(DataToPlot_Working[[CurrentMarker]], na.rm = TRUE)
    Min_Yaxis_Value <- min(DataToPlot_Working[[CurrentMarker]], na.rm = TRUE)
  }
  
  # Subset plotting dataframes based on CurrentCluster
  DataToPlot_Working <- DataToPlot_Working %>% 
    dplyr::filter(Cluster %in% CurrentCluster)
  DataToPlot_Working_Summary <- DataToPlot_Working_Summary %>% 
    dplyr::filter(Cluster %in% CurrentCluster)
  
  
  # Define colors for each cluster to keep them consistent
  if (ColorToUse == "CellID") { 
    # numclusters <- length(unique(thresholding.ct.df.category$Cluster)) # Calculate number of total clusters
    # color_mapping <- rainbow(numclusters, alpha = 1) # Generate color palette
    
    color_mapping <- setNames(scales::hue_pal()(length(unique(thresholding.ct.df.category[["Cluster"]]))),
                              sort(unique(thresholding.ct.df.category[["Cluster"]])))
    
    colors <- color_mapping[as.numeric(as.character(unique(DataToPlot_Working[[PLOT_BY_linegraph]])))] # CurrentCluster
  } else if (ColorToUse == "Default") {
    color_mapping <- hcl(h = seq(15, 375, 
                                 length = length(unique(DataToPlot_Working[[PLOT_BY_linegraph]])) + 1), 
                         l = 65, 
                         c = 100)[1:length(unique(DataToPlot_Working[[PLOT_BY_linegraph]]))]
    colors <- color_mapping[unique(DataToPlot_Working[[PLOT_BY_linegraph]])]
  } else {
    color_mapping <- viridis(n = length(unique(DataToPlot_Working[[PLOT_BY_linegraph]])),
                             option = ColorToUse)
    colors <- color_mapping[unique(DataToPlot_Working[[PLOT_BY_linegraph]])] # CurrentCluster
  }
  
  
  ## Save plot into LINEGRAPH
  LINEGRAPH <- DataToPlot_Working %>% 
    ggplot(aes(x = TimePoint, y = .data[[CurrentMarker]], 
               color = !!sym(PLOT_BY_linegraph))
           ) +
    geom_point(size = 1, alpha = 0.5) +
    # # If you prefer error bars to ribbon for plotting StdDev:
    # geom_errorbar(data = DataToPlot_summary,
    #               aes(y = !!sym(CurrentMarker), ymin = !!sym(CurrentMarker) - StdDev, ymax = !!sym(CurrentMarker) - StdDev),
    #               linewidth = 0.2
    #               ) +
    # If you prefer ribbon to error bars for plotting StdDev:
    geom_ribbon(data = DataToPlot_Working_Summary,
                aes(x = TimePoint, y = .data[[CurrentMarker]], 
                    ymin = !!sym(CurrentMarker) - StdDev, ymax = !!sym(CurrentMarker) + StdDev, # Ribbon represents StdDev
                    fill = !!sym(PLOT_BY_linegraph), 
                    group = !!sym(PLOT_BY_linegraph)
                    ), alpha = 0.2, color = NA, inherit.aes = FALSE
                ) + 
    geom_point(data = DataToPlot_Working_Summary,
               aes(y = .data[[CurrentMarker]]),
               size = 2
               ) + 
    geom_line(data = DataToPlot_Working_Summary, 
              aes(y = .data[[CurrentMarker]], group = !!sym(PLOT_BY_linegraph)),
              linewidth = 0.5
              ) + 
    labs(color = PLOT_BY_linegraph,
         fill = PLOT_BY_linegraph
         ) +
    ylab(ifelse(test = isFALSE(Normalize_LineGraph),
                yes = CurrentMarker, 
                no = paste0(Normalize_LineGraph, "(", CurrentMarker, " / 0hr)"))
         ) + 
    
    
    # geom_hline(yintercept = if (Normalize_LineGraph == "FoldChange") {
    #   c(1/1.5, 1.5)
    #   } else if (Normalize_LineGraph == "log2FC") {
    #     log2(c(1/1.5, 1.5))
    #     },
    #            linetype = 'longdash') + 
    
    
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    ylim(Min_Yaxis_Value, Max_Yaxis_Value) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(FacetBy)) {
    LINEGRAPH <- LINEGRAPH + facet_wrap(facets = vars(!!sym(FacetBy)), 
                                        labeller = labeller(.default = function(x) paste0(FacetBy, ": ", x))
                                        )
  }

  #######
  
  stats_df <- NULL
  if (RunStatistics) {
    stats_df <- MakeStats_ConditionByTime(
      DataToPlot_Working,
      CurrentMarker,
      ref_treatment = "NoTreat",
      ref_timepoint = 0,
      test = Statistical_Test
      
    )
  }
  
  # FIX: Add stats annotations OUTSIDE the ggplot '+' chain
  if (!is.null(stats_df) && nrow(stats_df) > 0 && any(stats_df$Significance != "")) {
    stats_df$TimePoint <- factor(
      as.character(stats_df$TimePoint),
      levels = levels(DataToPlot_Working$TimePoint)
    )
    # stats_df$y <- Max_Yaxis_Value * 1.06
    for (i in 1:length(unique(stats_df[["TEST"]]))) {
      stats_df[which(stats_df[["Cluster"]] == unique(stats_df[["Cluster"]])[i]), "y"] <- (Max_Yaxis_Value * 0.9) - (0.15*(i-1)) 
      # (Max_Yaxis_Value * 1.06) - (i-1) / (Max_Yaxis_Value * 0.8) - (0.25*(i-1))  / (Max_Yaxis_Value * 1.06) - (2*(i-1)) / (Max_Yaxis_Value * 0.8) - (5*(i-1)) 
    }
    
    LINEGRAPH <- LINEGRAPH +
      geom_text(
        data = stats_df %>% dplyr::filter(Significance != ""),
        aes(x = TimePoint, y = y, label = Significance, color = .data[["Cluster"]]),
        inherit.aes = FALSE,
        size = 10
      )
  }
  
  if (List_p_values) {
    
    TablePlot <- grid.arrange(top = paste0(CurrentMarker, ": ",
                                           ifelse(test = length(CurrentCluster) > 0,
                                                  yes = paste0("Cluster ", paste0(CurrentCluster, collapse = ",")), 
                                                  no = "Pseudobulk")
                                           ), 
                              tableGrob(d = stats_df[, c("TimePoint", "Cluster", "p_value")],
                                        theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.1, fontsize = 10)),
                                                               rowhead = list(fg_params = list(fontface = "bold", fontsize = 10)))))
    
    Output_Table <- ggdraw() + draw_plot(TablePlot)
    
    LineGraph_Output[[paste0("Table: ", CurrentMarker)]] <- Output_Table
  }
  
  ########
  
  ## Return plot
  LineGraph_Output[[CurrentMarker]] <- ggdraw() + draw_plot(LINEGRAPH)
  
  return(LineGraph_Output)
}


#### Define variables for running functions ####################################
## Define how many times to apply the Plotting_Function (based on if data is PseudoBulk or clustered)
if (ClusterSelection == "PseudoBulk") {
  HowManyLoops <- "PseudoBulk"
} else if (ClusterSelection != "PseudoBulk") {
  HowManyLoops <- ClusterToPlot
}

## Define names for plots to generate
SelectedPlotsToGenerate <- c(PlotsToGenerate[which(!PlotsToGenerate == "Histogram" & !PlotsToGenerate == "BoxAndWhisker")], 
                             rep(unname(PLOT_BY_histogram), times = length(which(PlotsToGenerate == "Histogram" | PlotsToGenerate == "BoxAndWhisker")))
                             )
names(SelectedPlotsToGenerate) <- c(PlotsToGenerate[which(!PlotsToGenerate == "Histogram" & !PlotsToGenerate == "BoxAndWhisker")],
                                    sapply(PlotsToGenerate[which(PlotsToGenerate == "Histogram" | PlotsToGenerate == "BoxAndWhisker")], function (n) {
                                      paste0(n, "_", 1:length(PLOT_BY_histogram))
                                      }) %>% as.vector()
                                    ) %>% unlist()


#### Run functions and save plots ##############################################
sapply(names(SelectedPlotsToGenerate), function(PLOT_TYPE_name) {
  ## For testing
  # PLOT_TYPE_name <- "BoxAndWhisker" # "Biaxial" / "Percentage" / "BoxAndWhisker_1" / "BoxAndWhisker_2" / "Histogram_1" / "Histogram_2" / "LineGraph"
  
  ## Define which type of plot to set up
  PLOT_TYPE <- SelectedPlotsToGenerate[PLOT_TYPE_name] 
      # from names(SelectedPlotsToGenerate) / ("Biaxial", "Percentage", "BoxAndWhisker_1", "BoxAndWhisker_2", "Histogram_1", "Histogram_2", "LineGraph")
  
  ## Run plotting function for given PLOT_TYPE
  Plot_List <- list()
  if (names(PLOT_TYPE) == "Biaxial") {
    OUTPUT_FORMAT <- BIAXIAL_OUTPUT_FORMAT
    for (j in HowManyLoops) {
      Plot_List[[paste0("Cluster", j)]] <- Biaxial_Function(CurrentCluster = j)
      }
    } else if (grepl("Histogram", names(PLOT_TYPE))) {
      OUTPUT_FORMAT <- HISTOGRM_OUTPUT_FORMAT
      if (is.null(RidgePlot_By)) {
        for (j in HowManyLoops) {
          Plot_List[[paste0("Cluster", j)]] <- Histogram_Function(CurrentCluster = j,
                                                                  VaribaleToPlot = unname(PLOT_TYPE))
        }
      } else if (RidgePlot_By != "Cluster") {
        for (j in HowManyLoops) {
          Plot_List[[paste0("Cluster", j)]] <- Histogram_Function(CurrentCluster = j,
                                                                  VaribaleToPlot = unname(PLOT_TYPE))
        }
      } else {
        Plot_List[[paste0("Cluster", paste0(HowManyLoops, collapse = ","))]] <- Histogram_Function(CurrentCluster = HowManyLoops,
                                                                                                   VaribaleToPlot = unname(PLOT_TYPE))
      }
      } else if (names(PLOT_TYPE) == "Percentage") {
        OUTPUT_FORMAT <- PERCENT_OUTPUT_FORMAT
        for (j in HowManyLoops) {
          Plot_List[[paste0("Cluster", j)]] <- Percent_Function(CurrentCluster = j,
                                                                VaribaleToPlot = "ThreshCat_xy")
          }
        } else if (grepl("BoxAndWhisker", names(PLOT_TYPE))) {
          OUTPUT_FORMAT <- BOXandWHISKER_OUTPUT_FORMAT
          for (j in HowManyLoops) {
            Plot_List[[paste0("Cluster", j)]] <- BoxAndWhisher_Function(CurrentCluster = j,
                                                                        VaribaleToPlot = unname(PLOT_TYPE))
            }
          } else if (names(PLOT_TYPE) == "LineGraph") {
            OUTPUT_FORMAT <- LineGraph_OUTPUT_FORMAT
            if (length(MarkerToPlot) > 1 & length(ClusterToPlot_linegraph) > 1) {
              print("Check that either MarkerToPlot or ClusterToPlot has only one value")
              setwd(INPUT.FOLDER)
              } else if (length(MarkerToPlot) > 1 | PLOT_BY_linegraph == "Cluster") {
                for (j in MarkerToPlot) {
                  Plot_List[[j]] <- LineGraph_Function(CurrentCluster = ClusterToPlot_linegraph,
                                                       CurrentMarker = j)
                  }
                } else {
                  for (j in HowManyLoops) {
                    Plot_List[[paste0("Cluster", j)]] <- LineGraph_Function(CurrentCluster = j,
                                                                            CurrentMarker = MarkerToPlot)
                    }
                  }
            }
  
  ## Remove any blank plots
  Plot_List_FINAL <- Plot_List[nchar(names(Plot_List)) > 0]
  
  if (length(Plot_List_FINAL) == 1) { 
    PlotAsGrid <- FALSE
    }
  
  ## Output plots based on user-defined parameters
  if(PlotAsGrid) { # Generate grid, if selected, then save
    ## Arrange plots
    # Define plot dimensions for the outputted plot
    plot_height <- size_per_plot * ceiling(length(Plot_List_FINAL) / plot_column_length)
    plot_width <- size_per_plot * plot_column_length
    
    # Pad the plot list if needed (ex. if number of plots is not divisible by plot_column_length)
    remainder <- length(Plot_List_FINAL) %% plot_column_length
    if (remainder != 0) {
      NumberOfGridSlots <- ceiling(length(Plot_List_FINAL) / plot_column_length) * plot_column_length
      for (k in (length(Plot_List_FINAL) + 1):NumberOfGridSlots) {
        Plot_List_FINAL[[k]] <- ggplot() + theme_void()
      }
    }
    
    # Arrange and add title
    PlotsInGrid <- grid.arrange(grobs = Plot_List_FINAL,
                                nrow = ceiling(length(Plot_List_FINAL) / plot_column_length),
                                ncol = plot_column_length)
    
    grid_plot <- grid.arrange(textGrob(paste0(sub("(*_*)([^_]*)$", "", SELECTED.MARKER[["x"]]), 
                                              " -vs- ", sub("(*_*)([^_]*)$", "", SELECTED.MARKER[["y"]])), 
                                       gp = gpar(fontsize = 18, fontface = "bold")), 
                              PlotsInGrid, 
                              ncol = 1, heights = c(0.05, 0.9))
    
    ## Save arranged plots
    ggsave(filename = file.path(paste0(INPUT.FOLDER, "/", output_dir),
                                paste0(sub("_.", "", names(PLOT_TYPE)), "_",  # plot info
                                       paste0(sub("(*_*)([^_]*)$", "", SELECTED.MARKER[["x"]]), 
                                              "-vs-", 
                                              sub("(*_*)([^_]*)$", "", SELECTED.MARKER[["y"]])), # marker info
                                       "_Conditions-", paste0(ConditionsToPlot, collapse = ","), # Conditions info
                                       ifelse(test = ClusterSelection == "Individual",
                                              yes = paste0("_Cluster", paste0(ClusterToPlot, collapse = ",")),
                                              no = paste0("_", ClusterSelection)), # Cluster info
                                       "_ColorBy", SelectedPlotsToGenerate[[PLOT_TYPE_name]],  # Color info
                                       "-Grid", # Indicate grid
                                       ifelse(test = PercentsOnPlots, 
                                              yes = "-PlotPerc",
                                              no = ""),
                                       if (!is.null(RidgePlot_By) & grepl("Histogram", names(PLOT_TYPE))) {
                                         paste0("RigdePlot-", RidgePlot_By)
                                         }, # Info if plotting Ridge Plot
                                       "_", Statistical_Test,
                                       ".", LineGraph_OUTPUT_FORMAT) # file type info
                                ),
    plot = grid_plot,
    width = plot_width, height = plot_height,
    units = "in", dpi = 300
    )
    
    } else { # Or save each file
      ## Loop through all plots in Plot_List_FINAL and save individually
      
      if (length(unlist(Plot_List_FINAL)) > length(Plot_List_FINAL)) {
        ExportNames <- sub(paste0(names(Plot_List_FINAL), "."),  "", names(unlist(Plot_List_FINAL)))
        } else {
          ExportNames <- names(Plot_List_FINAL)
          }
      
      
      for (k in ExportNames) { # seq(length(Plot_List_FINAL)
        
        ggsave(filename = file.path(paste0(INPUT.FOLDER, "/", output_dir),
                                    paste0(sub("_.", "", names(PLOT_TYPE)), "_",  # plot info
                                           paste0(sub("(*_*)([^_]*)$", "", SELECTED.MARKER[["x"]]), 
                                                  "-vs-", 
                                                  sub("(*_*)([^_]*)$", "", SELECTED.MARKER[["y"]])), # marker info
                                           "_Conditions-", paste0(ConditionsToPlot, collapse = ","), # Conditions info
                                           ifelse(test = ClusterSelection == "Individual",
                                                  yes = paste0("_Cluster", paste0(ClusterToPlot, collapse = ",")),
                                                  no = paste0("_", ClusterSelection)), # Cluster info
                                           "_ColorBy", names(SelectedPlotsToGenerate[PLOT_TYPE_name]),  # Color info
                                           ifelse(test = PercentsOnPlots, 
                                                  yes = "-PlotPerc",
                                                  no = ""),
                                           if (!is.null(RidgePlot_By) & grepl("Histogram", names(PLOT_TYPE))) {
                                             paste0("RigdePlot-", RidgePlot_By)
                                             }, # Info if plotting Ridge Plot
                                           ifelse(test = !is.null(length(ClusterToPlot_linegraph)),
                                                  yes = paste0("Cl", paste0(ClusterToPlot_linegraph, collapse = ",")),
                                                  no = ""
                                                  ),
                                           "_", Statistical_Test,
                                           "_", k, # 
                                           ".", LineGraph_OUTPUT_FORMAT) # file type info
                                    ),
        plot = ifelse(test = length(unlist(Plot_List_FINAL)) > length(Plot_List_FINAL),
                      yes = Plot_List_FINAL[[1]][k],
                      no = Plot_List_FINAL[[k]]),
        width = OUTPUT_WIDTH, height = OUTPUT_HEIGHT,
        units = "in", dpi = 300
        )
        
        
        
        
      }
    }
  })


#### Reset working directory ###################################################
setwd(INPUT.FOLDER)

print("Complete! Have fun staring at your new plot!")





# ### Other ways to plot signaling ##############################################
#
# # Density Plot colored by a signaling marker
# DENSITY.PLOT <- DataToPlot %>%
#   group_by(TimePoint) %>%
#   ggplot(aes(x = .data[[SELECTED.MARKER[["x"]]]],
#              y = .data[[SELECTED.MARKER[["y"]]]]#,
#              # group = TimePoint,
#              # color = factor(TimePoint,
#              #                levels = sort(unique(TimePoint))),
#              # fill = factor(TimePoint,
#              #               levels = sort(unique(TimePoint))))
#              )
#          ) +
#   
#   geom_bin2d(aes(color = .data[[MarkerToPlot]]))
#   
#   
#   
#   geom_density_2d(aes(color = .data[[MarkerToPlot]]),
#                   alpha = 0.05) +
#   labs(color = "Legend",
#        fill = "Legend")
# 
# 
# output_density_filename <- paste0("DensityPlot_Conditions-", paste0(ConditionsToPlot, collapse = ","), "_Clusters-", paste0(ClusterToPlot, collapse = ","))
# 
# ggsave(filename = output_density_filename,
#        plot = DENSITY.PLOT,
#        device = OUTPUT_FORMAT,
#        height = OUTPUT_HEIGHT, width = OUTPUT_WIDTH)
# 
# setwd(INPUT.FOLDER)
# 


