# Updated script:
# Jonathon Sewell
# Version: 4
# 2026-02-23

# Source script: Figure_7A_View_3D_UMAP.R
# Amy Van Deusen, University of Virginia
# 2022
# Visualize 3D UMAP layout with plotly
# See https://plotly.com/r/3d-scatter-plots/#basic-3d-scatter-plot 
# This will allow for .html and other outputs to be interactive for sample metadata

# Updates to Source script:
# (1) Adjusted naming: "plot.df.final" to "plot.df.final.indexed" for subsetting (subsetted dataframe reverts back to "plot.df.final")
# (2) With adjusted naming, incorporated ability to subset dataframe if it is too large for plotly to plot
# (3) Puts values that were removed into a new dataframe that can be outputted, if necessary
# (4) Adjust output html name to indicate subsetting
#     ------------------------------------------------------------------------------------------------
# (5) Enables ability to "gray out" certain clusters while keeping the same colors for select clusters
# (6) Allows coloring by relevant metadata (made Cluster a factor and added 'layout' to plotly to adjust legend size)
# (7) Allows choice to make select clusters "red" to help indicate which events are being analyzed
# (8) Adjust output html name to indicate coloring
#     ------------------------------------------------------------------------------------------------
# (9) Enabled ability to save as a static image
# (10) Adjusted outputs when PLOT.COLOR.FACTOR is not NULL
#     ------------------------------------------------------------------------------------------------
# (11) Enables option for merged clusters
# (12) Includes ability to define a specific color palette

# NOTE: 
# • Incorporated subsetting because it seems that plotly doesn't like data frames that are too large (plots grid but no data points)
# • If there is an easy fix that I am just unaware of, the subsetting may be unnecessary
# • Updates were assisted by Microsoft Copilot AI to improve updated sections for readability and efficiency
      # Prompt: "Please improve the following code for readability and efficiency while maintaining the goal of the code"


print("Start View_3D_UMAP.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

## Load libraries ==============================================================
library(data.table)
# install.packages("plotly")
library(plotly)
library(reticulate)
library(umap)
library(dplyr)
library(RColorBrewer) # Optional: depends on how color plot
library(viridis) # Optional: depends on how color plot
library(htmlwidgets) # Optional: necessary if output HTML file
library(circlize)

### Input parameters ===========================================================
## Identify whether to subset dataframe and by what percentage
  DownSample <- FALSE # TRUE / FALSE
  DownSample_Percent <- 90

## Define which identified clusters to plot and which layout to use
  CLUSTER_ROUND <- 1
  CLUSTERS_MERGED <- FALSE # TRUE / FALSE
  LAYOUT_ROUND <- 1
  
## Identify clusters to color while all others are "grayed out"
  CLUSTERS.TO.COLOR <- sort(c()) # if coloring all clusters, use c() - unless coloring all red, use c(1:22)
  
# Do you want to use a single color to denote the selected cluster(s)
  CLUSTER.SINGLE.COLOR <- FALSE # TRUE / FALSE
  SINGLE.COLOR <- "red" # "lightblue4" "#440154FF"  "#21908CFF" / "red"
  
  # Define the order of clusters (useful if grouping by putative cell types)
  # ClusterOrder <- c(2, 4, 11, 1, 14, 15, 5, # Gradient A
  #                   12, 10, 6, # Gradient B
  #                   19, 21, # Astrocytic
  #                   8, 18, # Microglial
  #                   22, # Oligodendrocytic
  #                   16, 20, # Other Progenitor-like
  #                   17, # Misc. Neuron
  #                   7, 13, # Misc. Glia
  #                   3, 9) # Unknown
  # See "2025-08-12_cluster_R1_groups_JMS.csv" and "2025-08-12 Microsoft Copilot Discussion_Cluster ID.docx"
  # Directory: '/Users/jmsewell/Desktop/2025-05-13 BDNF_Rescue_FullSampleSet/2.2 Analysis_V3 (157Gd Gated_V2_Un-compensated)/
  # 4.1.0 Gated_157GdRemoved_BatchCorrection_2025-05-21/4.1 AllFiles_2025-05-21/quantile/2025-07-29 Separate Timepoints for clustering/
  # N. 0-16hr_SpecificitySamples'
  # c() if ClusterOrder is to go un-defined

## Metadata choices
# Choose which, if any, metadata variables by which to color 
  PLOT.COLOR.FACTOR <- NULL # NULL defaults to Cluster / Other metadata values: "Treatment" / "TimePoint"
# List of sample metadata to include in plot dataframe (added for each cell)
  METADATA.TO.INCLUDE <- c("Expt", "Treatment", "TimePoint") 
# Define the color scale to be used for selected metadata
  ColorPaletteToUse <- "Rainbow" # viridis / RColorBrewer / Rainbow / UserDefined
  UserDefinedPalette <- c("#F7FBFF", #  1  pale blue          (Blue gradient 1/6)
                          "#C6DBEF", #  2  mid-light blue     (Blue gradient 3/6)
                          "#89D0C3", #  3  light blue-green   (Light green -> light blue-green, step 2/2)  [CB-safe teal tint]
                          "#084594", #  4  dark blue          (Blue gradient 6/6)                           [CB-safe deep blue]
                          "#D9D9D9", #  5  gray (light)
                          "#00735A", #  6  blue-green (darker)    (Blue-green gradient 2/3)
                          "#DEEBF7", #  7  very light blue     (Blue gradient 2/6)
                          "#004B3C", #  8  blue-green (darkest)   (Blue-green gradient 3/3)
                          "#D55E00", #  9  red (vermillion)       [Okabe–Ito red; CB-safe alternative to pure red]
                          "#B3B3B3", # 10  gray (mid)
                          "#9ECAE1", # 11  medium blue         (Blue gradient 4/6)
                          "#009E73", # 12  blue-green (medium)    (Blue-green gradient 1/3)                  [Okabe–Ito bluish green]
                          "#6BAED6", # 13  deep blue           (Blue gradient 5/6)
                          "#7F7F7F", # 14  gray (dark)
                          "#B9E769", # 15  light green         (Light green -> light blue-green, step 1/2)  [yellow‑green for CB safety]
                          "#4D4D4D", # 16  gray (darker)
                          "#CC79A7", # 17  purple               [Okabe–Ito reddish purple]
                          "#E69F00", # 18  orange               [Okabe–Ito orange]
                          "#F0E442", # 19  yellow               [Okabe–Ito yellow; very distinct hue]
                          "#8C564B"  # 20  brown                [earth tone; distinct from others]
                          )
  
  COLOR_SCALE_ToUse <- c(12, "Paired")
  # With viridis, A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo
  # With RColorBrewer, use display.brewer.all(colorblindFriendly = TRUE) to identify palette codes
  # use c(# , name), where '#' is the number of colors from 'name' palette to be used for colorRampPalette interpolation
  # Ex. c(8, "Set2") / c(12, "Paired") / c(8, "Dark2") / c(11, "RdYlBu") / c(11, "PiYG") / c(11, "BrBG")
  
# Given the PLOT.COLOR.FACTOR, specify an order by which values will be reported
  if (is.null(PLOT.COLOR.FACTOR)) {
    METADATA_ORDER <- NULL
    PLOT.COLOR.FACTOR <- "Cluster"
    } else if (PLOT.COLOR.FACTOR == "Treatment") {
      METADATA_ORDER <- c("NoTreat", 
                          "ContDepr", "CompMed", "BDNF", 
                          "BDNFTest", "BDNF+K252aTest"
                          )
      } else if (PLOT.COLOR.FACTOR == "TimePoint") {
        METADATA_ORDER <- c("0", "0.083","0.25","1","4","16"
                            #,"72"
        )
        }
  # ^^ This needs to match with the values in metadata.csv in the PLOT.COLOR.FACTOR column
  # Otherwise, set METADATA_ORDER to NULL (below) and metadata will be sorted alphabetically
  # METADATA_ORDER <- NULL

## Make choices on how to plot data
# Do you want to plot the resulting graph in the RStudio viewer (good to make adjustments to orientation, etc.)?
  plot_in_viewer <- FALSE # TRUE / FALSE

# Has the inputted dataframe already been subsampled prior to running this script?
  SUBSAMPLE.IDS.FILENAME <- NULL # NULL if no subsample IDs ("Subsample_ID_by_Cluster.csv")
# Is there already a file with color information?
  CLUSTER.COLORS <- NULL # .csv file ("colors_R1.csv") or NULL will use default colors

## Input file information
  INPUT.FOLDER <- getwd()
  OUTOUT.FOLDER <- INPUT.FOLDER
  UMAP.3D.LAYOUT <- sub("X", LAYOUT_ROUND, "umap_xyz_3D_RX.csv")
  EXPRESSION.FILENAME <- "expression_matrix_analysis.csv" # "Concat_Transformed_by_Cluster.csv"
  CLUSTERS.FILENAME <- sub("X", CLUSTER_ROUND, "cluster_RX_assigns.csv") # "Clusters_Subsampled.csv"
  if (CLUSTERS_MERGED) {
    CLUSTERS.FILENAME <- sub(".csv", "_merged.csv", CLUSTERS.FILENAME)
  }
  METADATA.FILENAME <- "metadata.csv"
  FILENUMS.FILENAME <- "filenums.csv" # "filenums_Subsampled.csv"

## Output file information
  UMAP.FILENAME <- if (DownSample) { # Naming for downsampling
    paste0("UMAP_3D_", DownSample_Percent, "percSubset", 
           "_L", LAYOUT_ROUND, "_C", CLUSTER_ROUND)
  } else {
    paste0("UMAP_3D_All", 
           "_L", LAYOUT_ROUND, "_C", CLUSTER_ROUND)
  }
  OUTPUT.TYPE <- "pdf" # Choices are "svg" (preferred), "pdf", "png", etc.
  OUTPUT.Height <- 1200
  OUTPUT.Width <- 1600
  
  SAVE.HTML <- FALSE # TRUE = output .html file, FALSE = do not output file
  HTML.FILENAME <- if (DownSample) { # Naming for downsampling
    paste0("UMAP_3D_", DownSample_Percent, "percSubset", "_L", LAYOUT_ROUND, "_C", CLUSTER_ROUND, ".html")
    } else {
      paste0("UMAP_3D_All", "_L", LAYOUT_ROUND, "_C", CLUSTER_ROUND, ".html")
      }
  HTML.FILENAME <- if (length(CLUSTERS.TO.COLOR) > 0) { # Naming for coloring specific clusters
    
    if (CLUSTER.SINGLE.COLOR) {
      paste0(gsub(".html", "", HTML.FILENAME), "_ColorCluster", paste0(CLUSTERS.TO.COLOR, collapse = ","), "-", SINGLE.COLOR, ".html")
    } else {
      paste0(gsub(".html", "", HTML.FILENAME), "_ColorCluster", paste0(CLUSTERS.TO.COLOR, collapse = ","), ".html")
    }
  } else {
    HTML.FILENAME
  }
  
  FOLDER_BASENAME <- "Marker_Expression_VAR_PERpercSubset_XYZ"
  if (LAYOUT_ROUND > 1) {
    output_dir <- paste0("05.1-", # Script ID
                         format(Sys.time(), "%Y-%m-%d"), "_", # Date info
                         sub("VAR_PER", paste0("L", LAYOUT_ROUND,
                                               "_G", LAYOUT_GROUP,
                                               "_", DownSample_Percent),
                             FOLDER_BASENAME))
  } else {
    output_dir <- paste0("05.1-", # Script ID
                         format(Sys.time(), "%Y-%m-%d"), "_", # Date info
                         sub("VAR_PER", paste0("L", LAYOUT_ROUND,
                                               "_", DownSample_Percent),
                             FOLDER_BASENAME))
  }

## Read necessary files ========================================================
expression.filename <- paste0(INPUT.FOLDER, "/", EXPRESSION.FILENAME)
exprs.in <- fread(expression.filename, header=TRUE, sep=",", check.names=FALSE)
layout.filename <- paste0(INPUT.FOLDER, "/", UMAP.3D.LAYOUT)
layout.in <- read.table(layout.filename, header=TRUE, sep=",", check.names=FALSE)
clusters.filename <- paste0(INPUT.FOLDER, "/", CLUSTERS.FILENAME)
clusters.in <- read.table(clusters.filename, header=TRUE, sep=",", check.names=FALSE)
metadata.filename <- paste0(INPUT.FOLDER, "/", METADATA.FILENAME)
metadata.in <- read.table(metadata.filename, header=TRUE, sep=",", check.names=FALSE)
filenums.filename <- paste0(INPUT.FOLDER, "/", FILENUMS.FILENAME)
filenums.in <- read.table(filenums.filename, header=TRUE, sep=",", check.names=FALSE)
if (!is.null(SUBSAMPLE.IDS.FILENAME)) {
  subsamples.filename <- paste0(INPUT.FOLDER, "/", SUBSAMPLE.IDS.FILENAME)
  subsamples.in <- read.table(subsamples.filename, header=FALSE, sep=",", check.names=FALSE)
}

if (is.null(CLUSTERS.TO.COLOR)) {
  CLUSTERS.TO.COLOR <- sort(unique(clusters.in[["Cluster"]]))
  }

## Prepare plotting dataframe ==================================================
# Combine layout, clusters, and filenums into single dataframe
colnames(filenums.in) <- "File"
plot.df.init <- cbind(filenums.in, layout.in) # Combine file numbers and UMAP layout
colnames(clusters.in) <- "Cluster"
plot.df.init2 <- cbind(plot.df.init, clusters.in) # Add clusters to file numbers and UMAP layout

# Get cell IDs or subsample IDs and add to plotting dataframe
if (!is.null(SUBSAMPLE.IDS.FILENAME)) {
  cell.ids <- subsamples.in 
} else {
  cell.ids <- as.data.frame(1:length(clusters.in$Cluster))
}
colnames(cell.ids) <- "cell.id"
plot.df <- cbind(plot.df.init2, cell.ids)

# Prepare sample metadata (specified in METADATA.TYPE) and add to plotting dataframe
metadata.out <- vector(mode = "list") # Initialize empty list to collect metadata
filenums <- as.vector(filenums.in$File) # Convert file numbers to vector for subsetting below

for (i in METADATA.TO.INCLUDE) {
  metadata_key <- as.character(metadata.in[,i]) # Generate metadata for selected metadata
  names(metadata_key) <- 1:length(metadata.in[,i]) # Generate new index for selected metadata
  metadata.out[[i]] <- metadata_key[filenums] # Add selected metadata to compiled list
}
metadata.final <- as.data.frame(metadata.out) # Convert selected metadata into dataframe
plot.df.final.indexed <- cbind(plot.df, metadata.final) # Add metadata to plotting dataframe (quick)
plot.df.exprs <- cbind(plot.df.final.indexed, exprs.in) # Generate larger plotting dataframe with expression data



# Subset (as necessary) data to help plot datasets that are too large for plotly
set.seed(123) # Make downsampling reproducible
if (DownSample) {
  # Downsample
  plot.df.final <- plot.df.final.indexed %>%
    group_by(File) %>%
    slice_sample(prop = (DownSample_Percent/100)) %>%
    ungroup()
  
  # Create the removed rows dataframe
  plot.df.final.downsampled.removed <- anti_join(plot.df.final.indexed, plot.df.final, by = names(plot.df.final.indexed))
  
} else {
  plot.df.final <- plot.df.final.indexed
}


## Generate plotly plot ========================================================
dir.create(output_dir)
setwd(output_dir)

## Load necessary information to use save_image() in the loop
# Find path to python in Terminal or command prompt:
# which python
# Or in Windows:
# where python
# Or in R:
# py_config()

# Specify which python installation to use
reticulate::use_python("/Users/jmsewell/opt/anaconda3/bin/python")

# Ensure that the kaleido package in python is installed
# py_install("kaleido")
# Or from terminal:
# pip install kaleido==0.2.1 - I found this to work (py_intall couldn't find package)
# Note: kaleido verion 0.2.1 appears to be better suited for save_image()
# Verify installation
py_module_available("kaleido")


# Colored by clusters (Default)

if (PLOT.COLOR.FACTOR == "Cluster") {
  if (ColorPaletteToUse == "UserDefined") {
    numclusters <- length(unique(plot.df.final$Cluster)) # Calculate number of clusters to color
    colors.final <- UserDefinedPalette # Generate color palette: any RColorBrewer or viridis palette
    colors.final[-CLUSTERS.TO.COLOR] <- "gray"
    } else if (!is.null(CLUSTER.COLORS)) {
      colors.filename <- paste0(INPUT.FOLDER, "/", CLUSTER.COLORS) # Find colors.csv
      colors.in <- read.csv(colors.filename) # Read colors.csv
      colors.final <- as.vector(colors.in$color) # Convert colors to vector
      names(colors.final) <- unique(clusters.in)
      } else if (length(CLUSTERS.TO.COLOR) > 0) {
        if (CLUSTER.SINGLE.COLOR) {
          numclusters <- length(unique(plot.df.final$Cluster)) # Calculate number of clusters to color
          colors.final <- rep(SINGLE.COLOR, times = numclusters) # Generate a color palette with SINGLE.COLOR for each Cluster
          colors.final[-CLUSTERS.TO.COLOR] <- "gray" # Gray out clusters not defined in CLUSTERS.TO.COLOR
          } else {
            numclusters <- length(unique(plot.df.final$Cluster)) # Calculate number of clusters to color
            colors.final <- rainbow(numclusters, alpha = 1) # Generate color palette: any RColorBrewer or viridis palette
            
            color.palette <- viridis(n = numMetadata, 
                                     alpha = 1,
                                     option = COLOR_SCALE_ToUse) # Generate color palette: any RColorBrewer or viridis palette
            
            colors.final[-CLUSTERS.TO.COLOR] <- "gray"
            }
        } else {
          numclusters <- length(unique(plot.df.final$Cluster)) # Calculate number of clusters to color
          if (ColorPaletteToUse == "Rainbow") {
            colors.final <- rainbow(numclusters, alpha = 1) # Generate color palette: any RColorBrewer or viridis palette
            } else if (ColorPaletteToUse == "viridis") {
              colors.final <- viridis(n = numMetadata, 
                                      alpha = 1,
                                      option = COLOR_SCALE_ToUse) # Generate color palette: any RColorBrewer or viridis palette
              } else if (ColorPaletteToUse == "RColorBrewer") {
                colors.final <- colorRampPalette(color = brewer.pal(n = as.numeric(COLOR_SCALE_ToUse[1]),
                                                                    name = COLOR_SCALE_ToUse[2]))(length(unique(plot.df.final[["Cluster"]])))
                }
          }
  p <- plot_ly(plot.df.final, type = 'scatter3d', mode = 'markers', 
               x = ~umap_x, y = ~umap_y, z = ~umap_z, 
               marker = list(size = 0.3),
               alpha = 1, 
               color = ~factor(Cluster, levels = sort(unique(Cluster))), 
               colors = colors.final,
               hoverinfo = 'text', text = ~paste('</br> Cell ID: ', cell.id,
                                                 '</br> Cluster ID: ', Cluster,
                                                 '</br> Treatment: ', Treatment,
                                                 '</br> TimePoint: ', TimePoint
                                                 )
               ) %>%
    layout(
      scene = list(
        camera = list( # Adjust values to alter camera view: https://plotly.com/python/3d-camera-controls/
          up = list(x = 0, y = 0, z = 1), # Default: (x = 0, y = 0, z = 1)
          center = list(x = 0, y = 0, z = 0), # Default: (x = 0, y = 0, z = 0)
          eye = list(x = -0.85, y = -1.65, z = 2) # Default: (x = 1.25, y = 1.25, z = 1.25)
          )
        ),
      legend = list(itemsizing = 'constant',      # Prevent scaling with data marker size
                         font = list(size = 12),       # Adjust font size if needed
                         itemwidth = 30,               # Increase legend spacing
                         bgcolor = 'white',            # Optional: improve contrast
                         bordercolor = 'black',
                         borderwidth = 1
                         )
           )
  
  } else if (PLOT.COLOR.FACTOR == "Treatment") {
    # Define a color palette to use for the PLOT.COLOR.FACTOR
      numMetadata <- length(unique(plot.df.final[[PLOT.COLOR.FACTOR]]))
      if (ColorPaletteToUse == "Rainbow") {
        colors.final <- rainbow(numclusters, alpha = 1) # Generate color palette: any RColorBrewer or viridis palette
        
      } else if (ColorPaletteToUse == "viridis") {
        colors.final <- viridis(n = numMetadata, 
                                alpha = 1,
                                option = COLOR_SCALE_ToUse) # Generate color palette: any RColorBrewer or viridis palette
      } else if (ColorPaletteToUse == "RColorBrewer") {
        colors.final <- colorRampPalette(color = brewer.pal(n = as.numeric(COLOR_SCALE_ToUse[1]),
                                                            name = COLOR_SCALE_ToUse[2]))(length(unique(plot.df.final[["Cluster"]])))
      }
      # color.palette <- viridis(n = numMetadata, 
      #                          alpha = 1,
      #                          option = COLOR_SCALE_ToUse) # Generate color palette: any RColorBrewer or viridis palette
    # Change PLOT.COLOR.FACTOR column into a factor ordered by METADATA_ORDER
      plot.df.final[[PLOT.COLOR.FACTOR]] <- factor(plot.df.final[[PLOT.COLOR.FACTOR]],
                                                   levels = METADATA_ORDER)
    p <- plot_ly(plot.df.final, type = 'scatter3d', mode = 'markers', 
                 x = ~umap_x, y = ~umap_y, z = ~umap_z, 
                 marker = list(size = 0.5), # Incr size from "Cluster" to improve visibility of Treatments
                 alpha = 1, 
                 color = ~factor(Treatment, levels = sort(unique(Treatment))), 
                 colors = colors.final,
                 hoverinfo = 'text', text = ~paste('</br> Cell ID: ', cell.id,
                                                   '</br> Cluster ID: ', Cluster,
                                                   '</br> Treatment: ', Treatment,
                                                   '</br> TimePoint: ', TimePoint
                                                   )) %>%
      layout(
        scene = list(
          camera = list( # Adjust values to alter camera view: https://plotly.com/python/3d-camera-controls/
            up = list(x = 0, y = 0, z = 1), # Default: (x = 0, y = 0, z = 1)
            center = list(x = 0, y = 0, z = 0), # Default: (x = 0, y = 0, z = 0)
            eye = list(x = -0.85, y = -1.65, z = 2) # Default: (x = 1.25, y = 1.25, z = 1.25)
            )
          ),
        legend = list(itemsizing = 'constant',      # Prevent scaling with data marker size
                      font = list(size = 12),       # Adjust font size if needed
                      itemwidth = 30,               # Increase legend spacing
                      bgcolor = 'white',            # Optional: improve contrast
                      bordercolor = 'black',
                      borderwidth = 1
                      )
        )
    
    } else if (PLOT.COLOR.FACTOR == "TimePoint") {
      # Define a color palette to use for the PLOT.COLOR.FACTOR
        numMetadata <- length(unique(plot.df.final[[PLOT.COLOR.FACTOR]]))
        if (ColorPaletteToUse == "Rainbow") {
          colors.final <- rainbow(numclusters, alpha = 1) # Generate color palette: any RColorBrewer or viridis palette
          
        } else if (ColorPaletteToUse == "viridis") {
          colors.final <- viridis(n = numMetadata, 
                                  alpha = 1,
                                  option = COLOR_SCALE_ToUse) # Generate color palette: any RColorBrewer or viridis palette
        } else if (ColorPaletteToUse == "RColorBrewer") {
          colors.final <- colorRampPalette(color = brewer.pal(n = as.numeric(COLOR_SCALE_ToUse[1]),
                                                              name = COLOR_SCALE_ToUse[2]))(length(unique(plot.df.final[["Cluster"]])))
        }
        # color.palette <- viridis(n = numMetadata, 
        #                          alpha = 1,
        #                          option = COLOR_SCALE_ToUse) # Generate color palette: any RColorBrewer or viridis palette
      # Change PLOT.COLOR.FACTOR column into a factor ordered by METADATA_ORDER
        plot.df.final[[PLOT.COLOR.FACTOR]] <- factor(plot.df.final[[PLOT.COLOR.FACTOR]],
                                                     levels = METADATA_ORDER)
      p <- plot_ly(plot.df.final, type = 'scatter3d', mode = 'markers', 
                   x = ~umap_x, y = ~umap_y, z = ~umap_z, 
                   marker = list(size = 0.5), # Incr size from "Cluster" to improve visibility of timepoints
                   alpha = 1, 
                   color = ~factor(TimePoint, levels = sort(unique(TimePoint))), 
                   colors = colors.final,
                   hoverinfo = 'text', text = ~paste('</br> Cell ID: ', cell.id,
                                                     '</br> Cluster ID: ', Cluster,
                                                     '</br> Treatment: ', Treatment,
                                                     '</br> TimePoint: ', TimePoint
                                                     )) %>%
        layout(
          scene = list(
            camera = list( # Adjust values to alter camera view: https://plotly.com/python/3d-camera-controls/
              up = list(x = 0, y = 0, z = 1), # Default: (x = 0, y = 0, z = 1)
              center = list(x = 0, y = 0, z = 0), # Default: (x = 0, y = 0, z = 0)
              eye = list(x = -0.85, y = -1.65, z = 2) # Default: (x = 1.25, y = 1.25, z = 1.25)
              )
            ),
          legend = list(itemsizing = 'constant',      # Prevent scaling with data marker size
                        font = list(size = 12),       # Adjust font size if needed
                        itemwidth = 30,               # Increase legend spacing
                        bgcolor = 'white',            # Optional: improve contrast
                        bordercolor = 'black',
                        borderwidth = 1
                        )
          )
      }  
 

## Visualize 3D plot ===========================================================
# Visualize in R using plotly
# Note: Here is where to add/customize what is viewed in the plot 
if (plot_in_viewer) {
  p
}


### Save 3D plot outputs =======================================================
## Save as static image ========================================================
save_image(p,
           file = paste0(UMAP.FILENAME,
                         "_By", ifelse(is.null(PLOT.COLOR.FACTOR), "Cluster", PLOT.COLOR.FACTOR),
                         ifelse(test = CLUSTER.SINGLE.COLOR, 
                                yes = paste0("_", paste0(CLUSTERS.TO.COLOR, collapse = ","), "-", SINGLE.COLOR), 
                                no = ""),
                         ifelse(test = CLUSTER.SINGLE.COLOR,
                                yes = "",
                                no = paste0("_", ColorPaletteToUse)),
                         ".", OUTPUT.TYPE),
           width = OUTPUT.Width, height = OUTPUT.Height, scale = 3)

## Save as interactive UMAP ====================================================
if (isTRUE(SAVE.HTML)) {
  p <- p %>% # update browser-based configurations
    config(plot_ly(), toImageButtonOptions = list(format = OUTPUT.TYPE, 
                                                  filename = paste0(OUTPUT.BASENAME, 
                                                                    "_By", 
                                                                    ifelse(is.null(PLOT.COLOR.FACTOR), 
                                                                           "Cluster", 
                                                                           PLOT.COLOR.FACTOR)), 
                                                  height= 1200, width= 1600, scale = 3))
  
  p.bundle <- partial_bundle(p) # Significantly reduces the output file size
  saveWidget(widget = p.bundle, 
             file = paste0(UMAP.FILENAME,
                           "_By", ifelse(is.null(PLOT.COLOR.FACTOR), "Cluster", PLOT.COLOR.FACTOR),
                           ifelse(CLUSTER.SINGLE.COLOR, paste0("_", paste0(CLUSTERS.TO.COLOR, collapse = ","), "-red"), ""),
                           ".html"),
             selfcontained = T)
  # Note: If want single HTML file (e.g. simple output to share), use selfcontained = T
  # If want to generate multiple HTML files for single data library, consider using
  # selfcontained = F to save single "lib" for all .html outputs (way more efficient!)
  # e.g. saveWidget(p.bundle, HTML.FILENAME, selfcontained = F, libdir = "lib")
}

setwd(INPUT.FOLDER)

print("Finish View_3D_UMAP.R")
