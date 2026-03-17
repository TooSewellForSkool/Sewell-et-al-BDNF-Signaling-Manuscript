
# Updated script:
# Jonathon Sewell
# Version: 3
# 2025-09-12

# Updates:
# (1) Enabled re-scaling of color
# (2) Enabled selection of color scale
# (3) Enabled output as static image

# Source script: Figure_7A_View_3D_UMAP.R and "XY_Plot_By_Marker.R"
# Amy Van Deusen, University of Virginia and Corey Williams, University of Virginia
# 2022 and 2019
# Visualize 3D UMAP layout with plotly
# See https://plotly.com/r/3d-scatter-plots/#basic-3d-scatter-plot 
# This will allow for .html and other outputs to be interactive for sample metadata

# Updates to Source script:
# (1) Enables re-scaling selected marker values from 0-1 to better highlight marker expression and generating a scale comparable across markers
# (2) Color palette no longer discrete (switched to continuous "viridis")

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
  DownSample <- TRUE # TRUE / FALSE
  DownSample_Percent <- 90

## Define whether or not to re-scale marker values from 0-1
  NormalizeScale <- TRUE # TRUE (scale from 0-1) / FALSE (unique scale per marker)
  ShiftScaleBy <- 50 # NULL / Select a percent to shift what is the max value (helps negate skew of scale by "hot" cells)
  ColorAbove <- "cyan"
  
## Define the color scale to be used
  COLOR_SCALE_ToUse <- "A" # A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo

## Define which identified clusters to plot and which layout to use
  CLUSTER_ROUND <- 1
  LAYOUT_ROUND <- 1
  LAYOUT_GROUP <- 1 # Don't use group 0, not assigned XY coordinates
  
## Identify markers for which to generate individual UMAPs displaying their expression values
  MARKERS.TO.PLOT <- c(#"Tuj1_89Y", 
                       # "CXC3CR1_141Pr",
                       # "Vimentin_143Nd",
                       # "HB9_148Nd",
                       # "NeuN_150Nd",
                       # "Sox2_160Gd",
                       # "MAP2_164Dy",
                       # "BLBP_172Yb",
                       # "s100b_146Nd",
                       "OligoO4_145Nd"
                       ) # c("TrkB_i_174Yb") / "All" will loop through each input marker

  # All markers
  # "Tuj1_89Y" / "pPLK1_113In" / "H3K9ac_115In" / "pGSK3b_139La" / "CXC3CR1_141Pr" / "p-p38_142Nd" / 
  # "Vimentin_143Nd" / "pPLCG-2_144Nd" / "OligoO4_145Nd" / "s100b_146Nd" / "pSTAT5_147Sm" / "HB9_148Nd" /  
  # "CGRP_149Sm" / "NeuN_150Nd" / "p-p53_151Eu" / "CHAT_152Sm" / "pSTAT3_153Eu" / "pSrc_154Sm" / 
  # "nNOS_155Gd" / "TrkB_s_158Gd" / "pATM_159Tb" / "Sox2_160Gd" / "p75NTR_161Dy" / "p-cJUN_162Dy" / 
  # "DCX_163Dy" / "MAP2_164Dy" / "yH2AX_165Ho" / "pNFKB_166Er" / "Islet1_167Er" / "pAKt_168Er" / 
  # "GFAP_169Tm" / "pTrkB_170Er" / "pERK_171Yb" / "BLBP_172Yb" / "CC3_173Yb" / "TrkB_i_174Yb" / 
  # "pS6_175Lu" / "pCreb_176Lu" / "DNA_191Ir" / "DNA_193Ir"
  
## Metadata choices
# Choose which, if any, metadata varibles to color by (NULL = Color plot by "Cluster" as default)
  PLOT.COLOR.FACTOR <- NULL # Can also choose metadata (column IDs) - ex. "Cluster" / "Treatment" / "TimePoint"
# List of sample metadata to include in plot dataframe (added for each cell)
  METADATA.TO.INCLUDE <- c("Expt", "Treatment", "TimePoint") 
# How would you like to color metadata?
  METADATA.COLORS <- NULL # Can add c("color1", "color2", etc.) or NULL will use default palette (e.g. viridis(n))

## Make choices on how to plot data
# Has the inputted dataframe already been subsampled prior to running this script?
  SUBSAMPLE.IDS.FILENAME <- NULL # NULL if no subsample IDs ("Subsample_ID_by_Cluster.csv")
# Is there already a file with color information?
  CLUSTER.COLORS <- NULL # .csv file ("colors_R1.csv") or NULL will use default colors

## Input file information
  INPUT.FOLDER <- getwd()
  OUTOUT.FOLDER <- INPUT.FOLDER
  UMAP.3D.LAYOUT <- sub("X", CLUSTER_ROUND,"umap_xyz_3D_RX.csv")
  EXPRESSION.FILENAME.1 <- "expression_matrix_analysis.csv" # "Concat_Transformed_by_Cluster.csv"
  EXPRESSION.FILENAME.2 <- "expression_matrix_other.csv"
  CLUSTERS.FILENAME <- "cluster_R1_assigns.csv" # "Clusters_Subsampled.csv"
  METADATA.FILENAME <- "metadata.csv"
  FILENUMS.FILENAME <- "filenums.csv" # "filenums_Subsampled.csv"
  PANEL.FILENAME <- "panel.csv"

## Output file information
  OUTPUT.BASENAME <- "UMAP_3D_"
  OUTPUT.TYPE <- "pdf" # Choices are "svg" (preferred), "pdf", "png", etc.
  OUTPUT.Height <- 1200
  OUTPUT.Width <- 1600
  
  SAVE.HTML <- FALSE # TRUE = output .html file, FALSE = do not output file
  
  FOLDER_BASENAME <- "Marker_Expression_VAR_PERpercSubset_XYZ"
  if (LAYOUT_ROUND > 1) {
    output_dir <- paste0("05.2-", # Script ID
                         format(Sys.time(), "%Y-%m-%d"), "_", # Date info
                         sub("VAR_PER", paste0("L", LAYOUT_ROUND,
                                               "_G", LAYOUT_GROUP,
                                               "_", DownSample_Percent),
                             FOLDER_BASENAME),
                         if (NormalizeScale) {
                           "_0-1scale"
                           },
                         if (!is.null(ShiftScaleBy)) {
                           paste0("_ShiftScale", ShiftScaleBy)
                         }
                         )
  } else {
    output_dir <- paste0("05.2-", # Script ID
                         format(Sys.time(), "%Y-%m-%d"), "_", # Date info
                         sub("VAR_PER", paste0("L", LAYOUT_ROUND,
                                               "_", DownSample_Percent),
                             FOLDER_BASENAME),
                         if (NormalizeScale) {
                           "_0-1scale"
                           },
                         if (!is.null(ShiftScaleBy)) {
                           paste0("_ShiftScale", ShiftScaleBy)
                         }
                         )
  }

### Read necessary files =======================================================
expression.filename.1 <- paste0(INPUT.FOLDER, "/", EXPRESSION.FILENAME.1)
expression.filename.2 <- paste0(INPUT.FOLDER, "/", EXPRESSION.FILENAME.2)
exprs.in <- cbind(fread(expression.filename.1, header=TRUE, sep=",", check.names=FALSE),
                  fread(expression.filename.2, header=TRUE, sep=",", check.names=FALSE))
layout.filename <- paste0(INPUT.FOLDER, "/", UMAP.3D.LAYOUT)
layout.in <- read.table(layout.filename, header=TRUE, sep=",", check.names=FALSE)
clusters.filename <- paste0(INPUT.FOLDER, "/", CLUSTERS.FILENAME)
clusters.in <- read.table(clusters.filename, header=TRUE, sep=",", check.names=FALSE)
metadata.filename <- paste0(INPUT.FOLDER, "/", METADATA.FILENAME)
metadata.in <- read.table(metadata.filename, header=TRUE, sep=",", check.names=FALSE)
panel.filename <- paste0(INPUT.FOLDER, "/", PANEL.FILENAME)
panel.in <- read.table(panel.filename, header=TRUE, sep=",", check.names=FALSE)
filenums.filename <- paste0(INPUT.FOLDER, "/", FILENUMS.FILENAME)
filenums.in <- read.table(filenums.filename, header=TRUE, sep=",", check.names=FALSE)
if (!is.null(SUBSAMPLE.IDS.FILENAME)) {
  subsamples.filename <- paste0(INPUT.FOLDER, "/", SUBSAMPLE.IDS.FILENAME)
  subsamples.in <- read.table(subsamples.filename, header=FALSE, sep=",", check.names=FALSE)
}

### Prepare plotting dataframe =================================================
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
plot.vars <- panel.in[which(panel.in[,"Plot"]==1), "Fixed_Param"] # Identify markers to plot
plot.df.final.indexed <- cbind(plot.df, metadata.final) # Add metadata to plotting dataframe (quick)

plot.df.exprs <- cbind(plot.df.final.indexed, # Generate larger plotting dataframe with expression data
                       exprs.in[ , ..plot.vars])

# If user opts to rescale each marker from 0-1, apply range01 function
if (NormalizeScale) {
  # Set up and apply a function to normalize on 0-1 scale
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}

  plot.df.exprs[ , -c(1:which(colnames(plot.df.exprs) == "TimePoint"))] <- apply(plot.df.exprs[ , -c(1:which(colnames(plot.df.exprs) == "TimePoint"))],2,range01)

}

# Subset (as necessary) data to help plot datasets that are too large for plotly
set.seed(123) # Make downsampling reproducible
if (DownSample) {
  # Downsample
  plot.df.final <- plot.df.exprs %>%
    group_by(File) %>%
    slice_sample(prop = (DownSample_Percent/100)) %>%
    ungroup()

  # Create the removed rows dataframe
  plot.df.final.downsampled.removed <- anti_join(plot.df.exprs, plot.df.final, by = names(plot.df.exprs))

} else {
  plot.df.final <- plot.df.exprs
}


### Generate plotly plot =======================================================
dir.create(output_dir)
setwd(output_dir)

if ("All" %in% MARKERS.TO.PLOT) {
  MARKERS.TO.PLOT <- colnames(plot.df.final)[which(!colnames(plot.df.final) %in% 
                                                     colnames(plot.df.final)[c(1:which(colnames(plot.df.exprs) == "TimePoint"))])]
  }

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


## Loop through each marker and save the resulting image
for (i in MARKERS.TO.PLOT) {
  
  CurrentMarkerIndex <- which(plot.vars == i)

  # numclusters <- max(plot.df.final[ , plot.vars[[CurrentMarkerIndex]]]) # Calculate number of clusters to color
  # colors.final <- hcl.colors(numclusters, palette = "viridis", alpha = 1) # Generate color palette
  
  # Generate viridis colors
  viridis_colors <- viridis(100,
                            option = COLOR_SCALE_ToUse) # A = magma / B = inferno / C = plasma / D = viridis
  
  # Create a custom color scale
  # Map viridis from 0 to threshold, then user-defined color from threshold to 1
  color_scale <- lapply(0:99, function(i) {
    val <- i / 99 * (ShiftScaleBy / 100)
    list(val, viridis_colors[i + 1])
  })
  
  # Add user-defined color at threshold and above
  color_scale <- append(color_scale, list(
    list((ShiftScaleBy / 100), ColorAbove),
    list(1, ColorAbove)
  ))
  
  p <- plot_ly(plot.df.final, type = 'scatter3d', mode = 'markers', 
               x = ~umap_x, y = ~umap_y, z = ~umap_z, 
               marker = list(size = 0.3,
                             color = plot.df.final[[plot.vars[CurrentMarkerIndex]]],
                             colorscale = color_scale,
                             colorbar = list(
                               title = "Expression\nLevel",
                               tickvals = c(0, (ShiftScaleBy / 100), 1)
                               #,
                               #ticktext = c("Low", "Medium", "High", "Above Threshold"),
                               #len = 0.5,
                               #thickness = 20
                             )
                             
               ),
               hoverinfo = 'text', text = ~paste('</br> Cell ID: ', cell.id,
                                                 '</br> Cluster ID: ', Cluster,
                                                 '</br> Treatment: ', Treatment,
                                                 '</br> TimePoint: ', TimePoint,
                                                 '</br> Value: ', plot.df.final[[plot.vars[CurrentMarkerIndex]]]
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
    ) %>%
    config(plot_ly(), toImageButtonOptions = list(format = OUTPUT.TYPE, 
                                                  filename = paste0(OUTPUT.BASENAME, MARKERS.TO.PLOT), 
                                                  height= 1200, width= 1600, scale = 3))
  
  UMAP.FILENAME <- if (DownSample) { # Naming for downsampling
    paste0("UMAP_3D_", DownSample_Percent, "percSubset", 
           "_L", LAYOUT_ROUND, "_C", CLUSTER_ROUND, 
           "_", plot.vars[CurrentMarkerIndex],
           if (!is.null(ShiftScaleBy)) {
             paste0("_ShiftScale", ShiftScaleBy)
           }
           )
  } else {
    paste0("UMAP_3D_All", 
           "_L", LAYOUT_ROUND, "_C", CLUSTER_ROUND, 
           "_", plot.vars[CurrentMarkerIndex],
           if (!is.null(ShiftScaleBy)) {
             paste0("_ShiftScale", ShiftScaleBy)
           }
           )
  }
  
  paste0("UMAP_3D_", plot.vars[CurrentMarkerIndex])
  
  ### Save 3D plot outputs =======================================================
  ## Save as static image ========================================================
  save_image(p,
             file = paste0(UMAP.FILENAME, ".", OUTPUT.TYPE),
             width = OUTPUT.Width, height = OUTPUT.Height, scale = 3)
  
  
  ## Save as interactive UMAP ====================================================
  if (isTRUE(SAVE.HTML)) {
    p.bundle <- partial_bundle(p) # Significantly reduces the output file size
    saveWidget(widget = p.bundle, 
               file = paste0(UMAP.FILENAME, ".html"),
               selfcontained = T)
    # Note: If want single HTML file (e.g. simple output to share), use selfcontained = T
    # If want to generate multiple HTML files for single data library, consider using
    # selfcontained = F to save single "lib" for all .html outputs (way more efficient!)
    # e.g. saveWidget(p.bundle, HTML.FILENAME, selfcontained = F, libdir = "lib")
  }

}

setwd(INPUT.FOLDER)

print("Finish View_3D_UMAP.R")
