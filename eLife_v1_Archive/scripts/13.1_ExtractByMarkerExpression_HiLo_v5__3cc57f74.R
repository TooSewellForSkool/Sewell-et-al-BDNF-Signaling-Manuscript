
# Author: Jonathon Sewell
# Date: 2025-01-02
# Version: 5
# Update:
# (1) Added DelineateByQuantile
# (2) Cleaned up certain variables





# Source script:
#Austin Keeler, adapated from Corey Williams, University of Virginia
#01 August, 2020
#Plot cells over threshold MARKER EXPRESSION

# 2025-06-11 (13.1_v1)
# Edited by Jonathon Sewell
# Enables greater flexibility with user-defined variables
# Output files have better compatibility with downstream Pipeline2 scripts



print("Start PlotByThreshold.R")

rm(list = ls())
.libPaths( c( .libPaths(), "~/local/R_libs/") )

library(ggfortify)
library(data.table)
library(dplyr)
library(ggplot2)
# library(stringr)



INPUT.FOLDER <- getwd() # setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) to set source file location as working directory



#### User-defined input parameters #############################################
## Define values related to the selected marker(s)
# Marker(s) to threshold on: Example: "Tuj1_89Y", "CXC3CR1_141Pr", "s100b_146Nd", "OligoO4_145Nd", "TrkB_i_174Yb", "p75NTR_161Dy", "TrkB_s_158Gd"
    # SELECTED.MARKER <- c("p75NTR_161Dy", "TrkB_i_174Yb") # c("p75NTR_161Dy", "TrkB_i_174Yb")
    SELECTED.MARKER <- c(x = "TrkB_i_174Yb",
                         y = NA # NA to exclude from analysis
                         )
    # p75NTR_161Dy   / TrkB_i_174Yb
    # Sox2_160Gd     / BLBP_172Yb
    # Tuj1_89Y       / CXC3CR1_141Pr
    # Vimentin_143Nd


## Define how Threshold Categorizations (ex Lo, Med, Hi) will be defined
# User-defined threshold values:
    EXPRESSION.DELINEATION.VALUES <- c(x = list(c(0.15, 1.4, 3.7)),
                                       y = list(c()) # from QuantileValues <- c(1/3, 2/3)
                                       )
    # p75NTR_161Dy = c(0.45, 1.45, 3.95)  /   TrkB_i_174Yb = c(0.15, 1.4, 3.7)
    # Sox2_160Gd = 2.45   /   BLBP_172Yb = 1.95
    # Tuj1_89Y = 3.15     /   CXC3CR1_141Pr = 2.35
    # Vimentin_143Nd = 0.5
    # These values can be defined by looking at expression histograms/violin plots, or by using alternative justification

# Threshold values defined by quantiles:
  DelineateByQuantile <- FALSE # TRUE / FALSE
  QuantileValues <- c() # c(1/3, 2/3) / c(0.5)

# Operator(s): ">=" or "<=" to indicate "greater than" or "less than" the corresponding EXPRESSION.DELINEATION.VALUES
    # THRESHOLD.OPERATOR <- c(">=", ">=")
    THRESHOLD.OPERATOR <- c(x = list(c(">=", ">=", ">=")), # c(">=", ">=") / c(">=", ">=", ">=")
                            y = list(c(">=", ">=", ">=")) # c(">=", ">=") / c(">=", ">=", ">=")
                            )
    
# Have you justified your EXPRESSION.DELINEATION.VALUES values?
    DELINEATION.VALUES.justified <- TRUE # TRUE / FALSE

    
## Define cluster-related variables
    CLUSTER.ROUND <- 1 # Need to change this for the level of subsetting/subclustering!
    CLUSTERS.MERGED <- FALSE

## Define variables related to plotting parameters
    PlotPercentages <- TRUE # TRUE / FALSE
    # OUTPUT_HEIGHT <- 7
    # OUTPUT_WIDTH <- 7
    # Plot_Dimensions <- c(4 , 5) # Width x Height (in.)
    BIAXIAL_OUTPUT_FORMAT <- "png" # "png" or "pdf"
    HISTOGRAM_OUTPUT_FORMAT <- "pdf" # "png" or "pdf"
    PERCENT_OUTPUT_FORMAT <- "pdf" # "png" or "pdf"
    
#### Defined parameters ########################################################
## Input basenames
CLUSTERS.BASENAME <- "/cluster_RX_assigns.csv"
LAYOUT.BASENAME <- "/umap_xy_RX.csv"

## Input filenames
CLUSTERS.FILENAME <- sub("RX", paste0("R", CLUSTER.ROUND), CLUSTERS.BASENAME) #Cluster file name adjusted for round of clustering
if (CLUSTERS.MERGED) { #Adjust cluster file name if clusters were merged
  CLUSTERS.FILENAME <- sub(".csv", "_merged.csv", CLUSTERS.FILENAME)
  }
LAYOUT.FILENAME <- sub("RX", paste0("R", CLUSTER.ROUND), LAYOUT.BASENAME) # UMAP_layout.csv
PANEL.FILENAME <- "/panel.csv"
METADATA.FILENAME <- "/metadata.csv"
FILENUMS.FILENAME <- "/filenums.csv"
EXPRS.ANALYSIS.FILENAME <- "/expression_matrix_analysis.csv"
EXPRS.OTHER.FILENAME <- "/expression_matrix_other.csv"




print("Input parameters loaded")



#### Read in necessary files ###################################################
## Validate lengths of user-defined variables match 
# (i.e. each SELECTED.MARKER has a corresponding EXPRESSION.DELINEATION.VALUES and THRESHOLD.OPERATOR)
# if (!(length(SELECTED.MARKER) == length(unique(substr(names(EXPRESSION.DELINEATION.VALUES), 0, 1))) &
#      length(EXPRESSION.DELINEATION.VALUES) == length(THRESHOLD.OPERATOR))) {
#   stop("Check that SELECTED.MARKERS, EXPRESSION.DELINEATION.VALUES, and THRESHOLD.OPERATORS are appropriately matched.")
# }
## Read in files
layout.in <- read.csv(paste0(INPUT.FOLDER, LAYOUT.FILENAME))
clusters.in <- read.csv(paste0(INPUT.FOLDER, CLUSTERS.FILENAME))
panel.in <- read.csv(paste0(INPUT.FOLDER, PANEL.FILENAME))
metadata.in <- read.csv(paste0(INPUT.FOLDER, METADATA.FILENAME))
filenums.in <- read.csv(paste0(INPUT.FOLDER, FILENUMS.FILENAME))
exprs.analysis.in <- fread(paste0(INPUT.FOLDER, EXPRS.ANALYSIS.FILENAME), stringsAsFactors = FALSE, data.table = FALSE)
exprs.other.in <- fread(paste0(INPUT.FOLDER, EXPRS.OTHER.FILENAME), stringsAsFactors = FALSE, data.table = FALSE)

print("Necessary files read in")



#### Threshold data based on thresholding parameters ###########################
## Combine exprs_mat data frames to prep thresholding
# thresholding.ct.df <- merge(layout.clusters.id, layout.concat.transform.id)
thresholding.ct.df <- cbind(exprs.analysis.in, 
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

OUTPUT.FOLDER.NAME <- paste0("13.1-", # Script ID
                             format(Sys.time(), "%Y-%m-%d"), # Date info
                             "_ExtractByMarkerExpression_", 
                             ALL.THRESHOLD.NAME)


# Generate a histogram to depict thresholding
HISTOGRAM_LIST <- list()
for (i in names(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))])){
MAX.SELECTED.MARKER <- max(thresholding.ct.df[[SELECTED.MARKER[[i]]]], na.rm = TRUE)

HISTOGRAM <- thresholding.ct.df %>% 
  ggplot(aes(x = .data[[SELECTED.MARKER[[i]]]])) +
  geom_histogram(aes(y = after_stat(density)),
                 fill = "lightgray",
                 bins = 175) +
  scale_x_continuous(breaks = seq(0, MAX.SELECTED.MARKER, by = 1)) +
  geom_density(lwd = 1.2,
               linetype = 2,
               colour = "magenta") +
  geom_vline(xintercept = EXPRESSION.DELINEATION.VALUES[[i]],
             linetype = 'longdash')

HISTOGRAM_LIST[[i]] <- HISTOGRAM

print(HISTOGRAM_LIST[[i]])
}


if (length(names(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))])) == 2) {
  BIAXIAL_PLOT <- thresholding.ct.df %>% 
    ggplot(aes(x = .data[[SELECTED.MARKER[["x"]]]],
               y = .data[[SELECTED.MARKER[["y"]]]])) + 
    geom_point(size = 0.01, alpha = 0.75, color = "black") + 
    geom_density_2d(#aes(color = after_stat(level)),
                    bins = 10,
                    color = "white"
                    ) +
    geom_vline(xintercept = EXPRESSION.DELINEATION.VALUES[["x"]],
               linetype = 'longdash', color = "magenta") +
    geom_hline(yintercept = EXPRESSION.DELINEATION.VALUES[["y"]],
               linetype = 'longdash', color = "cyan")
  }


## Generate a ThresholdCategory column to categorize each event
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
                        probs = QuantileValues) %>% unname(), # delination value(s)
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


## Paste all category values across rows into a new column called ThresholdCategory
thresholding.ct.df.category <- thresholding.ct.df.category %>% 
  # create column to denote co-expression levels
  mutate(ThresholdCategory = if (length(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))]) == 1) {
    thresholding.ct.df.category[ , setdiff(colnames(thresholding.ct.df.category), colnames(thresholding.ct.df))] %>% as.character()
    } else {
      apply(X = thresholding.ct.df.category[ , setdiff(colnames(thresholding.ct.df.category), colnames(thresholding.ct.df))], 
            MARGIN = 1, FUN = paste0, collapse = "")
      }
    ) %>% 
  mutate(ThresholdCategory = factor(ThresholdCategory,
                                    levels = sort(unique(ThresholdCategory))))


## Generate plot for percentage of each ThresholdCategory
if (PlotPercentages) {
  Percentages_Category <- thresholding.ct.df.category %>% 
    group_by(ThresholdCategory) %>% 
    summarize(Percent = (n() / nrow(thresholding.ct.df.category) * 100),
              .groups = "drop")
  
  PercPlot <- ggplot(data = Percentages_Category, 
                     aes(x = ThresholdCategory, y = Percent)) + 
    geom_segment( aes(x = ThresholdCategory, xend = ThresholdCategory, 
                      y = 0, yend = Percent)) +
    geom_point(aes(color = ThresholdCategory), 
               size = 3) + 
    geom_text(aes(label = round(Percent, digits = 2)), 
              nudge_y = 0,
              nudge_x = 0.325,
              size = 4)
  }


## Update Histograms with ThresholdCategory values
# HISTOGRAM_LIST

for (i in names(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))])) {
  Percentages_Category_ByMarker <- thresholding.ct.df.category %>% 
    group_by(thresholding.ct.df.category[ , paste0("ThreshCat_", SELECTED.MARKER[[i]])]) %>% 
    summarize(Percent = (n() / nrow(thresholding.ct.df.category) * 100),
              .groups = "drop") %>% 
    rename(ThreshCat_ByMarker = 'thresholding.ct.df.category[, paste0("ThreshCat_", SELECTED.MARKER[[i]])]')
  
  HISTOGRAM_LIST[[i]] <- HISTOGRAM_LIST[[i]] + 
    annotate(geom = "text",
             x = c(0, EXPRESSION.DELINEATION.VALUES[[i]]),
             y = max(density(thresholding.ct.df[[SELECTED.MARKER[[i]]]])[["y"]], na.rm = TRUE) + 
               (max(density(thresholding.ct.df[[SELECTED.MARKER[[i]]]])[["y"]], na.rm = TRUE)*0.1) * 
               1:length(Percentages_Category_ByMarker[["ThreshCat_ByMarker"]]),
             label = paste0(Percentages_Category_ByMarker$ThreshCat_ByMarker, "=\n" , round(Percentages_Category_ByMarker$Percent, digits = 1), "%")
             )
  }


## Generate a biaxial plot to depict thresholding
# EXPRESSION.DELINEATION.VALUES[grepl(i, names(EXPRESSION.DELINEATION.VALUES))]
if (length(names(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))])) == 2 |
    mean(table(SELECTED.MARKER)) %% 2 == 0) {
  BIAXIAL_PLOT <- thresholding.ct.df.category %>% 
    ggplot(aes(x = .data[[SELECTED.MARKER[["x"]]]],
               y = .data[[SELECTED.MARKER[["y"]]]])) + 
    geom_point(aes(color = ThresholdCategory), 
               size = 0.01, alpha = 0.75) + 
    geom_density_2d(bins = 10,
                    color = "white"
                    ) +
    
    scale_x_continuous(breaks = seq(0, max(thresholding.ct.df[[SELECTED.MARKER[["x"]]]], na.rm = TRUE), by = 1)) +
    scale_y_continuous(breaks = seq(0, max(thresholding.ct.df[[SELECTED.MARKER[["y"]]]], na.rm = TRUE), by = 1)) +
    
    geom_vline(xintercept = EXPRESSION.DELINEATION.VALUES[["x"]],
               linetype = 'longdash', color = "black") +
    geom_hline(yintercept = EXPRESSION.DELINEATION.VALUES[["y"]],
               linetype = 'longdash', color = "black") +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  
  if (PlotPercentages) {
    # Define min/max values for each axis
    MaxValue_y <- ceiling(max(thresholding.ct.df.category[[SELECTED.MARKER[["y"]]]], na.rm = TRUE))
    MinValue_y <- min(thresholding.ct.df.category[[SELECTED.MARKER["y"]]], na.rm = TRUE) + (MaxValue_y * -0.025)
    
    MaxValue_x <- ceiling(max(thresholding.ct.df.category[[SELECTED.MARKER[["x"]]]], na.rm = TRUE))
    MinValue_x <- min(thresholding.ct.df.category[[SELECTED.MARKER["x"]]], na.rm = TRUE) + (MaxValue_y * -0.025)
    
    # Define Percentage label coordinates
    label_x_coordinates <- c("11" = min(EXPRESSION.DELINEATION.VALUES[["x"]], na.rm = TRUE)*0.9, 
                             "12" = min(EXPRESSION.DELINEATION.VALUES[["x"]], na.rm = TRUE)*0.9,
                             "13" = min(EXPRESSION.DELINEATION.VALUES[["x"]], na.rm = TRUE)*0.9, 
                             "14" = min(EXPRESSION.DELINEATION.VALUES[["x"]], na.rm = TRUE)*0.9,
                             "21" = MaxValue_x*0.25, "22" = MaxValue_x*0.25, "23" = MaxValue_x*0.25, "24" = MaxValue_x*0.25,
                             "31" = MaxValue_x*0.7, "32" = MaxValue_x*0.7, "33" = MaxValue_x*0.7, "34" = MaxValue_x*0.7,
                             "41" = MaxValue_x*0.9, "42" = MaxValue_x*0.9, "43" = MaxValue_x*0.9, "44" = MaxValue_x*0.9
                             )
    
    label_y_coordinates <- c("11" = 0, "12" = MaxValue_y*0.25, "13" = MaxValue_y*0.75, "14" = MaxValue_y*0.95,
                             "21" = 0, "22" = MaxValue_y*0.25, "23" = MaxValue_y*0.75, "24" = MaxValue_y*0.95,
                             "31" = 0, "32" = MaxValue_y*0.25, "33" = MaxValue_y*0.75, "34" = MaxValue_y*0.95,
                             "41" = 0, "42" = MaxValue_y*0.25, "43" = MaxValue_y*0.75, "44" = MaxValue_y*0.95
                             )
    
    BIAXIAL_PLOT <- BIAXIAL_PLOT +
      annotate(geom = "text",
               x = unname(label_x_coordinates[which(names(label_x_coordinates) %in% Percentages_Category$ThresholdCategory)]),
               y = unname(label_y_coordinates[which(names(label_y_coordinates) %in% Percentages_Category$ThresholdCategory)]),
               label = paste0(Percentages_Category$ThresholdCategory, "=\n" , round(Percentages_Category$Percent, digits = 1), "%")
               )
    }
  }


print("Threshold has been applied")



#### Write output files ########################################################
## Create and set ouput directory
dir.create(OUTPUT.FOLDER.NAME) #Create output directory
OUTPUT.FOLDER <- paste0(INPUT.FOLDER, "/", OUTPUT.FOLDER.NAME)
setwd(OUTPUT.FOLDER) #Set directory to output folder


## Output histogram to depict thresholding
if (DELINEATION.VALUES.justified) {
  for (i in 1:length(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))])) {
    if (HISTOGRAM_OUTPUT_FORMAT == "png") {
      png(paste0("Histogram_Thresholds_", SELECTED.MARKER[i], ".", HISTOGRAM_OUTPUT_FORMAT))
      } else if (HISTOGRAM_OUTPUT_FORMAT == "pdf") {
        pdf(paste0("Histogram_Thresholds_", SELECTED.MARKER[i], ".", HISTOGRAM_OUTPUT_FORMAT))
        }
    print(HISTOGRAM_LIST[[i]])
    dev.off()
    }
  
  if (length(SELECTED.MARKER[which(!is.na(SELECTED.MARKER))]) == 2) {
    if (BIAXIAL_OUTPUT_FORMAT == "png") {
      png(paste0("Biaxial_Thresholds_", SELECTED.MARKER[1], "-vs-", SELECTED.MARKER[2], ".", BIAXIAL_OUTPUT_FORMAT))
      } else if (BIAXIAL_OUTPUT_FORMAT == "pdf") {
        pdf(paste0("Biaxial_Thresholds_", SELECTED.MARKER[1], "-vs-", SELECTED.MARKER[2], ".", BIAXIAL_OUTPUT_FORMAT))
        }
    print(BIAXIAL_PLOT)
    dev.off()
    }
  
  if (PlotPercentages) {
    if (PERCENT_OUTPUT_FORMAT == "png") {
      png(paste0("Percentages_ThresholdCategories", ".", PERCENT_OUTPUT_FORMAT))
      } else if (PERCENT_OUTPUT_FORMAT == "pdf") {
        pdf(paste0("Percentages_ThresholdCategories", ".", PERCENT_OUTPUT_FORMAT))
        }
    print(PercPlot)
    dev.off()
    }
  }


## Output new files that will (or may) be useful for downstream analysis of thresholded data
fwrite(as.data.frame(thresholding.ct.df),
       file = "expression_matrix_all.csv", 
       row.names = FALSE)

fwrite(as.data.frame(thresholding.ct.df.category[which(thresholding.ct.df.category$ThresholdCategory == "lo"), ]),
       file = "expression_matrix_all_lo.csv", 
       row.names = FALSE)

fwrite(as.data.frame(thresholding.ct.df.category[which(thresholding.ct.df.category$ThresholdCategory == "hi"), ]),
       file = "expression_matrix_all_hi.csv", 
       row.names = FALSE)

fwrite(as.data.frame(thresholding.ct.df.category[ , "ThresholdCategory"]),
       file = "ThresholdCategory.csv", 
       row.names = FALSE)

## Output original files into OUTPUT.FOLDER so they are available for downstream analysis
fwrite(exprs.analysis.in,
       file = sub("/", "", EXPRS.ANALYSIS.FILENAME), 
       row.names = FALSE)

fwrite(exprs.other.in,
       file = sub("/", "", EXPRS.OTHER.FILENAME), 
       row.names = FALSE)

fwrite(filenums.in,
       file = sub("/", "", FILENUMS.FILENAME), #Remove slash at beginning of FILENUMS.FILENAME
       row.names = FALSE)

fwrite(metadata.in,
       file = sub("/", "", METADATA.FILENAME), # Remove slash at beginning of METADATA.FILENAME
       row.names = FALSE)

fwrite(panel.in,
       file = sub("/", "", PANEL.FILENAME), # Remove slash at beginning of METADATA.FILENAME
       row.names = FALSE)

fwrite(layout.in,
       file = sub("/", "", sub("/", "", LAYOUT.FILENAME)), # Remove slash at beginning of LAYOUT.FILENAME
       row.names = FALSE)

fwrite(clusters.in,
       file = sub("/", "", sub("/", "", CLUSTERS.FILENAME)), # Remove slash at beginning of CLUSTERS.FILENAME
       row.names = FALSE)



#### Return to original working directory ######################################
setwd(INPUT.FOLDER)







