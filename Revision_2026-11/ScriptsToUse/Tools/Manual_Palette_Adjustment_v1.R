################################################################################
# Manual_Palette_Adjustment.R
################################################################################
# Author:   Jonathon Sewell
# Created:  2026-10-31
# Edited:   NA
# Version:  1
# Update:   NA
#
# Purpose:
# * Manually test color palettes to confirm readibility in manuscript plots
# * Update Palette_Color.csv in the Working Directory and Helper Directory
#   > Updated file then updates the color palette used in every figure



#### Define the parameters of which their colors will be adjusted ##############
## From Palette_Color define the Set and ID that will be adjsuted
Set_to_Modify <- "Cluster"
ID_to_Modify <- seq(20)

## Define the Palette Directories were the csv will be read and modified
pallette_file_WD <- "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/CellsToAnalyze_157Gated/2026-10-29 Working Directory/A. Cell ID Clusters/Palette_Colors.csv"

pallette_file_Helper <- "/Users/jmsewell/Documents/BIMS Stuff/Lab Work/Experiments/Mass Cytometry_Temp/2025-05-13 BDNF_Rescue_FullSampleSet/2026-09-15 Figures_Revisions_v1/ScriptsToUse/Helpers/Palette_Colors.csv"

## Published palette
my_colors_pub <- rainbow(20, alpha = 1)



#### REFERENCE: Which combos are in the manuscript that need comparison? #######
## Cluster groups on line graphs:
ClustersInFigures <- list(Fig3C = c(1,2,4,7,13), # Problems: 1,4,7
                          Fig3D = c(1,3,15,17), # Problems: 15,17
                          Fig3E = c(9,17,18), # Problems: NA
                          Fig5F = c(1,15), # Problems: NA
                          Fig5F = c(2,4,7,13), # Problems: 4,7
                          Fig6D = c(1,2,4,18), # Problems: 1,4 (18)
                          Fig6E = c(1,2,4,18) # Problems: 1,4 (18)
                          )

# Clusters that are on line graphs some point:
# [1]  1  2  3  4  7  9 13 15 17 18 (10 clusters)
sort(unique(unlist(ClustersInFigures)))

# Clusters that are on multiple line graphs some point:
# [1]  1  2  4  7 13 15 17 18
sort(unique(unlist(ClustersInFigures)[
  duplicated(unlist(ClustersInFigures))]))

# Clusters that are on only one line graph set:
# [1] 3 9
sort(unique(unlist(ClustersInFigures))[!unique(unlist(ClustersInFigures)) %in%
                                         unique(unlist(ClustersInFigures)[
                                           duplicated(unlist(ClustersInFigures))])])



#### Manually test different color combinations ################################
## REFERENCE: Viewable published color palette
published_palette <- list(
  "1" =  "#FF0000FF",
  "2" =  "#FF4D00FF",
  "3" =  "#FF9900FF",
  "4" =  "#FFE500FF",
  "5" =  "#CCFF00FF",
  "6" =  "#80FF00FF",
  "7" =  "#33FF00FF",
  "8" =  "#00FF19FF",
  "9" =  "#00FF66FF",
  "10" = "#00FFB2FF",
  "11" = "#00FFFFFF",
  "12" = "#00B3FFFF",
  "13" = "#0066FFFF",
  "14" = "#001AFFFF",
  "15" = "#3300FFFF",
  "16" = "#7F00FFFF",
  "17" = "#CC00FFFF",
  "18" = "#FF00E6FF",
  "19" = "#FF0099FF",
  "20" = "#FF004DFF"
)

## Confirm that published_palette is the same as my_colors_pub
if (!identical(my_colors_pub, unname(unlist(published_palette)))) {
  stop("Check that my_colors_pub and published_palette are the same")
}

## REFERENCE: Colors that look ok on white background
list(
  "1" =  "#FF0000FF",
  "2" =  "#FF4D00FF",
  "3" =  "#FF9900FF",
  # "4" =  "#FFE500FF",
  # "5" =  "#CCFF00FF",
  # "6" =  "#80FF00FF",
  # "7" =  "#33FF00FF",
  "8" =  "#00FF19FF",
  "9" =  "#00FF66FF",
  # "10" = "#00FFB2FF",
  # "11" = "#00FFFFFF",
  "12" = "#00B3FFFF",
  "13" = "#0066FFFF",
  "14" = "#001AFFFF",
  "15" = "#3300FFFF",
  "16" = "#7F00FFFF",
  "17" = "#CC00FFFF",
  "18" = "#FF00E6FF",
  "19" = "#FF0099FF",
  "20" = "#FF004DFF"
  )

## MANUAL ADJUSTMENT: Modify the proposed palette
# Tip:  Look at the UMAP to determine which colors to switch around
#       Update palette, run UMAP, and iterate until satisfied with colors
proposed_palette <- list(
  "1" =  "#FF9900FF",
  "2" =  "#3300FFFF",
  "3" =  "#FF4D00FF",
  "4" =  "#00B3FFFF",
  "5" =  "#00FFFFFF",
  "6" =  "#FFE500FF",
  "7" =  "#FF0000FF",
  "8" =  "#00FF19FF",
  "9" =  "#FF00E6FF",
  "10" = "#80FF00FF",
  "11" = "#00FF66FF",
  "12" = "#00FFB2FF",
  "13" = "#CC00FFFF",
  "14" = "#33FF00FF",
  "15" = "#0066FFFF",
  "16" = "#001AFFFF",
  "17" = "#7F00FFFF",
  "18" = "#FF004DFF",
  "19" = "#FF0099FF",
  "20" = "#CCFF00FF"
)

## Check if proposed_palette is missing any colors from published_palette
if (length(unlist(published_palette)[!unlist(published_palette) %in% unlist(proposed_palette)]) > 0) {
  stop(paste0("proposed_palette is missing:\n", unlist(published_palette)[!unlist(published_palette) %in% unlist(proposed_palette)]))
}



#### Visualize color palette ###################################################
## Define colors to visualize
my_colors <- unlist(proposed_palette)
names(my_colors) <- names(proposed_palette)


## Visualize the defined combinations
FigureComparisons <- list()
for (i in names(ClustersInFigures)) {
  
  barplot_input <- my_colors[names(my_colors) %in% ClustersInFigures[[i]]]
  
  FigureComparisons[[i]] <- ggplot(data = as.data.frame(barplot_input),
                                   mapping = aes(x = .data[["barplot_input"]], 
                                                 fill = .data[["barplot_input"]])) +
    geom_bar() + 
    labs(x = i,
         y = "") +
    theme_void() + 
    labs(title = i) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  }

gridExtra::grid.arrange(grobs = FigureComparisons, 
                       nrow = length(FigureComparisons) / 2)



#### Make adjustments to palette file based on manual color choices ############
## Read in palette file to adjust
pal.in.adjust <- read.csv(pallette_file_WD, stringsAsFactors = FALSE, colClasses = "character")


## Make adjustments to the data frame
pal.in.adjust[pal.in.adjust$Set %in% Set_to_Modify & pal.in.adjust$ID %in% ID_to_Modify, "Color"] <- my_colors


## Write adjusted palette file to Palette Directories
# Write to Working Directory
write.csv(x = pal.in.adjust,
          file = pallette_file_WD)

# Write to Helper directory
write.csv(x = pal.in.adjust,
          file = pallette_file_Helper)


