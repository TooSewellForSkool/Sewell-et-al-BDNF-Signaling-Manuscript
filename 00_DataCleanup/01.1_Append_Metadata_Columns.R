# 2025-06-17
# Jonathon Sewell
# Determine events per FCS file and append to metadata.csv with additional columns containing relevant information

#Adapted from:
#Corey Williams, University of Virginia
#17 May, 2019
#Subsample and make csv file containing all points asinh transformed

print("Start 01.1_Append_Event_Count.R")

rm(list = ls())
.libPaths( c("/project/zunderlab/R/4.0.3_Pipeline2", .libPaths()) )

library(flowCore)
#library(ZunderPipelineFunctions)
#source("./00_Pipeline_Helper_Functions.R")
library(data.table)

print("Libraries loaded")

## Input parameters ===============================================================================

METADATA_FILENAME <- "metadata.csv"
PANEL_FILENAME <- "panel.csv"
#SUBSAMPLE_IDXS_FILENAME <- "subsample_idxs.csv" #unused if no subsampling
#EXPRS_MAT_ANALYSIS_FILENAME <- "expression_matrix_analysis.csv"
#EXPRS_MAT_OTHER_FILENAME <- "expression_matrix_other.csv"
#FILENUMS_FILENAME <- "filenums.csv"
#ASINH_FACTOR <- NULL #set to NULL if set in panel
#SUBSAMPLE_BOOL <- TRUE
#NUM_SUBSAMPLES <- 10000 #change this value if you want subsampling. Number of cells per file

# Get path to fcs files
args = commandArgs(trailingOnly=TRUE) # Get command line arguments
if (length(args) == 0) {
  INPUT_FOLDER <- getwd() # If there's no command line argument, set to current directory
} else {
  INPUT_FOLDER <- args[1] # Use filepath specified by command line argument
}
# Get path to where you want to save csv files
if (length(args) < 2) {
  OUTPUT_FOLDER <- INPUT_FOLDER # If there's no command line argument, make this the same as FCS_FILEPATH
} else {
  OUTPUT_FOLDER <- args[2] # Use filepath specified by command line argument
}

print("Got input parameters, loading fcs files")

## Read metadata ==================================================================================
metadata_path <- paste(INPUT_FOLDER, METADATA_FILENAME, sep="/")
metadata <- read.table(metadata_path, header=TRUE, sep=",", check.names=FALSE)

## Read fcs files =================================================================================
fs <- read.flowSet(metadata$File_Name, path=INPUT_FOLDER, transformation=FALSE,
                   truncate_max_range=FALSE)

## Determine event counts for each fcs file =======================================================
ColumnsToAdd <- c()
for (i in 1:nrow(metadata)) {
  SimpleSampleIdentifier <- strsplit(x = sub(".*normalized_(.*?)_Cells.*", "\\1",
                                             metadata$File_Name[i]),
                                     split = "_")[[1]]
  
  ColumnsToAdd <-rbind(ColumnsToAdd, 
                       data.frame(File_Name = metadata$File_Name[i],
                                  Reordered = "User-Defined",
                                  BarcodeSet = as.numeric(sub(".*Sample([0-9]+).*", "\\1", metadata$File_Name[i])),
                                  SampleCollection = "User-Defined",
                                  Expt = "User-Defined",
                                  Condition = if (length(SimpleSampleIdentifier) == 3) {
                                    paste0(SimpleSampleIdentifier[1], "_", SimpleSampleIdentifier[2])
                                    } else if (length(SimpleSampleIdentifier) > 3) {
                                      gsub(pattern = "\\.", replacement = "",
                                           x = paste0(SimpleSampleIdentifier[1], SimpleSampleIdentifier[2], "_", SimpleSampleIdentifier[3]))
                                      } else if (length(SimpleSampleIdentifier) < 3) {
                                        gsub(pattern = "\\.", replacement = "",
                                             x = paste0(SimpleSampleIdentifier[1], SimpleSampleIdentifier[2]))
                                        },
                                  Treatment = if (length(SimpleSampleIdentifier) <= 3) {
                                    paste0(SimpleSampleIdentifier[1])
                                    } else if (length(SimpleSampleIdentifier) > 3) {
                                      gsub(pattern = "\\.", replacement = "", 
                                           x = paste0(SimpleSampleIdentifier[1], SimpleSampleIdentifier[2]))
                                      },
                                  TimePoint = if (grepl("min", paste0(SimpleSampleIdentifier, collapse = ""))) {
                                    round(x = as.numeric(gsub(pattern = "min", replacement = "", 
                                                              x = SimpleSampleIdentifier[grepl("min", SimpleSampleIdentifier)])) / 60,
                                          digits = 3)
                                    } else if (grepl("hr", paste0(SimpleSampleIdentifier, collapse = ""))) {
                                      as.numeric(gsub(pattern = "hr", replacement = "", 
                                                      x = SimpleSampleIdentifier[grepl("hr", SimpleSampleIdentifier)]))
                                      } else {
                                        NA
                                      },
                                  Replicate = if (length(SimpleSampleIdentifier) > 4) {
                                    if (grepl("min", paste0(SimpleSampleIdentifier, collapse = ""))) {
                                      as.numeric(sub(".*min_(.*?)_.*", "\\1", paste0(SimpleSampleIdentifier, collapse = "_")))
                                      } else if (grepl("hr", paste0(SimpleSampleIdentifier, collapse = ""))) {
                                        as.numeric(sub(".*hr_(.*?)_.*", "\\1", paste0(SimpleSampleIdentifier, collapse = "_")))
                                        }
                                    } else if (grepl("Stock", paste0(SimpleSampleIdentifier, collapse = ""))) {
                                      as.numeric(sub(".*Stock_(.*?)_.*", "\\1", paste0(SimpleSampleIdentifier, collapse = "_")))
                                      } else {
                                      as.numeric(last(SimpleSampleIdentifier))
                                      },
                                  Event_Count = nrow(fs[[i]])
                                  )
                       )
  }






## Append ColumnsToAdd to metadata, then save the updated metadata.csv =============================
fwrite(merge(x = metadata, 
             y = ColumnsToAdd, 
             by = "File_Name"),
       file = METADATA_FILENAME, 
       row.names=FALSE)


print(paste0("Minimum Event Count = ", min(ColumnsToAdd$Event_Count)))
print("Finish 01.1_Append_Event_Count.R")


