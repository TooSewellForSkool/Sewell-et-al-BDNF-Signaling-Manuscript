
# Author: Jonathon Sewell
# Date: 2025-09-18
# Version: 1
# Update:
# (1) Incorporated ability to filter on File_Nums from previous clustering rounds to generate file

# Original Script: 01_Get_Metadata.R
#smgoggin
#05.14.19

print("Start 01_Get_Metadata.R")

rm(list = ls())
.libPaths( c("/project/zunderlab/R/4.0.3_Pipeline2", .libPaths()) )

library(flowCore)
library(data.table)

print("Libraries loaded")


Original_WD <- getwd() # setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) to set source file location as working directory

# Define the Cluster round that these files will be used for
INTENDED_CLUSTER_ROUND <- 2

##Set filepaths ===================================================================================
# Get path to fcs files
args = commandArgs(trailingOnly=TRUE) # Get command line arguments
if (length(args) == 0) {
  FCS_FILEPATH <- getwd() # If there's no command line argument, set to current directory
} else {
  FCS_FILEPATH <- args[1] # Use filepath specified by command line argument
}
# Get path to where you want to save csv files
if (length(args) < 2) {
  WRITE_PATH <- FCS_FILEPATH # If there's no command line argument, make this the same as FCS_FILEPATH
} else {
  WRITE_PATH <- args[2] # Use filepath specified by command line argument
}
METADATA_FILENAME <- "metadata.csv"
PANEL_FILENAME <- "panel.csv"

##read in fcs files =============================================================================== 
file_list <- list.files(FCS_FILEPATH, full.names = TRUE, pattern = "*.fcs")

if (INTENDED_CLUSTER_ROUND > 1) {
  # Read in file numbers associated with the expression file
  FILE_NUM <- paste0("filenums", ".csv")
  FileNum <- read.table(FILE_NUM, header=TRUE, sep=",",
                        check.names=FALSE)
  
  file_list <- file_list[unique(FileNum[[1]])]
  }

fs <- read.flowSet(file_list, transformation=FALSE)

##read in metadata ================================================================================ 
#file names
metadata_df <- as.data.frame(cbind(sampleNames(fs), 1:length(sampleNames(fs))))
colnames(metadata_df) <- c("File_Name", "File_Num")

#NOTE: THIS ASSUMES YOU HAVE THE SAME PANEL FOR ALL FCS FILES!!!!!!!
panel_name <- pData(parameters(fs[[1]]))$name
panel_desc <- pData(parameters(fs[[1]]))$desc
# Fix parameter names (remove spaces, hyphens, NAs, etc.)
fixed_params <- c()
for (n in 1:ncol(fs[[1]])) {
  pname <- pData(parameters(fs[[1]]))[n, "name"]
  pdesc <- pData(parameters(fs[[1]]))[n, "desc"]
  if (is.na(pname) || is.na(pdesc)) {
    if (is.na(pname)) {
      fparam <- pdesc
    } else {
      fparam <- pname
    }
  } else {
    fparam <- pdesc
  }
  if (grepl("^\\d+\\w+_", fparam)){
    mname <- regmatches(fparam, regexpr("^[^_]+", fparam))
    fparam <- paste0(sub("^[^_]+_", "", fparam), "_", mname)
  }
  fparam <- gsub("[ /]", "_", fparam)
  fixed_params <- c(fixed_params, fparam)
}
panel_fixed <- fixed_params

panel_analysis <- rep(0, length(panel_name)) #fill with zeroes to start
panel_plot <- rep(0, length(panel_name)) #fill with zeroes to start

panel_trans <- rep("Linear", length(panel_name))
panel_trans[grep("Di$|bc_neg|bc_pos", panel_name)] <- "Asinh"

panel_factor <- rep(1, length(panel_name)) #fill with 1 to start
panel_factor[grep("Di$|bc_neg|bc_pos", panel_name)] <- 5

panel_df <- as.data.frame(cbind(unlist(panel_name), unlist(panel_desc),
                                unlist(panel_fixed), panel_analysis,
                                panel_plot, panel_trans, panel_factor))
colnames(panel_df) <- c("Name", "Desc", "Fixed_Param", "Analysis", "Plot",
                        "Transform", "Factor")

##write metadata files ============================================================================
if (!dir.exists(WRITE_PATH)) {
  dir.create(WRITE_PATH)
}
fwrite(metadata_df, file=paste(WRITE_PATH, METADATA_FILENAME, sep='/'), row.names=FALSE)
fwrite(panel_df, file=paste(WRITE_PATH, PANEL_FILENAME, sep='/'), row.names=FALSE)

print("Metadata files output complete")
print("Finish 01_Get_Metadata.R")
