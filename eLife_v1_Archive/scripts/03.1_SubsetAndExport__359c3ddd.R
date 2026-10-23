
# Author: Jonathon Sewell
# Date: 2025-09-15
# Version: 1
# PURPOSE: 
# • Subset dataset based on user-defined parameters and export relevant .csv files into new folder
# • Specific use to generate files that are filtered based on R round of clustering for use in R+1 clustering that is not influenced by R clustering input


Original_WD <- getwd() # setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) to set source file location as working directory



#### Variables to Change ######################################################
VariableToSubset <- "Condition" # "File_Num" /  "Expt" /  "Condition" / "Treatment" / "TimePoint" / "Cluster" / "Group"
ValuesToSubset <- c("NoTreat_0hr",
                    "ContDepr_15min", "ContDepr_16hr", "ContDepr_1hr", "ContDepr_4hr", "ContDepr_5min", 
                    # "ContDepr_72hr",
                    "CompMed_15min", "CompMed_16hr", "CompMed_1hr", "CompMed_4hr", "CompMed_5min", 
                    # "CompMed_72hr",
                    "BDNF_15min", "BDNF_16hr", "BDNF_1hr", "BDNF_4hr", "BDNF_5min", 
                    # "BDNF_72hr",
                    "BDNFTest_1hr", "BDNF+K252aTest_1hr", "ContDepr_Test" #,
                    # "DIV4Stock_1", "DIV4Stock_2", "DIV4Stock_3", "DIV4Stock_4", "DIV4Stock_5", "DIV4Stock_6", "DIV4Stock_9",
                    # "Universal1", "Universal2", "Universal3", "Universal4"
                    )

# For Expt:
# c(1, 2, 3, 4 ,5, 6, 7)

# For Condition
# c("NoTreat_0hr", # Expt1
#   "ContDepr_15min", "ContDepr_16hr", "ContDepr_1hr", "ContDepr_4hr", "ContDepr_5min", "ContDepr_72hr", # Expt2
#   "CompMed_15min", "CompMed_16hr", "CompMed_1hr", "CompMed_4hr", "CompMed_5min", "CompMed_72hr", # Expt3
#   "BDNF_15min", "BDNF_16hr", "BDNF_1hr", "BDNF_4hr", "BDNF_5min", "BDNF_72hr", # Expt4
#   "BDNFTest_1hr", "BDNF+K252aTest_1hr", "ContDepr_Test", # Expt5
#   "DIV4Stock_1", "DIV4Stock_2", "DIV4Stock_3", "DIV4Stock_4", "DIV4Stock_5", "DIV4Stock_6", "DIV4Stock_9", # Expt6
#   "Universal1", "Universal2", "Universal3", "Universal4") # Expt7

# For Treatment:
# c("BDNF", "CompMed", "ContDepr", "NoTreat", "Universal", "DIV4Stock", "BDNFTest", "BDNF+K252aTest")

# For TimePoint:
# c(0.000, 0.083, 0.250, 1.000, 4.000, 16.000, 72.000, NA)



# For Clusters:
  # c(2, 4, 11, 1, 14, 15, 5, # Path A
  #   12, 10, 6, # Path B
  #   19, 21, # Astrocytic (Glia)
  #   8, 18, # Microglial (Glia)
  #   22, # Oligodendrocytic (Glia)
  #   16, 20, # Other Progenitor-like (Misc.)
  #   17, # Misc. Neuron (Misc.)
  #   7, 13, # Misc. Glia (Glia)
  #   3, 9 # Unknown
  #   )

# For Groups:
  # c(1, 2, 3, 4, 5)
  # 1 = Path A / 2 = Path B / 3 = Misc. / 4 = Glia / 5 = Unknown


CLUSTER_ROUND <- NULL # Need to change this for the level of subsetting/subclustering!
CLUSTERS_MERGED <- FALSE


#### Variables that may need to be changed #####################################
FILE_NUM <- paste0("filenums", ".csv")
PANEL_INFILE <- paste0("panel", ".csv")
METADATA_FILENAME <- paste0("metadata", ".csv")

EXPRS_MAT_INFILE1 <- "expression_matrix_analysis.csv"
EXPRS_MAT_INFILE2 <- "expression_matrix_other.csv"

CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)

GROUPS_BASENAME <- "cluster_RX_groups.csv"
groups_filename <- sub("RX", paste0("R", CLUSTER_ROUND), GROUPS_BASENAME)



#### Read in necessary files ##################################################
# Read in file numbers associated with the expression file
FileNum <- read.table(FILE_NUM, header=TRUE, sep=",",
                      check.names=FALSE)

# Read in panel (will show which expression markers will be included)
panel <- read.table(PANEL_INFILE, header=TRUE, sep=",",
                    check.names=FALSE)

# Read in metadata
metadata <- read.table(METADATA_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)

# Read in cluster file
if (!is.null(CLUSTER_ROUND)) {
  clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                            check.names=FALSE)
  }

# Read in group file
if (!is.null(CLUSTER_ROUND)) {
  groups_in <- read.table(groups_filename, header=TRUE, sep=",",
                          check.names=FALSE)
  }

# Read in expression matrix and perform initial clean up of data
exprs_mat_IN1 <- fread(paste0(Original_WD, "/", EXPRS_MAT_INFILE1), 
                       stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)

exprs_mat_IN2 <- fread(paste0(Original_WD, "/", EXPRS_MAT_INFILE2), 
                       stringsAsFactors = FALSE, header=TRUE, data.table = FALSE)


#### Define variables for subsetting based on user-defined inputs ##############
if (VariableToSubset %in% colnames(metadata)) {
  IndexFromDataframeToSubset <- metadata[[VariableToSubset]]
  
  SubsetOnValues <- metadata$File_Num[which(IndexFromDataframeToSubset %in% ValuesToSubset)]
  
  SubsetIndices <- which(FileNum$V1 %in% SubsetOnValues)
  
  } else if (VariableToSubset == "Cluster") {
    IndexFromDataframeToSubset <- clusters_in$R1
    
    SubsetOnValues <- ValuesToSubset
    
    SubsetIndices <- which(IndexFromDataframeToSubset %in% SubsetOnValues)
    
    } else if (VariableToSubset == "Group") {
      IndexFromDataframeToSubset <- groups_in$Group_Assigns
  
      SubsetOnValues <- which(IndexFromDataframeToSubset %in% ValuesToSubset)
      
      SubsetIndices <- which(clusters_in$R1 %in% SubsetOnValues)
      
      }


#### Create folder for output files ############################################
output_dir_All <- paste0("03_", # Script ID
                         format(Sys.time(), "%Y-%m-%d"), # Date info
                         "_SubsetOn-",
                         VariableToSubset,
                         if (!VariableToSubset == "Condition") {
                           paste0(ValuesToSubset, collapse = ",")
                         }
                         )
dir.create(output_dir_All)
setwd(paste0(Original_WD, "/", output_dir_All))



#### Write and save output files ###############################################
# Write panel file into output directory
  write.csv(x = panel, 
            file = PANEL_INFILE,
            row.names = FALSE)

# Write filtered metadata file into output directory
write.csv(x = as.data.frame(metadata[which(metadata[[VariableToSubset]] %in% ValuesToSubset), ]), 
          file = METADATA_FILENAME,
          row.names = FALSE)

# Write filtered FileNum file into output directory
  write.csv(x = as.data.frame(FileNum[SubsetIndices, ]) %>% 
              rename("V1" = "FileNum[SubsetIndices, ]"), 
            file = FILE_NUM,
            row.names = FALSE)

# If cluster data is read in, write relevant files
if (!is.null(CLUSTER_ROUND)) {
  # Write filtered Clusters (input round) file into output directory
  write.csv(x = as.data.frame(clusters_in[SubsetIndices , ]) %>% 
              rename(!!paste0("R", CLUSTER_ROUND) := !!sym("clusters_in[SubsetIndices, ]")), 
            file = clusters_filename,
            row.names = FALSE)
  
  # Write groups file into output directory
  write.csv(x = groups_in, 
            file = groups_filename,
            row.names = FALSE)
  }

# Write filtered expression_matrix_analysis file into output directory
  write.csv(x = as.data.frame(exprs_mat_IN1[SubsetIndices, ]), 
            file = EXPRS_MAT_INFILE1,
            row.names = FALSE)

# Write filtered expression_matrix_other file into output directory
  write.csv(x = as.data.frame(exprs_mat_IN2[SubsetIndices, ]), 
            file = EXPRS_MAT_INFILE2,
            row.names = FALSE)



#### Return to original directory once everything is completed ################
setwd(Original_WD)


print("Complete! Have fun with your new spreadsheets, ya nerd ;)")




