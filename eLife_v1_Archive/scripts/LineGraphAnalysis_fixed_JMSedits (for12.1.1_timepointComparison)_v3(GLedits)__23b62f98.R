# Title: LineGraphAnalysis
# Version: 11 (all bugs fixed)
# Author: Jonathon Sewell (edited for robustness)
# Date: 2025-11-10
# Version 1 - Set up from LineGraphAnalysis and updated to include stats by Grace Lee (LineGraphAnalysis_fixed);
#     JMS Updated for flexibility
# Version 2 - JMS updated to output stats into a table
# Version 3 - JMS updated minor bugs in columning naming convention

# Update

library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(cowplot)
library(viridis)
library(gridExtra)
library(grid)
# install.packages("multcomp")
library(multcomp)
library(car) # For Levene's Test
# install.packages("FSA")
library(FSA) # For Dunn's Test

Original_WD <- getwd()

#### User-defined Variables ####################################################
### Statistical comparison settings ===========================================
RunStatistics <- TRUE   # TRUE / FALSE
List_p_values <- TRUE   # TRUE / FALSE
stats_yAxis_scaler <- 1 # 1 / 2 / 5

StatisticalComparison_Type  <- "TimePoint"   # "Condition", "TimePoint", "Cluster"
StatisticalComparison_Value <- 0             # e.g. 0 for 0hr
Statistical_Test <- "ANOVA"                 # "t.test" / "wilcox" / "ANOVA" / "DUNN" / "auto"
Alpha <- 0.05

### Marker-related variables ==================================================
# Are different markers plotted on the same graph?
PlotBy_marker <- FALSE # TRUE/FALSE

# Which markers are going to be included?
MarkerSelect <- "Individual" # "Individual" / "All"
if (PlotBy_marker) MarkerSelect <- "Individual"

IndividualMarkers <- c("pERK_171Yb", "pAKt_168Er", "pCreb_176Lu", "pS6_175Lu", "pGSK3b_139La",
                       "pPLK1_113In", "yH2AX_165Ho",#"TrkB_s_158Gd", 
                       "pATM_159Tb", "pPLCG-2_144Nd",
                       "pSrc_154Sm", "H3K9ac_115In", "p-p38_142Nd", "CC3_173Yb", "pSTAT3_153Eu",
                       "pSTAT5_147Sm", "p-p53_151Eu", "p-cJUN_162Dy", "pNFKB_166Er"#, "pTrkB_170Er"
                       )
  # c(
  # "pERK_171Yb",
  # "pAKt_168Er",
  # "pCreb_176Lu",
  # "pS6_175Lu",
  # "pGSK3b_139La"#,
  # # "p-cJUN_162Dy"
  # )

# Are there defined categories based on a marker threhold?
IsThereThresholdCategory <- FALSE # TRUE / FALSE
WhichToPlot <- c()

### Cluster-related variables =================================================
AnalyzeAllClusters <- FALSE # TRUE / FALSE (if selecting specific clusters or if using pseudobulk data)
ClustersToAnalyze <- c() # c(1, 2, 4, 7, 13) / c(9, 17, 18) / c(1, 3, 15, 17) / c(1, 2, 4, 18) / c()

if (AnalyzeAllClusters) {
  ClusterName <- "AnalyzeAllClusters"
} else if (!AnalyzeAllClusters & length(ClustersToAnalyze) == 0) {
  ClusterName <- "PseudoBulk"
} else {
  ClusterName <- paste(ClustersToAnalyze, collapse = ",")
}

### Treatments ================================================================
TreatmentToAnalyze <- c("NoTreat", "ContDepr", "BDNF") # "NoTreat" / "ContDepr" / "CompMed" / "BDNF"
TreatmentName <- paste(TreatmentToAnalyze, collapse = ",")

### Define PlotBy =============================================================
if (PlotBy_marker) {
  if (length(IndividualMarkers) > 1) {
    PlotBy <- "Marker"
  } else {
    stop("Select more IndividualMarkers to compare")
  }
} else if (IsThereThresholdCategory) {
  if (length(WhichToPlot) == 1 & isFALSE(AnalyzeAllClusters) & length(ClustersToAnalyze) > 0) {
    PlotBy <- "Cluster"
  } else {
    PlotBy <- "ThreshCat" # "ThreshCat" 
  }
} else if (length(TreatmentToAnalyze) > 3) { 
  if (length(ClustersToAnalyze) > 1) {
    stop("Only select one ClustersToAnalyze at a time until Jonathon can adjust this script :P")
  } else {
    PlotBy <- "Treatment"
  }
} else if (length(ClustersToAnalyze) == 0) {
  PlotBy <- "Treatment"
} else {
  PlotBy <- "Cluster"
}

### Time handling =============================================================
TimePoint_VariableType <- "Categorical" # "Categorical" / "Continuous"
EarlyTimePoint_Cutoff <- 5

### Output info ===============================================================
PlotAsGrid <- FALSE
plot_column_length <- 5

UserDefinedColor <- FALSE # TRUE / FALSE
if (PlotBy == "Cluster") {
  ColorToUse <- "Rainbow"
} else if (length(TreatmentToAnalyze) > 3 | UserDefinedColor) { 
  # user-defined colors for each treatment to keep them consistent
  ColorToUse <- c(NoTreat = "lightblue4", 
                  ContDepr = "#30123BFF", 
                  CompMed = "#3BAF29FF", 
                  BDNF = "#B63679FF")
  # Select only those that are included in TreatmentToAnalyze
  ColorToUse <- ColorToUse[names(ColorToUse) %in% TreatmentToAnalyze]
  } else {
  ColorToUse <- "A" # With viridis, A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo
}

Plot_Width <- 14
Plot_Height <- 10
Output_File_Type <- "pdf" # "pdf", "png", etc.

#### Read in and prep data #####################################################
input_file_transformation <- "_FC" # "_FC" / "hr_RAWvalues"
input_file_type <- "\\.csv$"

input_filenames <- list.files(Original_WD, pattern = input_file_type) %>%
  (\(x) x[grepl(input_file_transformation, x)])()

if (length(input_filenames) == 0) {
  stop(
    paste0(
      "No input files found.\n",
      "Working directory: ", Original_WD, "\n",
      "Expected: .csv files containing '", input_file_transformation, "' in the filename.\n",
      "Example match: Something_FC-<Category>.csv"
    )
  )
}

read_input_file <- function(path) {
  ext <- tools::file_ext(path)
  switch(ext,
         csv = read_csv(path, show_col_types = FALSE),
         tsv = read_tsv(path, show_col_types = FALSE),
         txt = read_table(path),
         stop("Unsupported file type: ", ext)
  )
}

## Read in file
input_dataframe_in <- read_input_file(file.path(Original_WD, 
                                                input_filenames[which(!grepl("RepsAvg", input_filenames))]) # Exclude file containing "RepsAvg"
                                      ) %>% 
  dplyr::select(-'...1') # Remove first row (contains row names)


if (is.null(input_dataframe_in) || nrow(input_dataframe_in) == 0) {
  stop("Read-in failed: input_dataframe_in is empty. Check your file contents/format.")
}


## Add a "Treatment" column to help with subsetting and make sure column naming is consistent
input_dataframe <- input_dataframe_in %>% 
  rename("Condition" = names(input_dataframe_in)[grep("Condition" , names(input_dataframe_in))],
         "TimePoint" = names(input_dataframe_in)[grep("TimePoint" , names(input_dataframe_in))]) %>% # Rename columns for consistency in script
  mutate("Treatment" = sub(pattern = "_[^_]+$", # Everything after last "_"
                           replacement = "", # Remove pattern
                           x = .data[["Condition"]]), # data from "Condition" column
         .after = "Condition") # Mutated column after "Condition"

## Add a "ThresholdCat" column if not already present
if (is.null(input_dataframe[["ThreshCat"]])) {
  input_dataframe <- input_dataframe %>% # Add a ThreshCat column for consistency in script
    mutate("ThreshCat" = NA, # Add NA values, since they were not provided in the input
           .after = "Cluster") # Mutated column after "Condition"
}

## Marker selection
if (MarkerSelect == "All") {
  SelectMarker <- colnames(input_dataframe)[!colnames(input_dataframe) %in%
                                              c("InputFile","Condition","Treatment","TimePoint","Replicate","Cluster","ThreshCat")]
} else {
  SelectMarker <- IndividualMarkers
}

## Update ClustersToAnalyze as necessary
if (AnalyzeAllClusters & is.null(ClustersToAnalyze)) {
  ClustersToAnalyze <- sort(unique(input_dataframe[["Cluster"]]))
}



# --- QUICK DIAGNOSTICS ---
cat("\n==== DIAGNOSTICS ====\n")
cat("Rows:", nrow(input_dataframe), "\n")
cat("Treatments found:", paste(unique(input_dataframe$Treatment), collapse=", "), "\n")
cat("Clusters found:", paste(sort(unique(input_dataframe$Cluster)), collapse=", "), "\n")
cat("Markers requested:", paste(SelectMarker, collapse=", "), "\n")
cat("Markers present?:", all(SelectMarker %in% colnames(input_dataframe)), "\n")
cat("Marker classes:\n")
print(sapply(SelectMarker, function(m) class(input_dataframe[[m]])))
cat("Marker NA counts:\n")
print(sapply(SelectMarker, function(m) sum(is.na(input_dataframe[[m]]))))
cat("=====================\n\n")


#### Output directory ##########################################################
output_dir <- paste0(
  format(Sys.time(), "%Y-%m-%d"),
  " Linegraphs", input_file_transformation, "_",
  ifelse (test = IsThereThresholdCategory & PlotBy == "ThreshCat", # plot by ThreshCat, or other
          yes = paste0("By", PlotBy,"-", WhichToPlot, "_"), 
          no = paste0("By", PlotBy, "_")),
  paste0(TreatmentName, "_"),
  ifelse(ClusterName == "AnalyzeAllClusters" | ClusterName == "PseudoBulk",
         ClusterName,
         paste0("Cl", ClusterName)),
  ifelse(test = length(WhichToPlot) == 1 & isFALSE(AnalyzeAllClusters) & length(ClustersToAnalyze) > 0,
         yes = paste0("-ThreshCat", paste0(WhichToPlot, collapse = ",")),
         no = ""),
  ifelse(test = RunStatistics,
         yes = paste0("_",Statistical_Test, "-against-", StatisticalComparison_Type, "-",StatisticalComparison_Value),
         no = "")
  )

output_path <- file.path(Original_WD, output_dir)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

cat("Current working directory:", getwd(), "\n")
cat("Original_WD:", Original_WD, "\n")
cat("output_path:", output_path, "\n")
cat("input_filenames:\n")
print(input_filenames)



cat("Does output folder exist? ", dir.exists(output_path), "\n")

## Update variables as necessary for downstream analysis
# Axis ranges for non-marker-compare mode
Max_Values_ByMarker <- sapply(SelectMarker, function(column) max(input_dataframe[[column]], na.rm = TRUE))
Min_Values_ByMarker <- sapply(SelectMarker, function(column) ceiling(sd(input_dataframe[[column]], na.rm = TRUE)) * -1)


# Change WhichToPlot as necessary
#if (IsThereThresholdCategory & length(WhichToPlot) == 0) {
#  WhichToPlot <- sort(unique(DataToPlot[["ThreshCat"]]))
#}


#### Generate functions for analysis ###########################################
### Stats Function for PlotBy_marker == TRUE ==================================
Run_Statistics <- function(data, marker_column) {
  
  if (!RunStatistics) return(NULL)
  
  df <- data
  
  if (StatisticalComparison_Type == "TimePoint") {
    
    df$CompareGroup <- ifelse(df$TimePoint == StatisticalComparison_Value,
                              "Reference", "Other")
    
    results <- df %>%
      group_by(Treatment, Cluster) %>%
      filter(n_distinct(CompareGroup) > 1) %>%
      summarise(
        p_value = tryCatch(
          {
            if (Statistical_Test == "t.test") {
              t.test(.data[[marker_column]] ~ CompareGroup, paired = FALSE)$p.value
            } else {
              wilcox.test(.data[[marker_column]] ~ CompareGroup, paired = FALSE)$p.value
            }
          },
          error = function(e) NA_real_
        ),
        .groups = "drop"
      )
    
  } else if (StatisticalComparison_Type == "Condition") {
    
    df$CompareGroup <- ifelse(df$Treatment == StatisticalComparison_Value,
                              "Reference", "Other")
    
    results <- df %>%
      group_by(Cluster, TimePoint) %>%
      filter(n_distinct(CompareGroup) > 1) %>%
      summarise(
        p_value = tryCatch(
          {
            if (Statistical_Test == "t.test") {
              t.test(.data[[marker_column]] ~ CompareGroup, paired = FALSE)$p.value
            } else {
              wilcox.test(.data[[marker_column]] ~ CompareGroup, paired = FALSE)$p.value
            }
          },
          error = function(e) NA_real_
        ),
        .groups = "drop"
      )
    
  } else if (StatisticalComparison_Type == "Cluster") {
    
    df$CompareGroup <- ifelse(df$Cluster == StatisticalComparison_Value,
                              "Reference", "Other")
    
    results <- df %>%
      group_by(Treatment, TimePoint) %>%
      filter(n_distinct(CompareGroup) > 1) %>%
      summarise(
        p_value = tryCatch(
          {
            if (Statistical_Test == "t.test") {
              t.test(.data[[marker_column]] ~ CompareGroup, paired = FALSE)$p.value
            } else {
              wilcox.test(.data[[marker_column]] ~ CompareGroup, paired = FALSE)$p.value
            }
          },
          error = function(e) NA_real_
        ),
        .groups = "drop"
      )
  }
  
  results <- results %>%
    mutate(Significance =
             case_when(
               p_value < 0.001 ~ "***",
               p_value < 0.01  ~ "**",
               p_value < 0.05  ~ "*",
               TRUE ~ ""
             ))
  
  return(results)
}


### Stats Function ============================================================
MakeStats_ConditionByTime <- function(df, marker, ref_treatment = "NoTreat", ref_timepoint = 0, test = "auto") {
  
  ## Testing
  # df <- DataToPlot
  # marker <- b
  # ref_treatment <- "NoTreat"
  # ref_timepoint <- 0
  # test <- Statistical_Test
  
  # 1. Prep IDs
  df_mod <- df %>% mutate(Group_ID = factor(paste(Treatment, TimePoint, sep = "_")))
  ref_id <- paste(ref_treatment, ref_timepoint, sep = "_")
  
  # 2. QC Diagnostics
  # Shapiro-Wilk for Normality
  shapiro_res <- shapiro.test(residuals(lm(as.formula(paste0("`", marker, "` ~ Group_ID")), data = df_mod)))
  is_normal <- shapiro_res$p.value > 0.05
  
  # Levene's for Homogeneity of Variance
  levene_res <- car::leveneTest(as.formula(paste0("`", marker, "` ~ Group_ID")), data = df_mod)
  is_homogenous <- levene_res$`Pr(>F)`[1] > 0.05
  
  # 3. Automated Test Selection Logic
  if (test == "auto") {
    # First level: Parametric vs Non-Parametric
    test_type <- if (is_normal && is_homogenous) "ANOVA" else "DUNN"
    
    # Second level: If ANOVA, decide One-Way vs Two-Way
    if (test_type == "ANOVA") {
      # Check if we have multiple factors and they are appropriate (categorical or coerced factors)
      # We use 'is.numeric' check to decide if we need to force factor conversion later
      if (!is.null(df_mod$Treatment) && !is.null(df_mod$TimePoint)) {
        test <- "ANOVA_two"
      } else {
        test <- "ANOVA_one"
      }
    } else {
      test <- "DUNN"
    }
  }
  
  # 4. Warnings for mismatched tests
  if (test == "ANOVA" && (!is_normal || !is_homogenous)) {
    warning("CAUTION: You are using ANOVA, but data failed Normality or Variance checks. Results may be unreliable.")
  }
  
  # 5. Statistical Execution
  # Initialize results as NULL to prevent "object not found" errors
  results <- NULL
  
  if (grepl("ANOVA", test)) {
    
    if (test == "ANOVA_one" || test == "ANOVA") {
      # --- One-way ANOVA Logic ---
      df_mod <- mutate(df_mod, Group_ID = relevel(Group_ID, ref = ref_id))
      model <- aov(as.formula(paste0("`", marker, "` ~ Group_ID")), data = df_mod)
      post_hoc <- summary(glht(model, linfct = mcp(Group_ID = "Dunnett")))
      
      results <- data.frame(
        Group_ID = names(post_hoc$test$tstat),
        p_value = as.numeric(post_hoc$test$pvalues)
      ) %>%
        mutate(Group_ID = gsub(paste0(" - ", ref_id), "", Group_ID))
      
    } else if (test == "ANOVA_two") {
      # 1. Run the model
      df_mod <- df_mod %>%
        mutate(Treatment = factor(Treatment), TimePoint = factor(TimePoint))
      
      model <- aov(as.formula(paste0("`", marker, "` ~ Treatment * TimePoint")), data = df_mod)
      post_hoc <- TukeyHSD(model)
      
      # 2. Extract the interaction results (Treatment:TimePoint)
      tukey_df <- as.data.frame(post_hoc[[3]]) 
      
      # 3. CRITICAL: Reformat to match One-Way ANOVA expectations
      # We filter for only comparisons against the reference ID to keep the scale small
      results <- data.frame(
        Full_Comparison = rownames(tukey_df),
        p_value = tukey_df$`p adj`
      ) %>%
        # We only want rows where our reference group is being compared
        filter(grepl(ref_id, Full_Comparison)) %>%
        # Clean the name so it's just the 'Other' group name, matching your One-Way format
        mutate(Group_ID = gsub(paste0(":", ref_id, "|", ref_id, ":|-"), "", Full_Comparison),
               Group_ID = trimws(Group_ID)) %>%
        select(Group_ID, p_value)
      
    }
  
  # 6. Final Formatting with QC info
  results <- results %>%
    separate(Group_ID, into = c("Treatment", "TimePoint"), sep = "_", convert = TRUE) %>%
    mutate(
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ "ns"
      ),
      QC_Normality_P = shapiro_res$p.value,
      QC_Variance_P = levene_res$`Pr(>F)`[1],
      Test_Used = test
    )
  return (results)
}


# ## Validation of MakeStats_ConditionByTime_QC
# # 1. Generate Mock Data: This code creates a scenario where the 16hr time point has much higher variance, which often triggers a "failed" QC check for ANOVA.
# set.seed(42)
# mock_df <- data.frame(
#   Treatment = c(rep("NoTreat", 4), rep("BDNF", 12)),
#   TimePoint = c(rep(0, 4), rep(c(0.25, 1, 4, 16), each = 3)),
#   pS6_175Lu = c(
#     rnorm(4, mean = 1, sd = 0.2),   # NoTreat 0hr
#     rnorm(3, mean = 7, sd = 0.5),   # BDNF 15min
#     rnorm(3, mean = 6, sd = 0.5),   # BDNF 1hr
#     rnorm(3, mean = 8, sd = 0.8),   # BDNF 4hr
#     rnorm(3, mean = 3, sd = 2.0)    # BDNF 16hr (High variance!)
#   )
# )
# 
# # 2. Run the Function: set test = "auto" so it evaluates the data automatically
# stats_results <- MakeStats_ConditionByTime(mock_df, marker = "pS6_175Lu")
# # View the results
# print(stats_results)
# 
# # 3. Summary of QC Columns
# # When you inspect stats_results, look for these indicators to ensure your math is "clean":
# # Column            # What it tells you                         # Goal
# # QC_Normality_P    Is the data normally distributed?           > 0.05 (Fail means the data is skewed)
# # QC_Variance_P     Is the spread of data equal across groups?  > 0.05 (Fail means some groups are more "noisy" than others)
# # Test_Used         Which math was applied?                     Should be ANOVA if QC passes, DUNN if it fails.


### Plotting Function =========================================================
Plotting_Function <- function(InputDataframe, Treatments, Clusters, Marker) {

  ## Testing
  # InputDataframe <- input_dataframe
  # Treatments <- TreatmentToAnalyze
  # Clusters <- ClustersToAnalyze #4 # ClustersToAnalyze # From HowManyLoops
  # Marker <- SelectMarker
  # b <- "H3K9ac_115In" # "pS6_175Lu"
  
  
  df <- InputDataframe
  
  # Filter by Treatments
  if (length(Treatments) > 0) {
    DataToAnalyze_AllClusters <- df[df$Treatment %in% Treatments, ]
  } else {
    DataToAnalyze_AllClusters <- df
  }
  
  if (nrow(DataToAnalyze_AllClusters) == 0) {
    stop("After filtering by TreatmentToAnalyze, there are 0 rows. Check Treatment names.")
  }
  
  # Filter by Clusters
  if (AnalyzeAllClusters | length(Clusters) > 0) {
    DataToAnalyze <- DataToAnalyze_AllClusters[DataToAnalyze_AllClusters$Cluster %in% Clusters, ] %>%
      mutate(Cluster = factor(as.character(Cluster), levels = Clusters))
  } else {
    DataToAnalyze <- DataToAnalyze_AllClusters %>%
      mutate(Cluster = as.character(Cluster))
  }
  
  # TimePoint factor if categorical
  if (TimePoint_VariableType == "Categorical") {
    DataToAnalyze <- DataToAnalyze %>%
      mutate(TimePoint = factor(as.character(TimePoint),
                                 levels = as.character(sort(unique(as.numeric(TimePoint))))))
  }
  
  # Treatment factor
  DataToAnalyze <- DataToAnalyze %>%
    mutate(Treatment = factor(Treatment, levels = TreatmentToAnalyze))
  
  Plot_List <- list()
  # Table_List <- list()
  
  if (PlotBy_marker) {
    
    # Safe 0-1 normalization
    range01 <- function(x){
      r <- range(x, na.rm = TRUE)
      if (!is.finite(r[1]) || !is.finite(r[2]) || diff(r) == 0) return(rep(0, length(x)))
      (x - r[1]) / diff(r)
    }
    
    non_meta_cols <- colnames(DataToAnalyze)[!colnames(DataToAnalyze) %in%
                       c("InputFile","Condition","Treatment","Cluster","TimePoint","Replicate","ThreshCat")]
    
    if (IsThereThresholdCategory & length(WhichToPlot) > 0) {
      DataToPlot <- DataToPlot %>%
        dplyr::filter(ThreshCat %in% WhichToPlot)
    }
    
    
    
    DataToPlot_01 <- DataToAnalyze %>%
      mutate(across(all_of(non_meta_cols), ~ range01(.x)))
    
    # FIX: loop variable renamed from 'c' (shadows base R) to 'cluster_id'
    for (cluster_id in unique(DataToAnalyze$Cluster)) {
      
      DataToPlot <- DataToPlot_01[, c("Condition","Treatment","TimePoint","Cluster","Replicate", all_of(Marker), "ThreshCat")] %>%
        dplyr::filter(Cluster == cluster_id)
      
      DataToPlot_pivoted <- DataToPlot %>%
        pivot_longer(cols = all_of(Marker))
      
      DataToPlot_pivoted_summary <- DataToPlot_pivoted %>%
        group_by(Condition, Treatment, Cluster, TimePoint, ThreshCat, name) %>%
        summarize(sd = sd(value), value = mean(value), .groups = "drop")
      
      # FIX: compare TimePoint as character consistently
      treatments_with_0hr <- DataToPlot_pivoted_summary %>%
        dplyr::filter(as.character(TimePoint) == "0") %>%
        pull(Treatment) %>% unique()
      
      treatments_without_0hr <- setdiff(unique(DataToPlot_pivoted_summary$Treatment), treatments_with_0hr)
      
      append_0hr <- lapply(treatments_without_0hr, function(trt) {
        DataToPlot_pivoted_summary %>%
          dplyr::filter(as.character(TimePoint) == "0") %>%
          mutate(Condition = trt, Treatment = trt)
      }) %>% bind_rows()
      
      DataToPlot_pivoted_summary <- bind_rows(DataToPlot_pivoted_summary, append_0hr) %>%
        mutate(Treatment = factor(Treatment, levels = TreatmentToAnalyze))
      
      # Define colors to use
      if (length(TreatmentToAnalyze) > 3) { 
        colors <- ColorToUse
        } else if (PlotBy == "Treatment" | 
          PlotBy == "ThreshCat") {
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
      
      # FIX: Run statistics using correct loop variable (cluster_id, not b)
      stats_df <- NULL
      if (RunStatistics) {
        stats_df <- Run_Statistics(DataToPlot_pivoted, "value")
      }
      cat("\n--- STATS CHECK: cluster", cluster_id, "---\n")
      print(head(stats_df, 10))
      cat("nrow(stats_df) =", ifelse(is.null(stats_df), "NULL", nrow(stats_df)), "\n")
      cat("non-empty sig count =", ifelse(is.null(stats_df), "NULL", sum(stats_df$Significance != "")), "\n")
      
      # FIX: Build base plot first, then add stats layer separately (no 'if' inside '+' chain)
      All_TimePoint <- ggplot(DataToPlot_pivoted, aes(x = TimePoint, y = value, color = name)) +
        geom_point(size = 1, alpha = 0.5) +
        geom_ribbon(
          data = DataToPlot_pivoted_summary,
          aes(x = TimePoint, y = value, ymin = value - sd, ymax = value + sd, fill = name, group = name),
          alpha = 0.2, color = NA, inherit.aes = FALSE
        ) +
        geom_point(
          data = DataToPlot_pivoted_summary[DataToPlot_pivoted_summary$Treatment %in% treatments_without_0hr, ],
          aes(y = value), size = 2
        ) +
        geom_line(
          data = DataToPlot_pivoted_summary[DataToPlot_pivoted_summary$Treatment %in% treatments_without_0hr, ],
          aes(y = value, group = name), linewidth = 0.5
        ) +
        geom_point(
          data = DataToPlot_pivoted_summary[DataToPlot_pivoted_summary$Treatment %in% treatments_with_0hr, ],
          aes(y = value), size = 2
        ) +
        labs(
          color = paste0("Signaling Marker\nCluster ", cluster_id),
          fill  = paste0("Signaling Marker\nCluster ", cluster_id),
          x = "Time (hr)",
          y = if (input_file_transformation == "_FC") "Fold Change (Marker / 0hr) - 0-1 Normalized" else "log2FC(Marker / 0hr) - 0-1 Normalized"
        ) +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        # FIX: use coord_cartesian so stats annotations above 1 are not clipped
        coord_cartesian(ylim = c(0, 1.15)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      # FIX: Add stats annotations OUTSIDE the ggplot '+' chain
      if (!is.null(stats_df) && nrow(stats_df) > 0 && any(stats_df$Significance != "")) {
        stats_df$TimePoint <- factor(
          as.character(stats_df$TimePoint),
          levels = levels(DataToPlot$TimePoint)
        )
        stats_df$y <- 1.08
        
        All_TimePoint <- All_TimePoint +
          geom_text(
            data = stats_df %>% dplyr::filter(Significance != ""),
            aes(x = TimePoint, y = y, label = Significance, color = PlotBy),
            inherit.aes = FALSE,
            size = 4
          )
      }
      
      Output_Plot <- ggdraw() + draw_plot(All_TimePoint)
      Plot_List[[as.character(cluster_id)]] <- Output_Plot
    }
    
  } else {
    
    for (b in Marker) {
      
      DataToPlot <- DataToAnalyze[, c("Condition","Treatment","TimePoint","Cluster","Replicate", b, "ThreshCat")] %>%
        dplyr::filter(Cluster %in% if (length(ClustersToAnalyze) > 1) ClustersToAnalyze else unique(DataToAnalyze$Cluster))
      
      
      if (IsThereThresholdCategory & length(WhichToPlot) > 0) {
        DataToPlot <- DataToPlot %>% 
          dplyr::filter(ThreshCat %in% WhichToPlot)
      }
      
      DataToPlot_summary <- DataToPlot %>%
        group_by(Condition, Treatment, Cluster, TimePoint, ThreshCat) %>%
        summarize(Mean = mean(.data[[b]]), SD = sd(.data[[b]]), .groups = "drop")
      
      treatments_with_0hr <- DataToPlot_summary %>%
        dplyr::filter(as.character(TimePoint) == "0") %>%
        pull(Treatment) %>% unique()
      treatments_without_0hr <- setdiff(unique(DataToPlot_summary$Treatment), treatments_with_0hr)
      
      append_0hr <- lapply(treatments_without_0hr, function(trt) {
        DataToPlot_summary %>%
          dplyr::filter(as.character(TimePoint) == "0") %>%
          mutate(Condition = trt, Treatment = trt)
      }) %>% bind_rows()
      
      DataToPlot_summary <- bind_rows(DataToPlot_summary, append_0hr) %>%
        mutate(Treatment = factor(Treatment, levels = TreatmentToAnalyze))
      
      # Define colors to use
      if (length(TreatmentToAnalyze) > 3) { 
        colors <- ColorToUse
      } else if (PlotBy == "Treatment" | 
          PlotBy == "ThreshCat") {
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
      
      y_upper <- Max_Values_ByMarker[[b]] + abs(Min_Values_ByMarker[[b]])
      y_lower <- Min_Values_ByMarker[[b]]
      
      # FIX: Build base plot first using coord_cartesian so annotation headroom is available
      All_TimePoint <- ggplot(DataToPlot, aes(x = TimePoint, y = .data[[b]], color = !!sym(PlotBy))) +
        geom_point(size = 1, alpha = 0.5) +
        geom_ribbon(
          data = DataToPlot_summary,
          aes(x = TimePoint, y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = !!sym(PlotBy), group = !!sym(PlotBy)),
          alpha = 0.2, color = NA, inherit.aes = FALSE
        ) +
        geom_point(
          data = DataToPlot_summary[DataToPlot_summary$Treatment %in% treatments_without_0hr, ],
          aes(y = Mean), size = 2
        ) +
        geom_line(
          data = DataToPlot_summary[DataToPlot_summary$Treatment %in% treatments_without_0hr, ],
          aes(y = Mean, group = !!sym(PlotBy)), linewidth = 0.5
        ) +
        geom_point(
          data = DataToPlot_summary[DataToPlot_summary$Treatment %in% treatments_with_0hr, ],
          aes(y = Mean), size = 2
        ) +
        labs(
          title = paste0(b, ": ",
                         ifelse(test = length(Clusters) > 0,
                         yes = paste0("Cluster ", paste0(Clusters, collapse = ",")), 
                         no = "Pseudobulk")
                         ),
          color = PlotBy,
          fill  = PlotBy,
          x = "Time (hr)",
          y = if (input_file_transformation == "_FC") paste0("Fold Change (", b, " / 0hr)") else paste0("log2FC(", b, " / 0hr)")
        ) +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        # FIX: coord_cartesian instead of ylim so geom_text stars above panel are not clipped
        coord_cartesian(ylim = c(y_lower, y_upper * 1.12)) #+
        # theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      # FIX: Run statistics ONCE using MakeStats_ConditionByTime, then add to plot
      # (previously stats_df was computed twice - first result was silently overwritten)
      stats_df <- NULL
      if (RunStatistics) {
        stats_df <- MakeStats_ConditionByTime(
          DataToPlot,
          b,
          ref_treatment = "NoTreat",
          ref_timepoint = 0,
          test = Statistical_Test
        
        )
      }
      
      # FIX: Add stats annotations OUTSIDE the ggplot '+' chain
      if (!is.null(stats_df) && nrow(stats_df) > 0 && any(stats_df$Significance != "")) {
        stats_df$TimePoint <- factor(
          as.character(stats_df$TimePoint),
          levels = levels(DataToPlot$TimePoint)
        )
        
        for (i in 1:length(unique(stats_df[[PlotBy]]))) {
          stats_df[which(stats_df[[PlotBy]] == unique(stats_df[[PlotBy]])[i]), "y"] <- (y_upper * 1.06) - (stats_yAxis_scaler*(i-1))
          # (y_upper * 1.06) - (i-1) / (y_upper * 1.06) - (2*(i-1)) / (y_upper * 0.8) - (5*(i-1)) 
        }
        
        # Add statistics to plot and add subtitle to include stats test
        All_TimePoint <- All_TimePoint +
          geom_text(
            data = stats_df %>% dplyr::filter(Significance != ""),
            aes(x = TimePoint, y = y, label = Significance, color = .data[[PlotBy]]),
            inherit.aes = FALSE,
            size = 10
            ) +
          labs(subtitle = if (unique(stats_df[["Test_Used"]]) == "ANOVA") {
            "One-way ANOVA with Dunnett's Post-hoc Test"
            } else if (unique(stats_df[["Test_Used"]]) == "DUNN") {
              "Dunn's Test with Bonferroni Post-hoc Test"
              } else {
                unique(stats_df[["Test_Used"]])
                }
            )
      }
      
      if (List_p_values && !is.null(stats_df)) {
        TablePlot <- grid.arrange(top = paste0(b, ": ",
                                               ifelse(test = length(Clusters) > 0,
                                                      yes = paste0("Cluster ", paste0(Clusters, collapse = ",")), 
                                                      no = "Pseudobulk"),
                                               "\n", if (unique(stats_df[["Test_Used"]]) == "ANOVA") {
                                                 "One-way ANOVA with Dunnett's Post-hoc Test"
                                               } else if (unique(stats_df[["Test_Used"]]) == "DUNN") {
                                                 "Dunn's Test with Bonferroni Post-hoc Test"
                                               } else {
                                                 unique(stats_df[["Test_Used"]])
                                               }
                                               ), 
                                  tableGrob(d = stats_df[, c("TimePoint", PlotBy, "p_value", "Significance")],
                                            theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.1, fontsize = 10)),
                                                                   rowhead = list(fg_params = list(fontface = "bold", fontsize = 10)))))
        
        Output_Table <- ggdraw() + draw_plot(TablePlot)
        
        Plot_List[[paste0("Table: ", b)]] <- Output_Table
      }
      
     
      Output_Plot <- ggdraw() + draw_plot(All_TimePoint)
      Plot_List[[b]] <- Output_Plot
    }
  }
  
  message("Done building plot list")
  return(Plot_List)
  # return()
}

#### Generate and save plots ###################################################
if (PlotAsGrid) {
  
  if (AnalyzeAllClusters & !PlotBy_marker) {
    HowManyLoops <- sort(unique(input_dataframe$Cluster))
  } else {
    HowManyLoops <- 1
  }
  
  lapply(HowManyLoops, function(clusterID) {
    
    if (length(HowManyLoops) == 1) {
      clusterID <- ClustersToAnalyze
    }
    
    Output_List <- Plotting_Function(
      InputDataframe = input_dataframe,
      Treatments = TreatmentToAnalyze,
      Clusters = clusterID,
      Marker = SelectMarker
    )
    
    
    Export_List <- list(Tables = Output_List[grepl("Table", names(Output_List))],
                        Plots = Output_List[!grepl("Table", names(Output_List))])
    
    fixed_ncol <- plot_column_length
    fixed_nrow <- round(length(Export_List[["Plots"]]) / fixed_ncol) # 4
    fixed_slots <- fixed_ncol * fixed_nrow
    
    
    # Add void plots to empty slots in the grid
    for (x in names(Export_List)) {
      if (length(Export_List[[x]]) < fixed_slots) {
        voidPlot <- ggplot() + theme_void()
        
        Export_List[[x]] <- append(Export_List[[x]], 
                                   replicate(n = fixed_slots - length(Export_List[[x]]), 
                                             expr = voidPlot))
        
      } else if (length(Export_List[[x]]) > fixed_slots) {
        Export_List[[x]] <- Export_List[[x]][1:fixed_slots]
      }
      }
    
    WhichToExport <- c()
    for (f in seq_along(Export_List)) {
      WhichToExport <- c(
        WhichToExport,
        if (is.null(Export_List[[f]]) || length(Export_List[[f]]) == 0) "" else names(Export_List)[f]
      )
    }
    WhichToExport <- WhichToExport[WhichToExport != ""]
    
    
    sapply(WhichToExport, function (export) {
      
      PlotsInGrid <- grid.arrange(
        grobs = Export_List[[export]],
        nrow = fixed_nrow,
        ncol = fixed_ncol
      )
      
      grid_plot <- grid.arrange(
        textGrob(
          ifelse(ClusterName != "PseudoBulk" & ClusterName != "AnalyzeAllClusters",
                 paste0("Cluster", paste0(clusterID, collapse = ","), " Plots"),
                 paste0(ClusterName, " Plots")),
          gp = gpar(fontsize = 18, fontface = "bold")
        ),
        PlotsInGrid,
        ncol = 1,
        heights = c(0.05, 0.95)
      )
      
      plot_width  <- 5 * fixed_ncol
      plot_height <- 5 * fixed_nrow
      
      ggsave(
        filename = file.path(
          output_path,
          paste0(
            "Linegraph_",
            ifelse(test = ClusterName == "PseudoBulk", 
                   yes = paste0(ClusterName, "_"), 
                   no = "Clusters_"),
            paste0("By", PlotBy),
            "-Grid",
            if (IsThereThresholdCategory) {
              paste0("-CompareCat", paste0(WhichToPlot, collapse = ","))
            } else {
              ""
            },
            # ifelse(test = ClusterName != "PseudoBulk",
            #        yes = paste0("_Cluster", paste0(clusterID, collapse = "-")),
            #        no = ""),
            "-", export,
            ".", Output_File_Type
          )
        ),
        plot = grid_plot,
        width = plot_width,
        height = plot_height,
        units = "in",
        dpi = 300
      )
      
    })
    
    message("Saved grid to: ", output_path)
  })
  
} else {
  stop("PlotAsGrid == FALSE branch not included in this paste. Set PlotAsGrid <- TRUE.")
}

  
  message("SUCCESS! File saved at: ", pdf_file)
 
#### Reset directory ###########################################################
setwd(Original_WD)
message("Done. Output folder: ", output_path)

  

