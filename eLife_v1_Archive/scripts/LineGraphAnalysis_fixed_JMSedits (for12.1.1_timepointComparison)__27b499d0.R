# Title: LineGraphAnalysis
# Version: 11 (all bugs fixed)
# Author: Jonathon Sewell (edited for robustness)
# Date: 2025-11-10

# Update

library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(cowplot)
library(viridis)
library(gridExtra)
library(grid)

Original_WD <- getwd()

### User-defined Variables ====================================================
### Statistical comparison settings ==========================================
RunStatistics <- TRUE   # TRUE/FALSE
List_p_values <- TRUE   # TRUE/FALSE

StatisticalComparison_Type  <- "TimePoint"   # "Condition", "TimePoint", "Cluster"
StatisticalComparison_Value <- 0             # e.g. 0 for 0hr
Statistical_Test <- "t.test"                 # "t.test" or "wilcox"
Alpha <- 0.05

## Marker-related variables
PlotBy_marker <- FALSE # TRUE/FALSE

MarkerSelect <- "Individual" # "Individual" / "All"
if (PlotBy_marker) MarkerSelect <- "Individual"

IndividualMarkers <- c(
  "pERK_171Yb",
  "pAKt_168Er",
  "pCreb_176Lu",
  "pGSK3b_139La",
  "pS6_175Lu"
  )

IsThereThresholdCategory <- FALSE # TRUE/FALSE

## Cluster-related variables
AnalyzeAllClusters <- FALSE # TRUE / FALSE (if selecting specific clusters or if using pseudobulk data)
ClustersToAnalyze <- c(1, 2, 4, 7, 13) # c(1, 2, 4, 7, 13) / c(9, 17, 18) / c(1, 3, 15, 17) / c()

if (AnalyzeAllClusters) {
  ClusterName <- "AnalyzeAllClusters"
} else if (!AnalyzeAllClusters & length(ClustersToAnalyze) == 0) {
  ClusterName <- "PseudoBulk"
} else {
  ClusterName <- paste(ClustersToAnalyze, collapse = ",")
}

## Treatments
TreatmentToAnalyze <- c("NoTreat", "BDNF") # "NoTreat" / "ContDepr" / "CompMed" / "BDNF"
TreatmentName <- paste(TreatmentToAnalyze, collapse = ",")

## Define PlotBy
if (PlotBy_marker) {
  if (length(IndividualMarkers) > 1) {
    PlotBy <- "Marker"
  } else {
    stop("Select more IndividualMarkers to compare")
  }
} else if (IsThereThresholdCategory) {
  PlotBy <- "ThresholdCategory"
} else if (length(ClustersToAnalyze) == 0) {
  PlotBy <- "Treatment"
} else {
  PlotBy <- "Cluster"
}

## Time handling
TimePoint_VariableType <- "Categorical" # "Categorical" / "Continuous"
EarlyTimePoint_Cutoff <- 5

## Output info
PlotAsGrid <- TRUE
plot_column_length <- 5

if (PlotBy == "Cluster") {
  ColorToUse <- "Rainbow"
} else {
  ColorToUse <- "A" # With viridis, A = magma / B = inferno / C = plasma / D = viridis / E = cividis / F = rocket / G = mako / H = turbo
}

Plot_Width <- 14
Plot_Height <- 10
Output_File_Type <- "pdf" # "pdf", "png", etc.

### Read in and prep data =====================================================
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

input_dataframe_all <- NULL
input_dataframe_all_avg <- NULL

for (i in seq_along(input_filenames)) {
  f <- input_filenames[i]
  
  if (!grepl("_RepsAvg", f)) {
    tmp <- read_input_file(file.path(Original_WD, f)) %>%
      mutate(
        `...1` = sub("_.*", "", f),
        ThresholdCategory = sub(".*FC-([^\\.]+)\\.csv", "\\1", f)
      ) %>%
      rename(InputFile = `...1`)
    
    input_dataframe_all <- bind_rows(input_dataframe_all, tmp)
    
  } else {
    tmp <- read_input_file(file.path(Original_WD, f)) %>%
      mutate(
        `...1` = sub("_.*", "", f),
        ThresholdCategory = sub(".*RepsAvg-([^\\.]+)\\.csv", "\\1", f)
      ) %>%
      rename(InputFile = `...1`)
    
    input_dataframe_all_avg <- bind_rows(input_dataframe_all_avg, tmp)
  }
}

if (is.null(input_dataframe_all) || nrow(input_dataframe_all) == 0) {
  stop("Read-in failed: input_dataframe_all is empty. Check your file contents/format.")
}

if (IsThereThresholdCategory) {
  factor_order <- c("hi",
                    unique(input_dataframe_all$ThresholdCategory)[-which(unique(input_dataframe_all$ThresholdCategory) %in% c("hi", "lo"))],
                    "lo")
  input_dataframe_all$ThresholdCategory <- factor(input_dataframe_all$ThresholdCategory, levels = factor_order)
  if (!is.null(input_dataframe_all_avg) && nrow(input_dataframe_all_avg) > 0) {
    input_dataframe_all_avg$ThresholdCategory <- factor(input_dataframe_all_avg$ThresholdCategory, levels = factor_order)
  }
}

RemoveOverlappingConditions <- function(DATAFRAME) {
  ColumnOfInterest <- "InputFile"
  reference_group <- unique(DATAFRAME[[ColumnOfInterest]])[1]
  reference_conditions <- dplyr::filter(DATAFRAME, .data[[ColumnOfInterest]] == reference_group) %>%
    pull(Condition)
  
  other_groups <- setdiff(unique(DATAFRAME[[ColumnOfInterest]]), reference_group)
  
  dplyr::filter(
    DATAFRAME,
    !( .data[[ColumnOfInterest]] %in% other_groups &
         .data[["Condition"]] %in% reference_conditions )
  )
}

# Create final dataframes explicitly (no assign())
input_dataframe <- RemoveOverlappingConditions(input_dataframe_all) %>%
  mutate(Treatment = sub("_[^_]+$", "", .data$Condition), .after = "Condition")

meta_cols <- c("InputFile","Condition","Treatment","TimePoint","Replicate","Cluster","ThresholdCategory")
marker_cols <- setdiff(colnames(input_dataframe), meta_cols)

input_dataframe <- input_dataframe %>%
  mutate(across(all_of(marker_cols), ~ readr::parse_number(as.character(.x))))

input_dataframe_avg <- NULL
if (!is.null(input_dataframe_all_avg) && nrow(input_dataframe_all_avg) > 0) {
  input_dataframe_avg <- RemoveOverlappingConditions(input_dataframe_all_avg) %>%
    mutate(Treatment = sub("_[^_]+$", "", .data$Condition), .after = "Condition") %>%
    mutate(across(all_of(setdiff(colnames(.), meta_cols)), ~ readr::parse_number(as.character(.x))))
}

# Marker selection
if (MarkerSelect == "All") {
  SelectMarker <- colnames(input_dataframe)[!colnames(input_dataframe) %in%
                                              c("InputFile","Condition","Treatment","TimePoint","Replicate","Cluster","ThresholdCategory")]
} else {
  SelectMarker <- IndividualMarkers
}

# Update ClustersToAnalyze as necessary
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

# Output directory
output_dir <- paste0(
  format(Sys.time(), "%Y-%m-%d"),
  " Linegraphs_",
  ifelse(IsThereThresholdCategory & PlotBy == "ThresholdCategory",
         paste0("By", PlotBy, "-CompareCat_"),
         paste0("By", PlotBy, "_")),
  paste0(TreatmentName, "_"),
  ifelse(ClusterName == "AnalyzeAllClusters" | ClusterName == "PseudoBulk",
         ClusterName,
         paste0("Cl", ClusterName)),
  ifelse(test = RunStatistics,
         yes = paste0("_",Statistical_Test, "-against-", StatisticalComparison_Type, "-",StatisticalComparison_Value),
         no = "")
  )

output_path <- file.path(Original_WD, output_dir)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Axis ranges for non-marker-compare mode
Max_Values_ByMarker <- sapply(SelectMarker, function(column) max(input_dataframe[[column]], na.rm = TRUE))
Min_Values_ByMarker <- sapply(SelectMarker, function(column) ceiling(sd(input_dataframe[[column]], na.rm = TRUE)) * -1)


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

Add_Stats_Ypos <- function(stats_df, y_max) {
  if (is.null(stats_df) || nrow(stats_df) == 0) return(NULL)
  stats_df %>%
    mutate(
      TimePoint = factor(as.character(TimePoint), levels = levels(stats_df$TimePoint)),
      y = y_max * 1.05
    )
}

# FIX: MakeStats_ConditionByTime now uses .data[[marker]] correctly inside summarise


MakeStats_ConditionByTime <- function(df, marker, ref_treatment = "NoTreat", ref_timepoint = 0, test = "wilcox") {
  
  # df <- DataToPlot
  # marker <- b
  # ref_treatment <- "NoTreat"
  # ref_timepoint <- 0
  # test <- Statistical_Test
    
  
  
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
    group_by(TimePoint, Cluster) %>%
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



### Plotting Function =========================================================
Plotting_Function <- function(InputDataframe, Treatments, Clusters, Marker) {

  
  ## Test
  # InputDataframe <- input_dataframe
  # Treatments <- TreatmentToAnalyze
  # Clusters <- ClustersToAnalyze #4 # ClustersToAnalyze # From HowManyLoops
  # Marker <- SelectMarker
  # b<- "pS6_175Lu"
  
  
  
  
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
  
  if (PlotBy_marker) {
    
    # Safe 0-1 normalization
    range01 <- function(x){
      r <- range(x, na.rm = TRUE)
      if (!is.finite(r[1]) || !is.finite(r[2]) || diff(r) == 0) return(rep(0, length(x)))
      (x - r[1]) / diff(r)
    }
    
    non_meta_cols <- colnames(DataToAnalyze)[!colnames(DataToAnalyze) %in%
                       c("InputFile","Condition","Treatment","Cluster","TimePoint","Replicate","ThresholdCategory")]
    
    DataToPlot_01 <- DataToAnalyze %>%
      mutate(across(all_of(non_meta_cols), ~ range01(.x)))
    
    # FIX: loop variable renamed from 'c' (shadows base R) to 'cluster_id'
    for (cluster_id in unique(DataToAnalyze$Cluster)) {
      
      DataToPlot <- DataToPlot_01[, c("Condition","Treatment","TimePoint","Cluster","Replicate", all_of(Marker), "ThresholdCategory")] %>%
        dplyr::filter(Cluster == cluster_id)
      
      DataToPlot_pivoted <- DataToPlot %>%
        pivot_longer(cols = all_of(Marker))
      
      DataToPlot_pivoted_summary <- DataToPlot_pivoted %>%
        group_by(Condition, Treatment, Cluster, TimePoint, ThresholdCategory, name) %>%
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
      if (PlotBy == "Treatment" | 
          PlotBy == "ThresholdCategory") {
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
      
      DataToPlot <- DataToAnalyze[, c("Condition","Treatment","TimePoint","Cluster","Replicate", b, "ThresholdCategory")] %>%
        dplyr::filter(Cluster %in% if (length(ClustersToAnalyze) > 1) ClustersToAnalyze else unique(DataToAnalyze$Cluster))
      
      DataToPlot_summary <- DataToPlot %>%
        group_by(Condition, Treatment, Cluster, TimePoint, ThresholdCategory) %>%
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
      if (PlotBy == "Treatment" | 
          PlotBy == "ThresholdCategory") {
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
        coord_cartesian(ylim = c(y_lower, y_upper * 1.12)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
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
        # stats_df$y <- y_upper * 1.06
        for (i in 1:length(unique(stats_df[["TEST"]]))) {
          stats_df[which(stats_df[["Cluster"]] == unique(stats_df[["Cluster"]])[i]), "y"] <- (y_upper * 1.06) - (i-1)
        }
        
        All_TimePoint <- All_TimePoint +
          geom_text(
            data = stats_df %>% dplyr::filter(Significance != ""),
            aes(x = TimePoint, y = y, label = Significance, color = Cluster),
            inherit.aes = FALSE,
            size = 10
          )
      }
      
      if (List_p_values) {
        All_TimePoint <- All_TimePoint + 
          annotation_custom(tableGrob(d = stats_df[, c("TimePoint", "Cluster", "p_value")],
                                      theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.1, fontsize = 10)),
                                                             rowhead = list(fg_params = list(fontface = "bold", fontsize = 10)))), 
                            xmin = nrow(stats_df) * 0.75, xmax = nrow(stats_df) * 1.1, 
                            ymin = 0, ymax = 0)
      }
      
      
      
      Output_Plot <- ggdraw() + draw_plot(All_TimePoint)
      Plot_List[[b]] <- Output_Plot
    }
  }
  
  message("Done building plot list")
  return(Plot_List)
}

### Generate and save plots ===================================================
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
    
    OutputPlots_List <- Plotting_Function(
      InputDataframe = input_dataframe,
      Treatments = TreatmentToAnalyze,
      Clusters = clusterID,
      Marker = SelectMarker
    )
    
    fixed_ncol <- plot_column_length
    fixed_nrow <- round(length(OutputPlots_List) / 5) # 4
    fixed_slots <- fixed_ncol * fixed_nrow
    
    if (length(OutputPlots_List) < fixed_slots) {
      for (k in (length(OutputPlots_List) + 1):fixed_slots) {
        OutputPlots_List[[k]] <- ggplot() + theme_void()
      }
    } else if (length(OutputPlots_List) > fixed_slots) {
      OutputPlots_List <- OutputPlots_List[1:fixed_slots]
    }
    
    PlotsInGrid <- grid.arrange(
      grobs = OutputPlots_List,
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
          if (IsThereThresholdCategory) "-CompareCat" else "",
          "-Grid",
          ifelse(test = ClusterName != "PseudoBulk",
                 yes = paste0("_Cluster", paste0(clusterID, collapse = ",")),
                 no = ""),
          ".", Output_File_Type
        )
      ),
      plot = grid_plot,
      width = plot_width,
      height = plot_height,
      units = "in",
      dpi = 300
    )
    
    message("Saved grid to: ", output_path)
  })
  
} else {
  stop("PlotAsGrid == FALSE branch not included in this paste. Set PlotAsGrid <- TRUE.")
}

setwd(Original_WD)
message("Done. Output folder: ", output_path)
