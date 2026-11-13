################################################################################
# 02.1_SubsetAndExport.R
################################################################################
# Layer:    DataCleanup -- OPTIONAL step, between clustering rounds
#           (run after 02_FCS_to_CSV, and after Python clustering if you want to
#            select cells by cluster)
# Author:   Jonathon Sewell
# Created:  2025-09-15
# Version:  5   (v4 = rebuild of 02_1_SubsetAndExport_v3; v5 = adds the selection log, see C below)
# Edited:   2026-10-22
#
# Purpose:
# * Take the cells you want to carry into the NEXT round of clustering, choose
#   which channels the next round clusters on, and write a new data folder:
#     round N (clustered)  --02.1-->  round N+1 data (to be clustered)
# * Cells are selected by
#     - FCS file (metadata: Condition, TimePoint, ...): a table with Keep Y/N, and/or
#     - cluster / cluster group from the clustering you already did (Python),
#   and are KEPT (v3 asked what to EXCLUDE).
# * The kept cells' earlier cluster assignments are carried over, so Python can
#   continue with CLUSTER_ROUND = N + 1 (sub-clustering the groups you defined).
# * Also works before any clustering, to drop files you do not want clustered.
#
# Origin:
# * 03.1_SubsetAndExport.R - Jonathon Sewell, University of Virginia, Sept 15, 2025
#
# ==============================================================================
# CHANGE LOG: v3 -> v4 (rebuilt; the approach was agreed with you)
# ==============================================================================
# Convention:  [RESULT]   = can change the content of the output files
#              [BEHAVIOR] = changes how the script behaves, not the data
#
# A. Needs your attention
#  A1 [RESULT] KEEP instead of EXCLUDE. v3 asked which values to exclude and
#     combined variables with OR (a cell dropped if it matched ANY exclusion).
#     v4 keeps cells that are in a kept file AND in a kept cluster/group.
#  A2 [RESULT, bug fix] v3 did not run (syntax error `default =  =`, undefined
#     names such as METADATA_FILENAME, FILE_NUM, EXPRS_MAT_INFILE1, SUBSAMPLE_BOOL,
#     and `fread` without data.table), and its cluster logic had these problems:
#       - cluster column hard-coded to R1 for "Cluster" and "Group";
#       - "Group" used the ROW POSITION in the groups file as the group number
#         instead of the Group_Assigns value;
#       - the groups file was copied whole, so it listed clusters with no cells.
#     v4 uses the round you choose, the Group_Assigns value, and a groups file
#     that lists only the clusters still present.
#  A3 [RESULT] New Analysis channels. You choose the channels for the next round
#     in a table (starts from the current flags). The expression matrices are
#     re-split: expression_matrix_analysis.csv gets the new Analysis channels,
#     expression_matrix_other.csv the rest. Values are NOT recomputed (the
#     transform from 02_FCS_to_CSV is kept), so the numbers for a channel are
#     identical to the input folder's.
#  A4 [RESULT] metadata.csv keeps only files that still have cells, and Event_Count
#     is recomputed as the number of cells in the new matrix. File_Num is NOT
#     renumbered: filenums.csv still holds the original File_Num, so files stay
#     traceable across rounds.
#  A5 [BEHAVIOR] Row-aligned files are subset to the kept cells: all
#     cluster_R*_assigns.csv, umap_xy_R*.csv and umap_xyz_3D*.csv found in the input
#     folder. cluster_R*_groups.csv keep the clusters that still have cells (your
#     Group_Assigns and labels untouched). subsample_idxs.csv is not carried over
#     (it refers to the original files, not to the new rows); the input folder
#     keeps it.
#  A6 [BEHAVIOR] Output folder: <date>_<time>_Subset_R<N>-to-R<N+1>[_label] (or
#     _Subset_Files before any clustering), directly in the FCS folder, with the
#     input folder as parent in run_parameters.json. v3 wrote into a /tmp link.
#  A7 [BEHAVIOR] The run record section is named 02.1_SubsetAndExport_R<N>toR<N+1>
#     so repeating this step on later rounds adds to the history instead of
#     replacing the earlier entry.
#  A8 [BEHAVIOR] Clustering scripts for your compute environment are copied with
#     the same rules as 02_FCS_to_CSV v6. You must still edit CLUSTER_ROUND at the
#     top of the Python scripts to N + 1 (this script tells you).
#  A9 [BEHAVIOR] Checks: matrices, filenums and cluster files must have the same
#     number of rows; the matrix columns must equal panel.csv's Fixed_Param; at
#     least one cell and one Analysis channel must remain.
#
# C. v4 -> v5: selection log (ADDITIVE; no existing output file changes)
#  C1 [BEHAVIOR] The file table has a new editable column "Reason". Type why a file
#     is left out (e.g. "72h", "Univ", "Stock"). It has no effect on which cells are kept.
#  C2 [BEHAVIOR] New output file selection_log_files.csv in the new folder: one row
#     per file with its metadata, Cells_Before, Cells_After, Keep, Reason. This is the
#     record of the analysis-set rules; later scripts do not read it.
#  C3 [BEHAVIOR] The run record gains Files_Left_Out_Detail (file, metadata, cells,
#     reason) and Cells_By_Metadata (files and cells before/after for every metadata
#     value, e.g. TimePoint, SampleType) so the exclusions can be read at a glance.
#  C4 Which cells are kept is decided by exactly the same code as v4, so the
#     expression matrices, filenums, metadata and cluster files are unchanged.
#     Check: run v4 and v5 on the same choices and compare with Compare_GoldenFiles.R
#     (the only difference should be the new selection_log_files.csv and the run record).
#
# B. Efficiency / readability / reproducibility
#  B1 Shared helper file; all choices are typed answers or a browser table (no
#     XQuartz). Packages: svDialogs, rstudioapi, data.table, jsonlite, shiny,
#     miniUI, rhandsontable (dplyr no longer needed).
#  B2 Run record: input/output folders, rounds, which files/clusters/groups were
#     kept, cells before/after per file and per cluster, channels before/after,
#     checksums, session.
#  B3 No random numbers are used, so no seed is needed.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
# (none: everything is chosen in the pop-ups and tables)

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "02.1_SubsetAndExport.R"
SCRIPT_VERSION           <- 5
MODULE_BASE_NAME         <- "02.1_SubsetAndExport"
PREREQUISITE_MODULES     <- c("01_Get_Metadata", "01.2_Update_Metadata", "02_FCS_to_CSV")
REQUIRED_HELPERS_VERSION <- 6
EXPR_ANALYSIS_NAME       <- "expression_matrix_analysis.csv"
EXPR_OTHER_NAME          <- "expression_matrix_other.csv"
FILENUMS_NAME            <- "filenums.csv"



#### Load shared helper functions (do not edit this block) #####################
# Looks for ScriptsToUse/Helpers/00_Helpers.R starting at this script's folder
# and moving up to 4 folders.
local({
  start_dir <- tryCatch({
    p <- rstudioapi::getSourceEditorContext()$path
    if (is.character(p) && nzchar(p)) dirname(p) else getwd()
  }, error = function(e) getwd())
  search_dir   <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  helpers_file <- NA_character_
  for (i in 1:5) {
    candidates <- c(file.path(search_dir, "Helpers", "00_Helpers.R"),
                    file.path(search_dir, "00_Helpers.R"))
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates) > 0) { helpers_file <- candidates[1]; break }
    search_dir <- dirname(search_dir)
  }
  if (is.na(helpers_file)) {
    stop("Could not find 00_Helpers.R.\n",
         "It should be in a folder called 'Helpers' inside ScriptsToUse.\n",
         "Looked upward from: ", start_dir,
         "\n(Is this script saved inside your ScriptsToUse folder, and is the ",
         "'rstudioapi' package installed? RStudio needs it to find this script.)",
         call. = FALSE)
  }
  assign("HELPERS_DIR", dirname(helpers_file), envir = .GlobalEnv)
  source(helpers_file, local = FALSE)
})
check_helpers_version(REQUIRED_HELPERS_VERSION)



#### Set up environment ########################################################
pkg_versions <- dc_load_packages(
  cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite",
           "shiny", "miniUI", "rhandsontable"))

cat("\n========================================================\n")
cat("[START] Subset expression matrices and supporting files:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

Original_WD <- get_working_dir()   # the folder with your original FCS files



#### Functions specific to this script #########################################
ask_text <- function(message, default = "") {
  res <- svDialogs::dlg_input(message = message, default = default)$res
  if (length(res) != 1 || is.na(res)) {
    stop("The pop-up window was closed without an answer. Please run the script again.",
         call. = FALSE)
  }
  trimws(res)
}

# Ask until `parse_fun(text)` returns list(ok = TRUE, value = ...); 3 tries
ask_until <- function(message, default, parse_fun) {
  for (attempt in 1:3) {
    parsed <- parse_fun(ask_text(message, default))
    if (isTRUE(parsed$ok)) return(parsed$value)
    message <- paste0(parsed$msg, "\n\n", message)
  }
  stop("No valid answer was entered. Please run the script again.", call. = FALSE)
}

# "1,3,5-7" -> c(1,3,5,6,7); NULL if it cannot be read
parse_number_list <- function(txt) {
  parts <- strsplit(gsub("\\s", "", txt), ",")[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(NULL)
  out <- numeric(0)
  for (p in parts) {
    if (grepl("^[0-9]+$", p)) {
      out <- c(out, as.numeric(p))
    } else if (grepl("^[0-9]+-[0-9]+$", p)) {
      lims <- as.numeric(strsplit(p, "-")[[1]])
      if (lims[1] > lims[2]) return(NULL)
      out <- c(out, lims[1]:lims[2])
    } else return(NULL)
  }
  sort(unique(out))
}

# Generic browser table. df: all columns text. readonly: columns that cannot be
# edited; dropdowns: named list column -> allowed values. Returns the edited
# table, or NULL if cancelled.
run_table_editor <- function(df, title, help_text, readonly, dropdowns) {
  ui <- miniPage(
    gadgetTitleBar(title),
    miniContentPanel(helpText(help_text),
                     rHandsontableOutput("grid", height = "100%"))
  )
  server <- function(input, output, session) {
    output$grid <- renderRHandsontable({
      d   <- if (!is.null(input$grid)) hot_to_r(input$grid) else df
      hot <- rhandsontable(d, rowHeaders = FALSE, height = 550, overflow = "visible",
                           columnSorting = TRUE, manualColumnResize = TRUE)
      hot <- hot_col(hot, readonly, readOnly = TRUE)
      for (nm in names(dropdowns)) {
        hot <- hot_col(hot, nm, type = "dropdown", source = dropdowns[[nm]], strict = TRUE)
      }
      hot_table(hot, highlightCol = TRUE, highlightRow = TRUE, fontSize = 11)
    })
    observeEvent(input$done, {
      final <- if (!is.null(input$grid)) hot_to_r(input$grid) else df
      stopApp(list(table = final))
    })
  }
  res <- tryCatch(runGadget(ui, server, viewer = browserViewer()),
                  error = function(e) {
                    if (grepl("cancel", conditionMessage(e), ignore.case = TRUE)) NULL else stop(e)
                  })
  if (is.null(res)) NULL else res$table
}



#### Choose the input folder ###################################################
dirs <- list.dirs(Original_WD, recursive = FALSE, full.names = TRUE)
dirs <- dirs[!grepl("^[_.]", basename(dirs))]
has_data <- vapply(dirs, function(d) file.exists(file.path(d, EXPR_ANALYSIS_NAME)) &&
                     file.exists(file.path(d, RUN_RECORD_NAME)), logical(1))
data_dirs <- sort(dirs[has_data], decreasing = TRUE)      # newest first
if (length(data_dirs) == 0) {
  stop("No folder with ", EXPR_ANALYSIS_NAME, " was found in:\n  ", Original_WD,
       "\nRun 02_FCS_to_CSV first.", call. = FALSE)
}
mem_dir <- get0("Compiled_WD", envir = .GlobalEnv, inherits = FALSE)
default_choice <- if (!is.null(mem_dir) && basename(mem_dir) %in% basename(data_dirs)) basename(mem_dir) else basename(data_dirs[1])
choice <- ask_choice("Which data folder should be subset?", basename(data_dirs), preselect = default_choice)
input_dir <- normalizePath(data_dirs[basename(data_dirs) == choice], winslash = "/", mustWork = TRUE)
cat("[INFO] Input folder:\n       ", input_dir, "\n", sep = "")

run_record <- require_modules(input_dir, PREREQUISITE_MODULES,
                              how_to_fix = "Run 02_FCS_to_CSV first.")



#### Round of clustering to select on ##########################################
assign_files <- list.files(input_dir, pattern = "^cluster_R[0-9]+_assigns\\.csv$")
rounds_found <- sort(as.integer(sub("^cluster_R([0-9]+)_assigns\\.csv$", "\\1", assign_files)))
if (length(rounds_found) == 0) {
  cat("[INFO] No cluster_R*_assigns.csv files in this folder: files can be selected, clusters cannot.\n")
  round_in <- 0L
} else {
  round_in <- ask_until(
    paste0("Cluster round to select cells from.\n\nFound: ",
           paste0("R", rounds_found, collapse = ", "),
           "\nType the number (e.g. 1), or 0 to ignore the cluster files."),
    as.character(max(rounds_found)),
    function(x) {
      v <- suppressWarnings(as.integer(x))
      if (!is.na(v) && (v == 0 || v %in% rounds_found)) list(ok = TRUE, value = v)
      else list(ok = FALSE, msg = "Not one of the rounds found.")
    })
}
use_clusters <- round_in >= 1
round_out    <- if (use_clusters) round_in + 1L else NA_integer_
Rcol         <- paste0("R", round_in)



#### Read the small files ######################################################
metadata <- read_metadata(input_dir)
panel    <- read_panel(input_dir)
filenums <- data.table::fread(file.path(input_dir, FILENUMS_NAME), data.table = FALSE)[[1]]
n_cells_in <- length(filenums)

if (anyDuplicated(panel$Fixed_Param)) stop("Fixed_Param is not unique in panel.csv.", call. = FALSE)
if (!all(unique(filenums) %in% metadata$File_Num)) {
  stop("filenums.csv contains File_Num values that are not in metadata.csv.", call. = FALSE)
}

if (use_clusters) {
  assigns_path <- file.path(input_dir, sprintf("cluster_R%d_assigns.csv", round_in))
  groups_path  <- file.path(input_dir, sprintf("cluster_R%d_groups.csv", round_in))
  assigns_in <- data.table::fread(assigns_path, data.table = FALSE)
  if (!Rcol %in% colnames(assigns_in)) {
    stop(basename(assigns_path), " has no column named ", Rcol, ".", call. = FALSE)
  }
  if (nrow(assigns_in) != n_cells_in) {
    stop(basename(assigns_path), " has ", nrow(assigns_in), " rows but filenums.csv has ",
         n_cells_in, ". The cluster file belongs to a different dataset.", call. = FALSE)
  }
  cell_cluster <- assigns_in[[Rcol]]
  groups_in <- if (file.exists(groups_path)) {
    data.table::fread(groups_path, data.table = FALSE, colClasses = "character")
  } else NULL
}



#### 1. Choose files (browser table) ###########################################
cells_per_file_in <- table(factor(filenums, levels = metadata$File_Num))
meta_cols <- intersect(c("SampleCollection", "Expt_ID", "BarcodeSet", "Treatment",
                         "Condition", "TimePoint", "Replicate", "SampleType"),
                       colnames(metadata))
present <- metadata$File_Num[cells_per_file_in[as.character(metadata$File_Num)] > 0]
file_tbl <- metadata[metadata$File_Num %in% present, c("File_Num", "File_Name", meta_cols), drop = FALSE]
file_tbl$Cells <- as.integer(cells_per_file_in[as.character(file_tbl$File_Num)])
file_tbl$Keep  <- "Y"
file_tbl$Reason <- ""
file_tbl[] <- lapply(file_tbl, as.character)
rownames(file_tbl) <- NULL

keep_files <- file_tbl$File_Num
file_reason <- setNames(rep("", nrow(file_tbl)), file_tbl$File_Num)
file_selection <- "all files kept"
if (ask_yes_no(paste0("Select cells by FCS file (", nrow(file_tbl), " files)?\n\n",
                      "Yes = a table opens in your browser: set Keep to N for files to leave out\n",
                      "No  = keep all files"))) {
  cat("\n[INFO] A table opens in your web browser. Set Keep (Y/N), then click 'Done'.\n")
  edited <- run_table_editor(
    file_tbl, "Subset: choose files to keep",
    "Set Keep to Y for files whose cells go into the next round, N to leave them out. Optionally type why in Reason (e.g. 72h, Univ, Stock): it is only recorded. You can sort by clicking a column header.",
    readonly = setdiff(colnames(file_tbl), c("Keep", "Reason")), dropdowns = list(Keep = c("Y", "N")))
  if (is.null(edited)) stop("Cancelled. No files were written.", call. = FALSE)
  edited <- edited[match(file_tbl$File_Num, edited$File_Num), , drop = FALSE]
  keep_files <- file_tbl$File_Num[edited$Keep == "Y"]
  file_reason[] <- ifelse(is.na(edited$Reason), "", trimws(as.character(edited$Reason)))
  file_selection <- "chosen in table"
}
keep_files <- as.integer(keep_files)
if (length(keep_files) == 0) stop("No files were kept.", call. = FALSE)



#### 2. Choose clusters / groups (typed) #######################################
cluster_selection <- list(Mode = "not used")
kept_clusters <- NULL
if (use_clusters) {
  cl_tab <- table(cell_cluster)
  cat("\n[INFO] Cells per cluster (", Rcol, ", all files):\n", sep = "")
  print(cl_tab)

  grp_num <- NULL
  if (!is.null(groups_in) && "Group_Assigns" %in% colnames(groups_in)) {
    grp_num <- suppressWarnings(as.numeric(groups_in$Group_Assigns))
    if (all(is.na(grp_num))) grp_num <- NULL
  }
  chosen <- "clusters"
  if (!is.null(grp_num)) {
    grp_df <- data.frame(Group_Assigns = grp_num, Cluster = groups_in[[Rcol]],
                         stringsAsFactors = FALSE)
    cat("\n[INFO] Cluster groups (Group_Assigns; 0 = carried through without sub-clustering):\n")
    print(aggregate(Cluster ~ Group_Assigns, grp_df, function(x) paste(x, collapse = ",")), row.names = FALSE)
    grp_choice <- ask_until(
      paste0("Keep which cluster GROUPS (Group_Assigns)?\n\n",
             "Type group numbers separated by commas (e.g. 1,3)\n",
             "'all' = keep every cluster\n",
             "'clusters' = choose individual clusters instead"),
      "all",
      function(x) {
        if (tolower(x) %in% c("all", "clusters")) return(list(ok = TRUE, value = tolower(x)))
        v <- parse_number_list(x)
        if (is.null(v)) list(ok = FALSE, msg = "Could not read that. Example: 1,3")
        else if (!all(v %in% grp_num)) list(ok = FALSE, msg = paste0("Not in the groups file: ", paste(setdiff(v, grp_num), collapse = ", ")))
        else list(ok = TRUE, value = v)
      })
    if (is.numeric(grp_choice)) {
      kept_clusters <- grp_df$Cluster[grp_df$Group_Assigns %in% grp_choice]
      cluster_selection <- list(Mode = "groups", Groups_Kept = as.list(grp_choice))
      chosen <- "groups"
    } else if (identical(grp_choice, "all")) {
      kept_clusters <- names(cl_tab)
      cluster_selection <- list(Mode = "all clusters")
      chosen <- "all"
    }
  } else {
    cat("[INFO] Group_Assigns is empty or missing in the groups file, so clusters are chosen individually.\n")
  }
  if (chosen == "clusters") {
    cl_choice <- ask_until(
      paste0("Keep which clusters of ", Rcol, "?\n\nType cluster numbers (e.g. 1,4,7-9) or 'all'."),
      "all",
      function(x) {
        if (tolower(x) == "all") return(list(ok = TRUE, value = "all"))
        v <- parse_number_list(x)
        if (is.null(v)) list(ok = FALSE, msg = "Could not read that. Example: 1,4,7-9")
        else if (!all(v %in% as.numeric(names(cl_tab)))) list(ok = FALSE, msg = paste0("No such cluster: ", paste(setdiff(v, as.numeric(names(cl_tab))), collapse = ", ")))
        else list(ok = TRUE, value = v)
      })
    if (is.numeric(cl_choice)) {
      kept_clusters <- as.character(cl_choice)
      cluster_selection <- list(Mode = "clusters", Clusters_Kept = as.list(cl_choice))
    } else {
      kept_clusters <- names(cl_tab)
      cluster_selection <- list(Mode = "all clusters")
    }
  }
}



#### Which cells are kept, and confirm #########################################
keep_cell <- filenums %in% keep_files
if (use_clusters) keep_cell <- keep_cell & (as.character(cell_cluster) %in% as.character(kept_clusters))
n_cells_out <- sum(keep_cell)
if (n_cells_out == 0) stop("No cells match your selection.", call. = FALSE)

cells_per_file_out <- table(factor(filenums[keep_cell], levels = metadata$File_Num))
summary_files <- data.frame(File_Num = metadata$File_Num[metadata$File_Num %in% present],
                            Before = as.integer(cells_per_file_in[as.character(present)]),
                            After  = as.integer(cells_per_file_out[as.character(present)]))
# Selection log (v5): per-file table and per-metadata-value summary
sel_log <- metadata[match(summary_files$File_Num, metadata$File_Num), c("File_Num", "File_Name", meta_cols), drop = FALSE]
sel_log$Cells_Before <- summary_files$Before
sel_log$Cells_After  <- summary_files$After
sel_log$Keep         <- ifelse(sel_log$File_Num %in% keep_files, "Y", "N")
sel_log$Reason       <- unname(file_reason[as.character(sel_log$File_Num)])
sel_log$Reason[is.na(sel_log$Reason)] <- ""
cells_by_metadata <- do.call(c, lapply(meta_cols, function(mc) {
  lv <- as.character(sel_log[[mc]]); lv[is.na(lv)] <- "NA"
  lapply(unique(lv), function(v) {
    i <- lv == v
    list(Column = mc, Value = v, Files_Before = sum(sel_log$Cells_Before[i] > 0),
         Files_After = sum(sel_log$Cells_After[i] > 0),
         Cells_Before = sum(sel_log$Cells_Before[i]), Cells_After = sum(sel_log$Cells_After[i]))
  })
}))
cat("\n[CHECK] Cells before / after: ", n_cells_in, " / ", n_cells_out,
    "   Files with cells: ", sum(summary_files$Before > 0), " / ", sum(summary_files$After > 0), "\n", sep = "")
if (n_cells_out < 1000) {
  cat("[WARNING] Only ", n_cells_out, " cells remain; clustering may not be meaningful.\n", sep = "")
}
if (!ask_yes_no(paste0("Keep ", n_cells_out, " of ", n_cells_in, " cells (",
                       round(100 * n_cells_out / n_cells_in, 1), "%) from ",
                       sum(summary_files$After > 0), " files?\n\n",
                       "Files: ", file_selection, "\n",
                       "Clusters: ", if (use_clusters) paste0(cluster_selection$Mode, " (", Rcol, ")") else "not used",
                       "\n\nContinue?"))) {
  stop("Cancelled. Nothing was written.", call. = FALSE)
}



#### 3. Choose Analysis channels for the next round (browser table) ############
panel_edit <- panel[, c("Name", "Fixed_Param", "Analysis", "Plot"), drop = FALSE]
panel_edit[] <- lapply(panel_edit, as.character)
panel_new <- panel
repeat {
  cat("\n[INFO] A table opens in your web browser. Set Analysis = 1 for the channels the next\n",
      "       round clusters on (Plot = 1 for channels shown in plots), then click 'Done'.\n", sep = "")
  edited <- run_table_editor(
    panel_edit, "Subset: choose Analysis channels for the next round",
    "Analysis = 1: used for clustering. Plot = 1: shown in plots. Name and Fixed_Param cannot be changed.",
    readonly = c("Name", "Fixed_Param"), dropdowns = list(Analysis = c("0", "1"), Plot = c("0", "1")))
  if (is.null(edited)) stop("Cancelled. Nothing was written.", call. = FALSE)
  edited <- edited[match(panel$Name, edited$Name), , drop = FALSE]
  if (sum(edited$Analysis == "1") == 0) {
    svDialogs::dlg_message("At least one channel must have Analysis = 1.", type = "ok")
    panel_edit <- edited
    next
  }
  panel_new$Analysis <- as.integer(edited$Analysis)
  panel_new$Plot     <- as.integer(edited$Plot)
  break
}
analysis_before <- panel$Fixed_Param[panel$Analysis == 1]
analysis_after  <- panel_new$Fixed_Param[panel_new$Analysis == 1]



#### Label and compute environment #############################################
label_txt <- ask_text("Optional label for the new folder (e.g. Neurons). Leave blank for none.", "")
label_txt <- gsub("[^A-Za-z0-9._+-]+", "-", label_txt)
cluster_env <- ask_choice("Will you be clustering locally or on an HPC?", c("locally", "HPC"), preselect = "locally")



#### Read the expression matrices and subset ###################################
cat("\n[INFO] Reading expression matrices...\n")
m_analysis <- data.table::fread(file.path(input_dir, EXPR_ANALYSIS_NAME), data.table = FALSE)
m_other    <- data.table::fread(file.path(input_dir, EXPR_OTHER_NAME),    data.table = FALSE)
if (nrow(m_analysis) != n_cells_in || nrow(m_other) != n_cells_in) {
  stop("The expression matrices and filenums.csv do not have the same number of rows.", call. = FALSE)
}
if (!identical(colnames(m_analysis), analysis_before) ||
    !identical(colnames(m_other), panel$Fixed_Param[panel$Analysis != 1])) {
  stop("The columns of the expression matrices do not match panel.csv in the input folder ",
       "(Analysis = 1 channels, then the other channels, in panel order).", call. = FALSE)
}
both <- cbind(m_analysis[keep_cell, , drop = FALSE], m_other[keep_cell, , drop = FALSE])
rm(m_analysis, m_other)
out_analysis <- both[, analysis_after, drop = FALSE]
out_other    <- both[, panel_new$Fixed_Param[panel_new$Analysis != 1], drop = FALSE]
rm(both)



#### Write the new folder ######################################################
stage_label <- paste0("Subset_",
                      if (use_clusters) sprintf("R%d-to-R%d", round_in, round_out) else "Files",
                      if (nzchar(label_txt)) paste0("_", label_txt))
Subset_WD <- make_stage_dir(Original_WD, stage_label)
cat("\n[INFO] Writing files to:\n       ", Subset_WD, "\n", sep = "")

data.table::fwrite(out_analysis, file.path(Subset_WD, EXPR_ANALYSIS_NAME))
data.table::fwrite(out_other,    file.path(Subset_WD, EXPR_OTHER_NAME))
data.table::fwrite(matrix(filenums[keep_cell], ncol = 1), file.path(Subset_WD, FILENUMS_NAME))

# metadata: files that still have cells; Event_Count = cells in the new matrix
meta_out <- metadata[metadata$File_Num %in% unique(filenums[keep_cell]), , drop = FALSE]
meta_out$Event_Count <- as.integer(cells_per_file_out[as.character(meta_out$File_Num)])
data.table::fwrite(meta_out,  file.path(Subset_WD, METADATA_NAME))
data.table::fwrite(panel_new, file.path(Subset_WD, PANEL_NAME))
data.table::fwrite(sel_log,   file.path(Subset_WD, "selection_log_files.csv"))

# Row-aligned files from earlier rounds
rowwise_written <- character(0); rowwise_skipped <- character(0)
rowwise_files <- list.files(input_dir, pattern = "^(cluster_R[0-9]+_assigns|umap_xy_R[0-9]+|umap_xyz_3D(_R[0-9]+)?)\\.csv$")
for (f in rowwise_files) {
  d <- data.table::fread(file.path(input_dir, f), data.table = FALSE)
  if (nrow(d) == n_cells_in) {
    data.table::fwrite(d[keep_cell, , drop = FALSE], file.path(Subset_WD, f))
    rowwise_written <- c(rowwise_written, f)
  } else rowwise_skipped <- c(rowwise_skipped, f)
}
# Groups files: keep the clusters that still have cells
groups_written <- character(0)
for (f in list.files(input_dir, pattern = "^cluster_R[0-9]+_groups\\.csv$")) {
  k <- sub("^cluster_(R[0-9]+)_groups\\.csv$", "\\1", f)
  g <- data.table::fread(file.path(input_dir, f), data.table = FALSE, colClasses = "character")
  a_file <- sprintf("cluster_%s_assigns.csv", k)
  if (a_file %in% rowwise_written && k %in% colnames(g)) {
    a <- data.table::fread(file.path(Subset_WD, a_file), data.table = FALSE)
    g <- g[g[[k]] %in% as.character(unique(a[[k]])), , drop = FALSE]
  }
  data.table::fwrite(g, file.path(Subset_WD, f))
  groups_written <- c(groups_written, f)
}

script_res <- copy_clustering_scripts(Subset_WD, cluster_env)



#### Keep objects in memory for the next scripts ###############################
Compiled_WD <- Subset_WD
# Also in memory: Original_WD (unchanged)



#### Record parameters used in this run ########################################
module_name <- paste0(MODULE_BASE_NAME, if (use_clusters) sprintf("_R%dtoR%d", round_in, round_out) else "_Files")
if (module_name %in% names(run_record)) {            # same step repeated: keep both entries
  k <- 2; while (paste0(module_name, "_", k) %in% names(run_record)) k <- k + 1
  module_name <- paste0(module_name, "_", k)
}
cluster_summary <- if (use_clusters) {
  kept_tab <- table(cell_cluster[keep_cell])
  lapply(names(cl_tab), function(cl) list(Cluster = cl, Cells_Before = as.integer(cl_tab[[cl]]),
                                          Cells_After = if (cl %in% names(kept_tab)) as.integer(kept_tab[[cl]]) else 0L))
} else NULL
csv_out <- file.path(Subset_WD, c(EXPR_ANALYSIS_NAME, EXPR_OTHER_NAME, FILENUMS_NAME, METADATA_NAME, PANEL_NAME, "selection_log_files.csv"))
Inputs <- list(
  Script_Version        = SCRIPT_VERSION,
  Completed_At          = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  Original_WD           = Original_WD,
  INPUT_FOLDER          = input_dir,
  OUTPUT_FOLDER         = Subset_WD,
  Input_Folder_Name     = basename(input_dir),
  Round_Selected_On     = if (use_clusters) round_in else "none",
  Round_For_Python      = if (use_clusters) round_out else "unchanged (no cluster files used)",
  File_Selection        = file_selection,
  Files_Kept            = as.list(sort(keep_files)),
  Files_Left_Out        = as.list(setdiff(present, keep_files)),
  Files_Left_Out_Detail = lapply(which(sel_log$Keep == "N"), function(i) as.list(sel_log[i, , drop = FALSE])),
  Cells_By_Metadata     = cells_by_metadata,
  Cluster_Selection     = cluster_selection,
  N_Cells_Before        = n_cells_in,
  N_Cells_After         = n_cells_out,
  Cells_Per_File        = lapply(seq_len(nrow(summary_files)), function(i)
    list(File_Num = summary_files$File_Num[i], Before = summary_files$Before[i], After = summary_files$After[i])),
  Cells_Per_Cluster     = cluster_summary,
  Analysis_Channels_Before = as.list(analysis_before),
  Analysis_Channels_After  = as.list(analysis_after),
  Plot_Channels_After      = as.list(panel_new$Fixed_Param[panel_new$Plot == 1]),
  Files_Subset_By_Rows     = as.list(rowwise_written),
  Groups_Files_Filtered    = as.list(groups_written),
  Files_Not_Carried_Over   = as.list(c(rowwise_skipped, if (file.exists(file.path(input_dir, "subsample_idxs.csv"))) "subsample_idxs.csv")),
  ComputeEnvironment       = cluster_env,
  Scripts_Copied           = as.list(script_res$copied),
  Scripts_Not_Found        = as.list(script_res$not_found),
  Random_Seed              = "not used in this script",
  Output_MD5               = file_md5(csv_out),
  Session                  = session_info_list(names(pkg_versions))
)
record_params(module_name, Inputs, write_dir = Subset_WD, parent_dir = input_dir)



#### Script completed message ##################################################
next_steps <- c(paste0("Files are in:  ", basename(Subset_WD)),
                if (use_clusters) paste0("Python: set CLUSTER_ROUND = ", round_out,
                                         " at the top of the clustering scripts (it sub-clusters the groups in cluster_R",
                                         round_in, "_groups.csv)")
                else "Python: CLUSTER_ROUND stays at 1.")
print_banner(" Subset complete ", next_steps)
