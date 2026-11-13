################################################################################
# 01.2_Update_Metadata_v5.R
################################################################################
# Layer:    DataCleanup -- SECOND script (run after 01_Get_Metadata)
# Author:   Jonathon Sewell
# Created:  2026-07-22
# Version:  5   (input version: 01.2_Update_Metadata_v4)
# Edited:   2026-10-21
#
# Purpose:
# * Open metadata.csv and panel.csv (made by 01_Get_Metadata) in an editor in
#   your web browser, so you can check and complete them
#     - metadata.csv : fill in the user-defined columns (Expt_ID, SampleType, ...)
#     - panel.csv    : choose the channels for clustering and plotting
#         'Analysis' = 1  -> channel is used for CLUSTERING
#         'Plot'     = 1  -> channel is included in PLOTS
#         (Transform / Factor control the asinh transformation in 02_FCS_to_CSV)
# * Save the edited files and record what you chose in run_parameters.json
#
# Origin: N/A
# * Previously required the user to edit the files in Excel
#
# ==============================================================================
# CHANGE LOG: v4 -> v5
# ==============================================================================
# Convention:  [RESULT]   = can change the content/order of the saved CSV files
#              [BEHAVIOR] = changes how the script behaves, not the data
#              [SAME]     = output unchanged
#
# A. Needs your attention
#  A1 [RESULT] Row order is restored before saving. The editor lets you sort by
#     clicking column headers, and v4 saved the rows in whatever order they were
#     on screen. Now: metadata.csv is saved in File_Num order, and panel.csv in
#     the ORIGINAL channel order from 01_Get_Metadata. This matters because
#     02_FCS_to_CSV names the data columns by panel row position -- a panel that
#     was sorted on screen would silently attach the wrong marker names to the
#     data. If you never sorted, nothing changes.
#  A2 [BEHAVIOR] The files on disk are always the starting point. v4 used the
#     copies in memory first, so edits made in Excel between scripts were
#     ignored. Now any edit you make outside R is picked up.
#  A3 [BEHAVIOR] Key columns cannot be edited: File_Name, File_Num, Event_Count
#     (metadata) and Name (panel). Changing them would break the link to the FCS
#     files. (To undo: delete the two readOnly lines in run_editor().)
#  A4 [BEHAVIOR] Checks before saving. Problems that must be fixed (editor reopens
#     with your changes): no Analysis channel; Analysis/Plot not 0 or 1; Transform
#     not Asinh/Linear; Factor missing or <= 0; Fixed_Param blank or duplicated.
#     Reminders you may override: SampleType not Expt/Ctrl, Expt_ID still
#     "User-Defined #" or blank, Condition/Treatment blank, no Plot channel.
#  A5 [SAME, format] Files are written with fwrite() like 01_Get_Metadata, so text
#     is no longer wrapped in quotes unless needed. Values are identical.
#  A6 [BEHAVIOR] A timestamped copy of the files is saved in _backups before they
#     are overwritten (the first copy is the untouched output of 01).
#
# B. Efficiency / readability / reproducibility
#  B1 Shared helper file (00_Helpers.R) loaded by a short block at the top.
#  B2 Sequence check: stops with a plain message if 01_Get_Metadata has not been
#     run in this folder (looks in run_parameters.json). Later scripts will use
#     the 'Reviewed' entry this script writes to check that 01.2 was run (this
#     will replace the "Has 01.2 been run?" question in 02_RunBatch).
#  B3 No automatic install.packages(); missing packages give install instructions.
#  B4 dplyr and the %>% pipe are no longer needed; code is plain stepwise calls.
#  B5 Panel editor: Transform is a dropdown (Asinh / Linear).
#  B6 Plain-language help text above each table; the v4 header described 'Analysis'
#     twice (the second one meant 'Plot') -- corrected.
#  B7 run_parameters.json now records WHAT YOU CHOSE (channels for Analysis and
#     Plot, SampleType counts, which columns you edited) and checksums of the
#     saved files, plus a 'Reviewed' time stamp. The v4 keys keep their names
#     (WRITE_PATH was removed with the project-wide change).
#  B8 Cancelling the editor stops cleanly with nothing saved.
################################################################################


#### USER SETTINGS -- edit this section only ###################################
# (none: everything is chosen in the editor)

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "01.2_Update_Metadata_v5.R"
SCRIPT_VERSION           <- 5
MODULE_NAME              <- "01.2_Update_Metadata"
PREREQUISITE_MODULE      <- "01_Get_Metadata"
REQUIRED_HELPERS_VERSION <- 1
SAMPLETYPE_CHOICES       <- c("Expt", "Ctrl")
TRANSFORM_CHOICES        <- c("Asinh", "Linear")



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
cat("[START] Modify files necessary for data clean up:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

Original_WD <- get_working_dir()

# Sequence check: 01_Get_Metadata must have been run in this folder
require_modules(Original_WD, PREREQUISITE_MODULE,
                how_to_fix = "Run 01_Get_Metadata first.")



#### Read the files from disk ##################################################
metadata_in <- read_metadata(Original_WD)
panel_in    <- read_panel(Original_WD)

if (!"SampleType" %in% colnames(metadata_in)) metadata_in$SampleType <- "Expt"
if (anyDuplicated(metadata_in$File_Num) || anyDuplicated(metadata_in$File_Name)) {
  stop("metadata.csv has duplicated File_Num or File_Name values. ",
       "Re-run 01_Get_Metadata.", call. = FALSE)
}
if (anyDuplicated(panel_in$Name)) {
  stop("panel.csv has duplicated channel Names. Re-run 01_Get_Metadata.", call. = FALSE)
}



#### Functions specific to this script #########################################
# Shiny gadget: two editable tables (metadata, panel) in the web browser.
# Returns list(metadata = ..., panel = ...) or NULL if cancelled.
run_editor <- function(md, pn) {
  pn$Analysis <- as.numeric(pn$Analysis)
  pn$Plot     <- as.numeric(pn$Plot)

  ui <- miniPage(
    gadgetTitleBar("Data Cleanup: Metadata & Panel Annotator"),
    miniTabstripPanel(

      # TAB 1: Edit metadata
      miniTabPanel("1. Metadata Editor", icon = icon("table"),
                   miniContentPanel(
                     helpText("Click any white cell to edit. Grey cells (File_Name, File_Num, Event_Count) cannot be changed. ",
                              "Fill in Expt_ID and SampleType (Expt or Ctrl) and check SampleCollection, Condition, TimePoint and Replicate. ",
                              "Click a column header to sort; the saved file keeps the original order."),
                     rHandsontableOutput("metadata_grid", height = "100%")
                   )
      ),

      # TAB 2: Edit panel annotations
      miniTabPanel("2. Panel Editor", icon = icon("sliders"),
                   miniContentPanel(
                     helpText("Analysis = 1: channel is used for CLUSTERING.  Plot = 1: channel is included in PLOTS.  ",
                              "Transform and Factor set the asinh transformation (Linear channels use Factor 1).  ",
                              "Fixed_Param is the marker name used downstream and must be unique.  ",
                              "Click a column header to sort; the saved file keeps the original channel order."),
                     rHandsontableOutput("panel_grid", height = "100%")
                   )
      )
    )
  )

  server <- function(input, output, session) {

    # 1. Metadata table
    output$metadata_grid <- renderRHandsontable({
      df  <- if (!is.null(input$metadata_grid)) hot_to_r(input$metadata_grid) else md
      hot <- rhandsontable(df,
                           rowHeaders = FALSE,
                           height = 550,
                           overflow = "visible",
                           columnSorting = TRUE,       # click headers to sort
                           manualColumnResize = TRUE)  # drag to resize columns
      hot <- hot_col(hot, c("File_Name", "File_Num", "Event_Count"), readOnly = TRUE)
      hot <- hot_col(hot, "SampleType", type = "dropdown",
                     source = SAMPLETYPE_CHOICES, strict = TRUE)
      hot_table(hot, highlightCol = TRUE, highlightRow = TRUE, fontSize = 11)
    })

    # 2. Panel table
    output$panel_grid <- renderRHandsontable({
      df  <- if (!is.null(input$panel_grid)) hot_to_r(input$panel_grid) else pn
      hot <- rhandsontable(df,
                           rowHeaders = FALSE,
                           height = 550,
                           overflow = "visible",
                           columnSorting = TRUE,
                           manualColumnResize = TRUE)
      hot <- hot_col(hot, "Name", readOnly = TRUE)
      hot <- hot_col(hot, "Analysis", type = "dropdown", source = c(0, 1), strict = TRUE)
      hot <- hot_col(hot, "Plot",     type = "dropdown", source = c(0, 1), strict = TRUE)
      hot <- hot_col(hot, "Transform", type = "dropdown",
                     source = TRANSFORM_CHOICES, strict = TRUE)
      hot_table(hot, highlightCol = TRUE, highlightRow = TRUE, fontSize = 11)
    })

    # Save and exit
    observeEvent(input$done, {
      final_metadata <- if (!is.null(input$metadata_grid)) hot_to_r(input$metadata_grid) else md
      final_panel    <- if (!is.null(input$panel_grid))    hot_to_r(input$panel_grid)    else pn
      stopApp(list(metadata = final_metadata, panel = final_panel))
    })
  }

  # Cancel -> NULL; any other error is passed on
  tryCatch(runGadget(ui, server, viewer = browserViewer()),
           error = function(e) {
             if (grepl("cancel", conditionMessage(e), ignore.case = TRUE)) NULL else stop(e)
           })
}

# Check the edited tables. Returns list(errors = ..., warnings = ...):
#   errors   = must be fixed before saving
#   warnings = reminders the user may override
find_edit_problems <- function(md, pn) {
  errors   <- character(0)
  warnings <- character(0)

  # ---- metadata ----
  stype <- as.character(md$SampleType)
  bad_stype <- is.na(stype) | !(stype %in% SAMPLETYPE_CHOICES)
  if (any(bad_stype)) {
    warnings <- c(warnings, paste0(
      "SampleType should be Expt or Ctrl, but ", sum(bad_stype),
      " file(s) have: ", paste(unique(ifelse(is.na(stype[bad_stype]), "(blank)",
                                             stype[bad_stype])), collapse = ", ")))
  }
  expt <- as.character(md$Expt_ID)
  expt[is.na(expt)] <- ""
  bad_expt <- grepl("User-Defined", expt, fixed = TRUE) | trimws(expt) == ""
  if (any(bad_expt)) {
    warnings <- c(warnings, paste0(
      "Expt_ID is still a placeholder or blank for ", sum(bad_expt), " file(s)."))
  }
  for (cn in c("Condition", "Treatment")) {
    v <- as.character(md[[cn]])
    if (any(is.na(v) | trimws(v) == "")) {
      warnings <- c(warnings, paste0(cn, " is blank for some files."))
    }
  }

  # ---- panel ----
  an <- suppressWarnings(as.numeric(pn$Analysis))
  pl <- suppressWarnings(as.numeric(pn$Plot))
  fa <- suppressWarnings(as.numeric(pn$Factor))
  if (any(is.na(an) | !(an %in% c(0, 1))) || any(is.na(pl) | !(pl %in% c(0, 1)))) {
    errors <- c(errors, "Analysis and Plot must be 0 or 1 for every channel.")
  }
  if (!any(an == 1, na.rm = TRUE)) {
    errors <- c(errors, "No channel has Analysis = 1. Choose the channels to cluster on.")
  }
  if (!any(pl == 1, na.rm = TRUE)) {
    warnings <- c(warnings, "No channel has Plot = 1.")
  }
  if (any(is.na(fa) | fa <= 0)) {
    errors <- c(errors, "Factor must be a number greater than 0 for every channel.")
  }
  if (any(is.na(pn$Transform) | !(pn$Transform %in% TRANSFORM_CHOICES))) {
    errors <- c(errors, "Transform must be Asinh or Linear for every channel.")
  }
  fp <- as.character(pn$Fixed_Param)
  if (any(is.na(fp) | trimws(fp) == "")) {
    errors <- c(errors, "Fixed_Param is blank for some channels.")
  }
  if (anyDuplicated(fp[!is.na(fp)])) {
    errors <- c(errors, paste0("Fixed_Param must be unique. Duplicated: ",
                               paste(unique(fp[duplicated(fp)]), collapse = ", ")))
  }
  list(errors = errors, warnings = warnings)
}

# Names of the columns whose values differ between `old` and `new`
# (rows matched on `key`)
changed_columns <- function(old, new, key) {
  old <- old[match(new[[key]], old[[key]]), , drop = FALSE]
  cols <- setdiff(intersect(colnames(old), colnames(new)), key)
  Filter(function(cn) {
    a <- as.character(old[[cn]]); a[is.na(a)] <- ""
    b <- as.character(new[[cn]]); b[is.na(b)] <- ""
    any(a != b)
  }, cols)
}



#### Execute the editor (repeats until the checks pass) ########################
cat("\n[INFO] A table editor will open in your web browser.\n",
    "       Edit both tabs, then click 'Done' at the top right.\n", sep = "")

current_md <- metadata_in
current_pn <- panel_in
repeat {
  edited <- run_editor(current_md, current_pn)
  if (is.null(edited)) {
    stop("Editing was cancelled. No files were changed.", call. = FALSE)
  }
  current_md <- edited$metadata
  current_pn <- edited$panel

  problems <- find_edit_problems(current_md, current_pn)
  if (length(problems$errors) > 0) {
    svDialogs::dlg_message(paste0(
      "These problems must be fixed before the files can be saved:\n\n- ",
      paste(problems$errors, collapse = "\n- "),
      if (length(problems$warnings) > 0) paste0(
        "\n\nAlso, for your information:\n- ",
        paste(problems$warnings, collapse = "\n- ")) else "",
      "\n\nThe editor will reopen with your changes."), type = "ok")
    next
  }
  if (length(problems$warnings) > 0) {
    if (!ask_yes_no(paste0(
      "Please check:\n\n- ", paste(problems$warnings, collapse = "\n- "),
      "\n\nGo back and fix these?\nYes = reopen the editor\nNo = save as it is"))) {
      break
    }
    next
  }
  break
}



#### Put rows back in their original order and fix column types ################
# (the editor lets you sort; the saved files must not depend on that)
if (!setequal(current_md$File_Name, metadata_in$File_Name) ||
    !setequal(current_pn$Name, panel_in$Name)) {
  stop("File_Name / channel Name values were changed or lost in the editor. ",
       "Nothing was saved.", call. = FALSE)
}
final_md <- current_md[order(current_md$File_Num), , drop = FALSE]
final_pn <- current_pn[match(panel_in$Name, current_pn$Name), , drop = FALSE]
final_md$File_Num <- as.integer(final_md$File_Num)
final_pn$Analysis <- as.integer(final_pn$Analysis)
final_pn$Plot     <- as.integer(final_pn$Plot)
final_pn$Factor   <- as.numeric(final_pn$Factor)
rownames(final_md) <- NULL
rownames(final_pn) <- NULL

edited_md_cols <- changed_columns(metadata_in, final_md, "File_Num")
edited_pn_cols <- changed_columns(panel_in,    final_pn, "Name")



#### Save ######################################################################
metadata_path <- file.path(Original_WD, METADATA_NAME)
panel_path    <- file.path(Original_WD, PANEL_NAME)
backup_file(metadata_path)
backup_file(panel_path)
fwrite(final_md, file = metadata_path)
fwrite(final_pn, file = panel_path)

# Keep the updated tables in memory for the next scripts
metadata <- final_md
panel    <- final_pn

cat("\nAnalysis channels (clustering): ", sum(final_pn$Analysis == 1),
    "\nPlot channels:                  ", sum(final_pn$Plot == 1),
    "\nSampleType: ", paste(names(table(final_md$SampleType)), table(final_md$SampleType),
                            sep = " = ", collapse = ", "), "\n", sep = "")



#### Record parameters used in this run ########################################
Inputs <- list(
  # --- keys that existed in v4 (names unchanged) ---
  Original_WD       = Original_WD,
  METADATA_FILENAME = METADATA_NAME,
  PANEL_FILENAME    = PANEL_NAME,
  # --- new in v5 ---
  Script_Version           = SCRIPT_VERSION,
  Reviewed                 = TRUE,
  Reviewed_At              = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  Metadata_Columns_Edited  = as.list(edited_md_cols),
  Panel_Columns_Edited     = as.list(edited_pn_cols),
  Analysis_Channels        = as.list(final_pn$Fixed_Param[final_pn$Analysis == 1]),
  Plot_Channels            = as.list(final_pn$Fixed_Param[final_pn$Plot == 1]),
  N_Analysis_Channels      = sum(final_pn$Analysis == 1),
  N_Plot_Channels          = sum(final_pn$Plot == 1),
  SampleType_Counts        = as.list(table(final_md$SampleType)),
  Unresolved_Warnings      = as.list(find_edit_problems(final_md, final_pn)$warnings),
  Random_Seed              = "not used in this script",
  Output_MD5               = file_md5(c(metadata_path, panel_path)),
  Session                  = session_info_list(names(pkg_versions))
)
record_params(MODULE_NAME, Inputs, write_dir = Original_WD)



#### Script completed message ##################################################
print_banner(" Cleanup complete ",
             "'metadata.csv' and 'panel.csv' successfully updated. If you will be performing batch correction, continue with 01.3_Update_ChannelsForBatchCorr; otherwise continue with 02_FCS_to_CSV.")
