################################################################################
# 01.3_Update_ChannelsForBatchCorr_v5.R
################################################################################
# Layer:    DataCleanup -- OPTIONAL step, only if you will do batch correction
#           (run after 01.2_Update_Metadata, before BatchCorrection/02_RunBatch)
# Author:   Jonathon Sewell
# Created:  2026-07-22
# Version:  5   (input version: 01.3_Update_ChannelsForBatchCorr_v3; v4 was an interim test version)
# Edited:   2026-10-21
#
# Purpose:
# * Confirm which channels will be adjusted by batch correction and save the
#   list as ChannelsToAdjust.txt (one channel name per line) in your working
#   folder. The starting list is every channel marked Plot = 1 in panel.csv
#   (set in 01.2). ChannelsToAdjust_Template.txt, if present, is only used to
#   tell you where your list differs from it.
# * You can accept the list as shown, or review it in a table and mark each
#   channel Y or N.
#
# Origin: N/A
# * Previously required the user to edit the files in Excel or TextEdit
#
# ==============================================================================
# CHANGE LOG: v3 -> v5
# ==============================================================================
# Convention:  [RESULT]   = can change which channels end up in ChannelsToAdjust.txt
#              [BEHAVIOR] = changes how the script behaves, not the data
#              [SAME]     = output unchanged
#
# v4 -> v5 (after your test run) -- you answered: "adjust Plot = 1 channels"
#  C1 [RESULT] The default list is now the channels with Plot = 1 in panel.csv
#     (dataset 2: 40 channels). The template no longer decides the list. Before,
#     Yes wrote the whole template (42 channels). If Plot = 1 is not set for any
#     channel, the template (if present) or all metal channels are used.
#  C2 [BEHAVIOR] The pop-up now also shows: Plot = 1 channels missing from the
#     template, and template channels left out because Plot = 0.
#  C3 [BEHAVIOR] The editor lists ALL metal channels in the panel (so you can add
#     a channel), Include starts as Y for the default list. Fixed_Param is shown
#     (read-only) as an aid.
#  Items A1-A5 and B1-B6 below describe v3 -> v4 and still apply, except A1's
#  question, which is answered by C1.
#
# A. Needs your attention
#  A1 [SAME output, QUESTION] The "Is this correct?" pop-up in v3 listed each
#     channel as "Name = TRUE/FALSE" under the heading "Channels that will be
#     adjusted". TRUE meant "Plot = 1 in panel.csv". But clicking Yes wrote ALL
#     listed channels to ChannelsToAdjust.txt, TRUE or FALSE (your current
#     ChannelsToAdjust.txt is identical to the template, i.e. all 42 channels).
#     v4 keeps exactly that behaviour but the pop-up now says what is true:
#     the whole list is adjusted, and channels are grouped by Plot = 1 / Plot = 0
#     for your information. Clicking No still lets you mark each channel Y/N,
#     where Y is pre-set for Plot = 1 channels, as in v3.
#     -> Is "adjust everything in the template" what you intend, or should the
#        default be "only Plot = 1 channels"? Say so and I will change the default.
#  A2 [BEHAVIOR] The editor is opened with browserViewer(), as in 01.2. v3 called
#     browserViewer("Data Cleanup Module", width = , height = ), which I believe
#     is not a valid use of that function (not tested -- the "No" path had
#     probably not been used).
#  A3 [BEHAVIOR] The Yes/No answer to "batch correction?" is read from
#     run_parameters.json (written by 01_Get_Metadata) instead of searching the
#     folder for any file containing "json". If you answered No in 01, the
#     script now asks whether you want to continue anyway (v3 stopped with
#     "object 'ChannelsToAdjust_in' not found").
#  A4 [BEHAVIOR] Channels in the list that are NOT in panel.csv are reported.
#     In v3 they were dropped silently if you clicked No (editor) and written to
#     the file if you clicked Yes. v4 does the same as v3 in both cases but
#     tells you.
#  A5 [BEHAVIOR] A time-stamped copy of an existing ChannelsToAdjust.txt is saved
#     in _backups before it is replaced.
#
# B. Efficiency / readability / reproducibility
#  B1 [SAME] File written is identical in format (one channel per line, no header,
#     no quotes).
#  B2 Shared helper file (00_Helpers.R) loaded by a short block at the top.
#  B3 Sequence check: stops with a plain message unless 01_Get_Metadata and
#     01.2_Update_Metadata have been run in this folder.
#  B4 panel.csv is read from disk (so edits made outside R are used).
#  B5 dplyr is no longer needed (rename/arrange/select/mutate/filter replaced by
#     base R). The memory object 'ChannelsToAdjust', which nothing created, is
#     gone. No automatic install.packages().
#  B6 Run record: lists the channels saved, where the starting list came from,
#     whether you accepted or edited it, and a checksum of the file. v3 keys keep
#     their names (WRITE_PATH removed with the project-wide change).
################################################################################


#### USER SETTINGS -- edit this section only ###################################
# (none: everything is chosen in the pop-ups)

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "01.3_Update_ChannelsForBatchCorr_v5.R"
SCRIPT_VERSION           <- 5
MODULE_NAME              <- "01.3_Update_ChannelsForBatchCorr"
PREREQUISITE_MODULES     <- c("01_Get_Metadata", "01.2_Update_Metadata")
REQUIRED_HELPERS_VERSION <- 2
CHANNELS_TEMPLATE_FILENAME <- "ChannelsToAdjust_Template.txt"
CHANNELS_FILENAME          <- "ChannelsToAdjust.txt"



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
cat("[START] Modify files necessary for batch correction:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

Original_WD <- get_working_dir()

# Sequence check: 01 and 01.2 must have been run in this folder
require_modules(Original_WD, PREREQUISITE_MODULES,
                how_to_fix = "Run 01_Get_Metadata and 01.2_Update_Metadata first.")



#### Was batch correction selected in 01? ######################################
run_record  <- read_run_record(Original_WD)
batch_flag  <- get_record_value(run_record, "01_Get_Metadata", "batchCorrection_res",
                                default = NA_character_)
if (!identical(batch_flag, "yes")) {
  if (!ask_yes_no(paste0(
    "In 01_Get_Metadata you answered that you will NOT perform batch correction.\n\n",
    "Set up batch correction anyway?"))) {
    stop("Nothing to do: batch correction was not selected. ",
         "Continue with 02_FCS_to_CSV.", call. = FALSE)
  }
}



#### Read panel.csv and the starting list of channels ##########################
panel <- read_panel(Original_WD)

read_channel_list <- function(path) {
  x <- trimws(readLines(path, warn = FALSE))
  x[nzchar(x)]
}

template_path <- file.path(Original_WD, CHANNELS_TEMPLATE_FILENAME)
template_list <- if (file.exists(template_path)) unique(read_channel_list(template_path)) else character(0)

metal_channels <- panel$Name[grepl("[1-9]", panel$Name)]
plot1_channels <- panel$Name[panel$Plot == 1]

if (length(plot1_channels) > 1) {
  start_source <- "panel.csv (Plot = 1)"
  start_list   <- plot1_channels
} else if (length(template_list) > 0) {
  start_source <- "template"
  start_list   <- template_list
  cat("\n[INFO] Fewer than 2 channels have Plot = 1; using", CHANNELS_TEMPLATE_FILENAME, "\n")
} else {
  start_source <- "panel.csv (all metal channels)"
  start_list   <- metal_channels
  cat("\n[INFO] Fewer than 2 channels have Plot = 1 and no template; using all metal channels.\n")
}
start_list <- unique(start_list)
if (length(start_list) == 0) {
  stop("The starting list of channels is empty. Check panel.csv.", call. = FALSE)
}
cat("\n[INFO] Starting list:", length(start_list), "channels from", start_source, "\n")
plot1_not_in_template <- if (length(template_list) > 0) setdiff(start_list, template_list) else character(0)
template_not_used     <- if (length(template_list) > 0) setdiff(template_list, start_list) else character(0)



#### Functions specific to this script #########################################
# Shiny gadget: one editable table (Name, Plot, Include).
# Returns the edited table or NULL if cancelled.
run_channel_editor <- function(df) {
  ui <- miniPage(
    gadgetTitleBar("Data Cleanup: Confirm Channels to Adjust with Batch Correction"),
    miniTabstripPanel(
      miniTabPanel("Channels to Adjust", icon = icon("table"),
                   miniContentPanel(
                     helpText("Set Include to Y for every channel that should be batch-adjusted, N for the rest. ",
                              "Name, Fixed_Param and Plot cannot be changed (Plot = 1 means the channel is plotted, from 01.2)."),
                     rHandsontableOutput("channels_grid", height = "100%")
                   )
      )
    )
  )

  server <- function(input, output, session) {
    output$channels_grid <- renderRHandsontable({
      d   <- if (!is.null(input$channels_grid)) hot_to_r(input$channels_grid) else df
      hot <- rhandsontable(d,
                           rowHeaders = FALSE,
                           height = 550,
                           overflow = "visible",
                           columnSorting = TRUE,
                           manualColumnResize = TRUE)
      hot <- hot_col(hot, c("Name", "Fixed_Param", "Plot"), readOnly = TRUE)
      hot <- hot_col(hot, "Include", type = "dropdown", source = c("Y", "N"), strict = TRUE)
      hot_table(hot, highlightCol = TRUE, highlightRow = TRUE, fontSize = 11)
    })

    observeEvent(input$done, {
      final <- if (!is.null(input$channels_grid)) hot_to_r(input$channels_grid) else df
      stopApp(list(channels = final))
    })
  }

  res <- tryCatch(runGadget(ui, server, viewer = browserViewer()),
                  error = function(e) {
                    if (grepl("cancel", conditionMessage(e), ignore.case = TRUE)) NULL else stop(e)
                  })
  if (is.null(res)) NULL else res$channels
}



#### Ask the user to confirm the list ##########################################
in_panel  <- start_list %in% panel$Name
plot_flag <- start_list %in% panel$Name[panel$Plot == 1]

list_summary <- paste0(
  "Channels that will be adjusted by batch correction (", length(start_list), "):\n\n",
  "Plot = 1 in panel (", sum(plot_flag), "): ",
  if (any(plot_flag)) paste(start_list[plot_flag], collapse = ", ") else "none", "\n\n",
  "Plot = 0 in panel (", sum(in_panel & !plot_flag), "): ",
  if (any(in_panel & !plot_flag)) paste(start_list[in_panel & !plot_flag], collapse = ", ") else "none", "\n\n",
  if (any(!in_panel)) paste0("NOT FOUND in panel.csv (", sum(!in_panel), "): ",
                             paste(start_list[!in_panel], collapse = ", "), "\n\n") else "",
  if (length(template_list) > 0) paste0(
    "Compared with ", CHANNELS_TEMPLATE_FILENAME, ":\n",
    "  in your list but not in the template (", length(plot1_not_in_template), "): ",
    if (length(plot1_not_in_template)) paste(plot1_not_in_template, collapse = ", ") else "none", "\n",
    "  in the template but not in your list (", length(template_not_used), "): ",
    if (length(template_not_used)) paste(template_not_used, collapse = ", ") else "none", "\n\n") else "",
  "Is this correct?\n",
  "Yes = adjust ALL channels listed above\n",
  "No  = review the list and choose channels in an editor")

accepted_as_listed <- ask_yes_no(list_summary)

if (accepted_as_listed) {
  channels_final <- start_list
  if (any(!in_panel)) {
    warning("These channels are not in panel.csv but are included in ",
            CHANNELS_FILENAME, ": ", paste(start_list[!in_panel], collapse = ", "),
            immediate. = TRUE, call. = FALSE)
  }
} else {
  # Editor: ALL metal channels in the panel; Include = Y for the starting list
  editor_df <- panel[panel$Name %in% metal_channels | panel$Name %in% start_list,
                     c("Name", "Fixed_Param", "Plot"), drop = FALSE]
  if (nrow(editor_df) == 0) {
    stop("No metal channels were found in panel.csv.", call. = FALSE)
  }
  editor_df$Include <- ifelse(editor_df$Name %in% start_list, "Y", "N")
  editor_df$Plot    <- as.character(editor_df$Plot)
  rownames(editor_df) <- NULL
  if (any(!in_panel)) {
    cat("[INFO] Not shown in the editor (not in panel.csv): ",
        paste(start_list[!in_panel], collapse = ", "), "\n", sep = "")
  }

  cat("\n[INFO] A table editor will open in your web browser.\n",
      "       Set Include (Y/N) for each channel, then click 'Done'.\n", sep = "")
  edited <- run_channel_editor(editor_df)
  if (is.null(edited)) {
    stop("Editing was cancelled. No files were changed.", call. = FALSE)
  }
  edited <- edited[match(editor_df$Name, edited$Name), , drop = FALSE]  # panel order
  channels_final <- edited$Name[edited$Include == "Y"]
}

if (length(channels_final) == 0) {
  stop("No channels were selected, so ", CHANNELS_FILENAME, " was not written.",
       call. = FALSE)
}



#### Save ######################################################################
channels_path <- file.path(Original_WD, CHANNELS_FILENAME)
backup_file(channels_path)
writeLines(channels_final, channels_path)   # one channel per line, no header, no quotes

cat("\n", length(channels_final), " channels written to ", CHANNELS_FILENAME, "\n", sep = "")



#### Record parameters used in this run ########################################
Inputs <- list(
  # --- keys that existed in v3 (names unchanged) ---
  Original_WD                = Original_WD,
  PANEL_FILENAME             = PANEL_NAME,
  CHANNELS_TEMPLATE_FILENAME = CHANNELS_TEMPLATE_FILENAME,
  CHANNELS_FILENAME          = CHANNELS_FILENAME,
  # --- new in v4/v5 ---
  Script_Version        = SCRIPT_VERSION,
  Completed_At          = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  Batch_Selected_In_01  = identical(batch_flag, "yes"),
  Starting_List_Source  = start_source,
  Template_Found        = length(template_list) > 0,
  Choice                = if (accepted_as_listed) "accepted list as shown" else "edited in editor",
  N_Channels            = length(channels_final),
  Channels_To_Adjust    = as.list(channels_final),
  Channels_Not_In_Panel = as.list(setdiff(channels_final, panel$Name)),
  Random_Seed           = "not used in this script",
  Output_MD5            = file_md5(channels_path),
  Session               = session_info_list(names(pkg_versions))
)
record_params(MODULE_NAME, Inputs, write_dir = Original_WD)



#### Script completed message ##################################################
print_banner(" Cleanup complete ",
             "'ChannelsToAdjust.txt' successfully updated. Continue with BatchCorrection/02_RunBatch.")
