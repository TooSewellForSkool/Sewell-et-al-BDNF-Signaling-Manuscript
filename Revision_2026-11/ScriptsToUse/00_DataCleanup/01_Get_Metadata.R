################################################################################
# 01_Get_Metadata_v8.R
################################################################################
# Layer:    DataCleanup -- FIRST script in the pipeline
# Author:   Jonathon Sewell
# Created:  2026-06-08
# Version:  8   (v5 carried "Version: 4" inside the file; numbering follows the
#                file names from here on)
# Edited:   2026-10-21
#
# Purpose:
# * Read the header of every FCS file in your working folder (no cell data is
#   loaded) and write two files for you to review/edit in the next script
#   (01.2_Update_Metadata):
#     - metadata.csv : one row per FCS file (parsed from the file name + event count)
#     - panel.csv    : one row per channel in the antibody panel
# * Record everything needed to repeat this run in run_parameters.json
#
# Origin:
# * Output metadata
#     > 01_Get_Metadata.R - smgoggin, 05.14.19
# * Determine events per FCS file and append to metadata.csv
#     > Cells_Per_Cluster.R - Corey Williams, University of Virginia, May 17, 2019
#     > 01.1_Append_Metadata_Columns.R - Jonathon Sewell, 2025-06-17
#       (implements logic from Cells_Per_Cluster.R)
#
# ==============================================================================
# CHANGE LOG: v5 -> v8
# ==============================================================================
# Convention:  [RESULT]   = can change the content of metadata.csv / panel.csv
#              [SAME]     = output unchanged (tested where noted)
#              [BEHAVIOR] = changes how the script behaves, not the data
#
# --- Changes you approved --------------------------------------------------
#  A1 [RESULT] File-name parser: Condition is now everything between "normalized"
#     and the time token (e.g. Cont._Depr., BDNF_Test, BDNF+K252a_Test), and
#     Treatment is Condition + time token (e.g. Cont._Depr._15min). A time token
#     is a WHOLE token like 8hr, 15min or 1.5hr. v5 used only the first token as
#     Condition, searched tokens loosely for "hr"/"min", and dropped decimals
#     (1.5hr -> 15). Files with no time token (Universal_1, HEK_Test_1,
#     DIV4_Stock_1_3) are parsed exactly as in v5.
#     Tested by re-implementing the rules in Python on both of your metadata
#     files: dataset 1 -> 0 of 40 files differ from v5 output; dataset 2 -> only
#     the 24 files with multi-token conditions change (Cont._Depr.* and *_Test*),
#     including their Replicate numbers (now counted per Treatment+time).
#  A2 [SAME]   File_Num is an integer in memory (v5 kept it as text). The file on
#     disk looks the same.
#  A3 [BEHAVIOR] FCS files are no longer all loaded into memory. Only the file
#     headers are read (event count + channel names for every file); the first
#     file is also read in full to cross-check the header numbers and to build
#     the panel. The 'fs' object is no longer created -- 02_FCS_to_CSV will read
#     the files itself.
#  A4 [BEHAVIOR] Removed WRITE_PATH and command-line arguments. Outputs always go
#     in your working folder (the folder with the FCS files).
#  A5 [BEHAVIOR] Removed INTENDED_CLUSTER_ROUND and filenums.csv handling.
#     Choosing which cells go into a later clustering round moves to the
#     rebuilt 02.1_SubsetAndExport.
#  A6 [BEHAVIOR] Backups now go into a "_backups" subfolder (fixes 01.3, which
#     searched the folder for any file containing "json").
#  A7 [BEHAVIOR] Shared helper file (00_Helpers.R) is loaded by a short block at
#     the top. The script still runs on its own as long as ScriptsToUse/Helpers
#     is in place.
#  A8 [BEHAVIOR] Run record: writes run_parameters.json as before, now also
#     merging with the file on disk, full numeric precision, and extra fields
#     (see "Record parameters"). All v5 keys keep their names, except
#     WRITE_PATH and INTENDED_CLUSTER_ROUND, which were removed (A4, A5).
#  A9 [BEHAVIOR] The batch-correction template is located by searching the
#     DataCleanup layer folder, so it no longer depends on a "Modules" folder.
#
# --- v7 -> v8 (bug fix after your test run) -----------------------------------
#  F1 [BEHAVIOR] v7 warned "Could not copy ChannelsToAdjust_Template.txt ...
#     found 0". The template search only looked inside a folder whose name
#     contains "DataCleanup" next to Helpers. Helpers v2 now also searches
#     everything under the folder that holds Helpers (ScriptsToUse), ignoring
#     _backups and hidden folders, and 01 uses the first match instead of
#     demanding exactly one. If no template exists anywhere it is now an INFO
#     message, not a warning: 01.3 builds its list from panel.csv.
#     Data in metadata.csv / panel.csv is unaffected.
#
# --- Efficiency / readability ----------------------------------------------
#  E1 [SAME] One USER SETTINGS section at the top; fixed names listed separately.
#  E2 [BEHAVIOR] Package check with install instructions (flowCore is Bioconductor).
#  E3 [SAME] File-name parsing, preview and checks rewritten as functions that
#     split each name once (v5 re-split it ~20 times per file).
#  E4 [BEHAVIOR] File names are checked against the naming convention BEFORE any
#     FCS data is read, and the preview now lists every Condition found (v5
#     showed only the first file).
#  E5 [SAME] Parameter-name cleanup vectorized (v5 looped over channels).
#  E6 [SAME] Replicate numbering with base R ave() (v5: dplyr). Event counts with
#     vapply() instead of rbind() in a loop. dplyr is no longer needed here.
#  E7 [SAME] metadata.csv rows stay in File_Num order. v5 used merge(), which
#     re-sorted rows by File_Name; the two orders match in practice (golden-file
#     check will confirm on your data).
#  E8 [BEHAVIOR] Existing metadata.csv / panel.csv / run_parameters.json are NOT
#     overwritten without asking; a timestamped backup is saved first.
#  E9 [BEHAVIOR] Warning if channel names/markers differ between FCS files.
# E10 [BEHAVIOR] Plain-language errors for: missing packages, no FCS files,
#     closed pop-up windows, unrecognised file names, unsaved script.
# E11 [SAME] panel.csv: Analysis and Plot are integers, Factor is numeric in
#     memory (v5: all text). Text on disk is identical.
# E12 [SAME] Cosmetic: banner function, file.path(), removed "<<-" at top level.
#
# --- Left as in v5 on purpose (flagged to you, not applied) ----------------
#  * list.files(pattern = "*.fcs") and its order (File_Num follows the order your
#    computer's locale sorts the file names; the locale is recorded in the JSON).
#  * Panel is taken from the FIRST file; a warning is shown if others differ.
#  * 5 min -> 0.083 h (rounded to TIMEPOINT_DIGITS decimals, as in v5).
################################################################################


#### USER SETTINGS -- edit this section only ###################################
# Number of decimal places kept when minutes are converted to hours
# (5 min = 0.083 h). Same as v5; change only if you need more precision.
TIMEPOINT_DIGITS <- 3

#### Fixed names -- do not change (later scripts expect these) ##################
SCRIPT_NAME              <- "01_Get_Metadata_v8.R"
SCRIPT_VERSION           <- 8
MODULE_NAME              <- "01_Get_Metadata"
REQUIRED_HELPERS_VERSION <- 2
CHANNELS_TEMPLATE_NAME   <- "ChannelsToAdjust_Template"
# A whole-token time label: 8hr, 16hr, 15min, 1.5hr ...
TIME_TOKEN_REGEX         <- "^[0-9]+(\\.[0-9]+)?(hr|min)$"



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
  cran = c("svDialogs", "rstudioapi", "data.table", "jsonlite"),
  bioc = c("flowCore"))

cat("\n========================================================\n")
cat("[START] Initialize data cleanup pipeline:\n",
    strrep(" ", nchar("[START] ")), "> Start ", SCRIPT_NAME, "\n", sep = "")
cat("========================================================\n")

# Working folder = the folder that holds the FCS files. Outputs go here too.
Original_WD <- get_working_dir()

if (!is.numeric(TIMEPOINT_DIGITS) || length(TIMEPOINT_DIGITS) != 1 ||
    TIMEPOINT_DIGITS < 0 || TIMEPOINT_DIGITS != round(TIMEPOINT_DIGITS)) {
  stop("TIMEPOINT_DIGITS must be a whole number of 0 or more.", call. = FALSE)
}



#### Functions specific to this script #########################################
# Split a file name into its '_'-separated parts and locate the key parts.
split_file_name <- function(file_name) {
  base     <- sub("\\.fcs$", "", file_name, ignore.case = TRUE)
  parts    <- strsplit(base, "_", fixed = TRUE)[[1]]
  norm_idx <- grep("norm", parts, fixed = TRUE)           # e.g. "normalized"
  after    <- if (length(norm_idx) == 1) {
    seq_along(parts)[seq_along(parts) > norm_idx]         # parts after 'normalized'
  } else {
    integer(0)
  }
  time_idx <- after[grepl(TIME_TOKEN_REGEX, parts[after])] # whole-token time label
  list(parts = parts, norm_idx = norm_idx, after = after, time_idx = time_idx)
}

# "" if the file name follows the naming convention, otherwise what is wrong
find_naming_problem <- function(file_name) {
  if (!grepl("\\.fcs$", file_name)) {
    return("does not end in '.fcs'")
  }
  s <- split_file_name(file_name)
  if (length(s$norm_idx) != 1) {
    return("needs exactly one '_'-separated part containing 'norm' (e.g. 'normalized')")
  }
  if (length(s$after) == 0) {
    return("has nothing after the 'normalized' part")
  }
  if (length(s$time_idx) > 1) {
    return("has more than one time label (like 8hr or 15min) after 'normalized'")
  }
  if (length(s$time_idx) == 1 && s$time_idx == s$norm_idx + 1) {
    return("has no treatment name between 'normalized' and the time label")
  }
  ""
}

# Build the metadata columns (all except File_Num and Event_Count) for one file.
# Placeholders ("User-Defined ...") are meant to be edited in 01.2.
parse_file_name <- function(file_name, digits = 3) {
  s        <- split_file_name(file_name)
  parts    <- s$parts
  has_univ <- any(grepl("Univ", parts, fixed = TRUE))   # sets Condition/Collection
  has_niv  <- any(grepl("niv",  parts, fixed = TRUE))   # sets TimePoint (as in v5)

  if (length(s$time_idx) == 1) {
    # Condition = everything between 'normalized' and the time label
    condition <- paste(parts[(s$norm_idx + 1):(s$time_idx - 1)], collapse = "_")
    treatment <- paste(condition, parts[s$time_idx], sep = "_")
    timepoint <- as.numeric(sub("(hr|min)$", "", parts[s$time_idx]))
    if (grepl("min$", parts[s$time_idx])) {
      timepoint <- round(timepoint / 60, digits = digits)   # minutes -> hours
    }
  } else {
    # No time label: first part = Condition, first two parts = Treatment (as in v5)
    condition <- parts[s$norm_idx + 1]
    treatment <- if (s$norm_idx + 2 <= length(parts)) {
      paste(parts[s$norm_idx + 1], parts[s$norm_idx + 2], sep = "_")
    } else {
      parts[s$norm_idx + 1]
    }
    timepoint <- suppressWarnings(as.numeric(parts[s$norm_idx + 2]))  # NA unless numeric
  }
  if (has_niv)  timepoint <- NA_real_      # universal samples have no time point
  if (has_univ) condition <- "Univ"

  data.frame(File_Name        = file_name,
             SampleCollection = if (has_univ) "Univ" else "User-Defined",
             Expt_ID          = "User-Defined #",
             BarcodeSet       = parts[1],
             Treatment        = treatment,
             Condition        = condition,
             TimePoint        = timepoint,
             Replicate        = NA_integer_,         # filled in below
             SampleType       = "Expt/Ctrl",
             stringsAsFactors = FALSE)
}

# Text for the confirmation pop-up: how the first file name was read, plus every
# Condition found
describe_parsing <- function(first_name, parsed) {
  s <- split_file_name(first_name)
  p <- s$parts
  if (length(s$time_idx) == 1) {
    time_txt <- p[s$time_idx]
    rep_txt  <- if (s$time_idx + 1 <= length(p)) p[s$time_idx + 1] else "(none)"
    rest_txt <- if (s$time_idx + 2 <= length(p)) {
      paste(p[(s$time_idx + 2):length(p)], collapse = "_")
    } else "(none)"
  } else {
    time_txt <- "(no time label found)"
    rep_txt  <- "(not detected)"
    rest_txt <- if (s$norm_idx + 2 <= length(p)) {
      paste(p[(s$norm_idx + 2):length(p)], collapse = "_")
    } else "(none)"
  }
  first_row <- parsed[1, ]
  conds <- unique(parsed$Condition)
  cond_lines <- vapply(conds, function(cn) {
    rows <- parsed$Condition == cn
    tp   <- sort(unique(parsed$TimePoint[rows]))
    sprintf("  %s: %d file(s); time points (h): %s", cn, sum(rows),
            if (length(tp) == 0) "none" else paste(tp, collapse = ", "))
  }, character(1))
  paste0("FIRST FILE: ", first_name, "\n\n",
         "BarcodeID: ", p[1], "\n",
         "Normalized: ", p[s$norm_idx], "\n",
         "Condition: ", first_row$Condition, "\n",
         "TimePoint label: ", time_txt, "\n",
         "Rep: ", rep_txt, "\n",
         "Remaining: ", rest_txt, "\n\n",
         "ALL CONDITIONS FOUND (", length(conds), "):\n",
         paste(cond_lines, collapse = "\n"))
}

# Clean up channel names for the panel: use Desc (marker) if available, else Name;
# move a leading mass tag to the end ("89Y_Tuj1" -> "Tuj1_89Y"); spaces and "/"
# become "_". (Hyphens are NOT changed, e.g. "p-p38".)
fix_param_names <- function(name, desc) {
  fixed <- ifelse(is.na(desc), name, desc)
  has_mass_prefix <- grepl("^\\d+\\w+_", fixed)
  fixed[has_mass_prefix] <- paste0(sub("^[^_]+_", "", fixed[has_mass_prefix]), "_",
                                   sub("_.*$", "", fixed[has_mass_prefix]))
  gsub("[ /]", "_", fixed)
}

# Read only the header of one FCS file. Returns a named character vector of
# keywords ("$TOT", "$PAR", "$P1N", ...) or NULL if that fails.
read_fcs_header <- function(file_path) {
  hdr <- tryCatch(
    flowCore::read.FCSheader(basename(file_path), path = dirname(file_path)),
    error = function(e) NULL)
  if (is.null(hdr) || length(hdr) == 0) return(NULL)
  hdr <- hdr[[1]]
  if (is.list(hdr)) hdr <- hdr[[1]]          # files holding several datasets
  hdr <- unlist(hdr)
  names(hdr) <- toupper(names(hdr))
  hdr
}



#### Find and check the FCS files ##############################################
file_list <- list.files(Original_WD, full.names = TRUE, pattern = "*.fcs")
if (length(file_list) == 0) {
  stop("No .fcs files were found in:\n  ", Original_WD, "\n",
       "Your working folder must be the folder that contains your FCS files.\n",
       "Start a new R session (Session > Restart R) or click No at the folder ",
       "question and pick the right folder.", call. = FALSE)
}
file_names <- basename(file_list)

# Check every name against the naming convention (before reading any data)
problems <- vapply(file_names, find_naming_problem, character(1), USE.NAMES = FALSE)
if (any(problems != "")) {
  bad   <- which(problems != "")
  shown <- head(bad, 10)
  stop("Some FCS file names do not follow the expected naming convention:\n",
       paste0("  - ", file_names[shown], ": ", problems[shown], collapse = "\n"),
       if (length(bad) > 10) paste0("\n  ... and ", length(bad) - 10, " more") else "",
       "\n\nExpected: [BarcodeID]_[Indication of Normalization]_[Treatment]_",
       "[TimePoint]_[Rep]_[Indication of Gating]_[Add'l Gating Info].fcs",
       call. = FALSE)
}

# Parse all names, then ask the user to confirm
parsed <- do.call(rbind, lapply(file_names, parse_file_name, digits = TIMEPOINT_DIGITS))
parsed$Replicate <- ave(seq_len(nrow(parsed)), parsed$Treatment, FUN = seq_along)

cat("[INFO] Note: you may need to adjust columns to fit the naming nomenclature",
    "of your dataset\n")
if (!ask_yes_no(paste0("Identified metadata components from input naming convention:\n\n",
                       describe_parsing(file_names[1], parsed),
                       "\n\nAre these correct?"))) {
  stop("\n=======================================================================\n",
       "Revisit the File_Name convention.\n",
       "For the current scripts, the file naming convention is:\n",
       "[BarcodeID]_[Indication of Normalization]_[Treatment]_[TimePoint]_[Rep]_",
       "[Indication of Gating]_[Add'l Gating Info].fcs\n",
       "=======================================================================",
       call. = FALSE)
}

# Protect files from earlier runs (they may contain your edits from 01.2)
metadata_path <- file.path(Original_WD, METADATA_NAME)
panel_path    <- file.path(Original_WD, PANEL_NAME)
record_path   <- file.path(Original_WD, RUN_RECORD_NAME)
existing      <- c(metadata_path, panel_path, record_path)
existing      <- existing[file.exists(existing)]
if (length(existing) > 0) {
  old_record <- read_run_record(Original_WD)
  later_steps <- setdiff(names(old_record), c("Folder_Info", MODULE_NAME))
  if (!ask_yes_no(paste0(
    "These files already exist in this folder and may contain your edits:\n\n",
    paste(basename(existing), collapse = "\n"),
    if (length(later_steps) > 0) paste0(
      "\n\nThe run record shows these later steps were already run:\n",
      paste(later_steps, collapse = "\n")) else "",
    "\n\nOverwrite them? (A timestamped copy of each is saved in _backups first.)"))) {
    stop("Stopped so that your existing files are not overwritten. ",
         "Nothing was changed.", call. = FALSE)
  }
}



#### Read FCS headers ##########################################################
cat("[INFO] Reading the headers of", length(file_list), "FCS files\n")
headers <- lapply(file_list, read_fcs_header)

# Read the FIRST file in full: it defines the panel, and it checks the header
# numbers (a safety net for how the headers are read)
cat("[INFO] Reading", file_names[1], "in full to define the panel\n")
first_fcs <- flowCore::read.FCS(file_list[1], transformation = FALSE)

event_counts <- vapply(headers, function(h) {
  if (is.null(h) || is.na(h["$TOT"])) NA_integer_ else
    suppressWarnings(as.integer(h["$TOT"]))
}, integer(1))

if (any(is.na(event_counts)) || event_counts[1] != nrow(first_fcs)) {
  warning("FCS header event counts could not be confirmed (first file: header ",
          event_counts[1], " vs data ", nrow(first_fcs), "). ",
          "Counting events by reading every file in full -- this may take a while.",
          immediate. = TRUE, call. = FALSE)
  event_counts <- vapply(file_list, function(f)
    nrow(flowCore::read.FCS(f, transformation = FALSE)), integer(1), USE.NAMES = FALSE)
}

# Panel consistency: compare each file's channel names/markers with the first file
param_info <- flowCore::pData(flowCore::parameters(first_fcs))
panel_name <- as.character(unname(param_info$name))
panel_desc <- as.character(unname(param_info$desc))
n_par      <- length(panel_name)

header_signature <- function(h) {
  if (is.null(h)) return(NULL)
  unname(h[c(paste0("$P", seq_len(n_par), "N"), paste0("$P", seq_len(n_par), "S"))])
}
signatures       <- lapply(headers, header_signature)
panel_consistent <- !any(vapply(signatures, is.null, logical(1))) &&
  all(vapply(signatures, identical, logical(1), signatures[[1]]))
if (!panel_consistent) {
  warning("Not every FCS file has the same channel names / markers as the first ",
          "file (or a header could not be read). panel.csv is based on the FIRST ",
          "file (", file_names[1], "). Check that this is what you want.",
          immediate. = TRUE, call. = FALSE)
}



#### Build panel.csv and metadata.csv ##########################################
is_mass_channel <- grepl("Di$|bc_neg|bc_pos", panel_name)

panel_df <- data.frame(Name        = panel_name,
                       Desc        = panel_desc,
                       Fixed_Param = fix_param_names(panel_name, panel_desc),
                       Analysis    = 0L,
                       Plot        = 0L,
                       Transform   = ifelse(is_mass_channel, "Asinh", "Linear"),
                       Factor      = ifelse(is_mass_channel, 5, 1),
                       stringsAsFactors = FALSE)

# File_Num = position of the file in the list of files (this is the number later
# scripts use to refer to a file)
metadata_df <- data.frame(File_Name   = file_names,
                          File_Num    = seq_along(file_names),
                          parsed[, setdiff(colnames(parsed), "File_Name")],
                          Event_Count = event_counts,
                          stringsAsFactors = FALSE)



#### Write and save output metadata files ######################################
backup_file(metadata_path)
backup_file(panel_path)
fwrite(metadata_df, file = metadata_path)
fwrite(panel_df,    file = panel_path)

cat("Minimum Event Count =", min(event_counts),
    "\nMean Event Count =", round(mean(event_counts), 1),
    "\nMedian Event Count =", median(event_counts),
    "\nMaximum Event Count =", max(event_counts), "\n")



#### Prompt user to decide move to batch correction or not #####################
batchCorrection_res <- if (ask_yes_no("Will you be performing batch correction on this data?")) "yes" else "no"

Output_Text <- "You're such a meticulous scientist!!! Be sure to modify metadata.csv and panel.csv to fit your needs."
if (batchCorrection_res == "yes") {
  template_src  <- find_support_file(paste0("^", CHANNELS_TEMPLATE_NAME, "\\.txt$"), "DataCleanup")
  template_dest <- file.path(Original_WD, paste0(CHANNELS_TEMPLATE_NAME, ".txt"))
  if (length(template_src) >= 1) {
    if (length(template_src) > 1) {
      message("[INFO] ", length(template_src), " copies of ", CHANNELS_TEMPLATE_NAME,
              ".txt were found; using the first:\n  ", template_src[1])
    }
    backup_file(template_dest)
    file.copy(from = template_src[1], to = template_dest, overwrite = TRUE)
    Output_Text <- "You're such a meticulous scientist!!! Be sure to modify ChannelsToAdjust_Template.txt, metadata.csv, and panel.csv to fit your needs."
  } else {
    message("[INFO] ", CHANNELS_TEMPLATE_NAME, ".txt was not found under ",
            dirname(HELPERS_DIR), ". This is fine: 01.3_Update_ChannelsForBatchCorr ",
            "will build its starting list from panel.csv (Plot = 1 channels).")
    Output_Text <- "No channel template was found; 01.3 will build the batch-correction channel list from panel.csv. Modify metadata.csv and panel.csv to fit your needs."
  }
}



#### Keep objects in memory for the next scripts ###############################
metadata <- metadata_df
panel    <- panel_df
# Also in memory: Original_WD, batchCorrection_res, HELPERS_DIR



#### Record parameters used in this run ########################################
Inputs <- list(
  # --- keys that existed in v5 (names unchanged) ---
  Original_WD         = Original_WD,
  FCS_FILEPATH        = Original_WD,
  METADATA_FILENAME   = METADATA_NAME,
  PANEL_FILENAME      = PANEL_NAME,
  batchCorrection_res = batchCorrection_res,
  # --- new in v7/v8 ---
  Script_Version                = SCRIPT_VERSION,
  Completed_At                  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  Naming_Convention_Confirmed   = TRUE,
  Time_Label_Rule               = TIME_TOKEN_REGEX,
  TIMEPOINT_DIGITS              = TIMEPOINT_DIGITS,
  Conditions_Found              = as.list(table(parsed$Condition)),
  N_FCS_Files                   = length(file_list),
  FCS_Files                     = data.frame(File_Num    = seq_along(file_list),
                                             File_Name   = file_names,
                                             Event_Count = event_counts,
                                             Size_Bytes  = file.size(file_list),
                                             stringsAsFactors = FALSE),
  Panel_Consistent_Across_Files = panel_consistent,
  Random_Seed                   = "not used in this script",
  Output_MD5                    = file_md5(c(metadata_path, panel_path)),
  Session                       = session_info_list(names(pkg_versions))
)

# reset = TRUE: this is the first script, so start a new record (old one backed up)
record_params(MODULE_NAME, Inputs, write_dir = Original_WD, reset = TRUE)



#### Script completed message ##################################################
print_banner(" Metadata files output complete ", Output_Text)
