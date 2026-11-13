################################################################################
# 00_Helpers.R -- shared helper functions for the DataCleanup layer
################################################################################
# Author:   Jonathon Sewell (functions consolidated with Claude)
# Version:  7   (HELPERS_VERSION below)
# Created:  2026-10-21
#
# WHAT THIS FILE IS
#   Functions that every DataCleanup script needs. Each script loads this file
#   automatically with the short "Load shared helper functions" block at its top,
#   so scripts can still be run on their own. You never run this file directly.
#
# WHERE IT LIVES
#   ScriptsToUse/Helpers/00_Helpers.R
#   Other layers can get their own file (01_Helpers.R, ...) in the same folder.
#
# FUNCTIONS (S = used by 01_Get_Metadata_v7; others are ready for later scripts
#            and have NOT been exercised yet -- they will be checked as each
#            script is converted)
#   dc_load_packages()   S  check + load packages, with install instructions
#   ask_yes_no()         S  yes/no pop-up that stops cleanly if the box is closed
#   print_banner()       S  centered title + wrapped text
#   timestamp_tag()      S  "2026-10-21_14-05-33" (safe in file names on any OS)
#   backup_file()        S  copy a file into <its folder>/_backups/ with a timestamp
#   get_working_dir()    S  Original_WD from memory, or ask the user
#   rel_path()              path relative to a base folder
#   make_stage_dir()        new timestamped output folder (sibling stage folders)
#   record_params()      S  write/merge run_parameters.json
#   read_run_record()    S  read run_parameters.json
#   require_modules()       sequence check: "has step X been run in this folder?"
#   get_record_value()      pull one value out of a run record
#   read_metadata()         read metadata.csv with the right column types
#   read_panel()            read panel.csv with the right column types
#   setup_seed()            choose/record a random seed
#   session_info_list()  S  R version, OS, locale, package versions for the record
#   file_md5()           S  checksums of output files
#   find_layer_dir()     S  locate a layer folder (e.g. "DataCleanup") by name
#   find_support_file()  S  locate a support file inside a layer folder
#   check_helpers_version() S
################################################################################

HELPERS_VERSION <- 7
RUN_RECORD_NAME <- "run_parameters.json"
METADATA_NAME   <- "metadata.csv"
PANEL_NAME      <- "panel.csv"


#### Basics ####################################################################
check_helpers_version <- function(required) {
  if (HELPERS_VERSION < required) {
    stop("This script needs 00_Helpers.R version ", required, " or newer, but ",
         "version ", HELPERS_VERSION, " was loaded from:\n  ", HELPERS_DIR,
         "\nCopy the newest 00_Helpers.R into that folder.", call. = FALSE)
  }
  invisible(TRUE)
}

# Check that packages are installed, then load them.
#   cran = packages from CRAN;  bioc = packages from Bioconductor
# Returns (invisibly) a named vector of the versions that were loaded.
dc_load_packages <- function(cran = character(), bioc = character()) {
  all_pkgs <- c(cran, bioc)
  is_missing <- !vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
  missing_cran <- all_pkgs[is_missing & all_pkgs %in% cran]
  missing_bioc <- all_pkgs[is_missing & all_pkgs %in% bioc]
  if (length(missing_cran) + length(missing_bioc) > 0) {
    msg <- paste0("Missing R package(s): ",
                  paste(c(missing_cran, missing_bioc), collapse = ", "), "\n")
    if (length(missing_cran) > 0) {
      msg <- paste0(msg, "\nInstall CRAN packages with:\n  install.packages(c(",
                    paste0("\"", missing_cran, "\"", collapse = ", "), "))\n")
    }
    if (length(missing_bioc) > 0) {
      msg <- paste0(msg, "\nInstall Bioconductor packages with:\n",
                    "  install.packages(\"BiocManager\")\n",
                    "  BiocManager::install(c(",
                    paste0("\"", missing_bioc, "\"", collapse = ", "), "))\n")
    }
    stop(msg, call. = FALSE)
  }
  for (p in all_pkgs) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  invisible(vapply(all_pkgs, function(p) as.character(utils::packageVersion(p)),
                   character(1)))
}

# Choice pop-up. Returns the chosen item of `choices` (stops if closed).
# [v4] Uses a plain typed-answer window (dlg_input) listing numbered options;
# svDialogs::dlg_list() needs XQuartz on macOS. Type the number (or the exact
# text). The window is shown again, up to 3 times, if the answer is not valid.
ask_choice <- function(message, choices, preselect = choices[1]) {
  # Only one option: ask for a Yes/No confirmation instead of a numbered list of one.
  if (length(choices) == 1) {
    if (ask_yes_no(paste0(message, "\n\nOnly one option was found:\n\n", choices))) return(choices)
    stop("Cancelled. Please run the script again.", call. = FALSE)
  }
  default_num <- match(preselect, choices)
  if (is.na(default_num)) default_num <- 1
  menu <- paste0(seq_along(choices), " = ", choices, collapse = "\n")
  for (attempt in 1:3) {
    res <- svDialogs::dlg_input(
      message = paste0(message, "\n\n", menu, "\n\nType the number of your choice:"),
      default = as.character(default_num))$res
    if (length(res) != 1 || is.na(res)) break
    res <- trimws(res)
    if (res %in% choices) return(res)
    num <- suppressWarnings(as.integer(res))
    if (!is.na(num) && num >= 1 && num <= length(choices)) return(choices[num])
    message <- paste0("'", res, "' is not one of the options.\n\n", message)
  }
  stop("No valid choice was entered. Please run the script again.", call. = FALSE)
}

# Yes/No pop-up. Returns TRUE for "yes", FALSE for "no".
# Stops with a plain message if the window is closed without answering.
ask_yes_no <- function(message) {
  res <- svDialogs::dlg_message(message, type = "yesno")$res
  if (length(res) != 1 || !(res %in% c("yes", "no"))) {
    stop("The pop-up window was closed without an answer. ",
         "Please run the script again.", call. = FALSE)
  }
  res == "yes"
}

# Centered title between lines of "=", with optional wrapped text underneath
print_banner <- function(title, lines = character(0), width = 72) {
  pad   <- max(width - nchar(title), 0)
  left  <- floor(pad / 2)
  right <- pad - left
  cat("\n", strrep("=", left), title, strrep("=", right), "\n", sep = "")
  if (length(lines) > 0) cat(paste0(strwrap(lines, width = width), "\n"), sep = "")
  cat(strrep("=", width), "\n", sep = "")
}

# "2026-10-21_14-05-33": date + time with hyphens (colons are not allowed in
# Windows file names), sorts correctly in a folder
timestamp_tag <- function() format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

# Copy a file to <its folder>/_backups/<name>_<timestamp>.<ext>
# (does nothing if the file does not exist)
backup_file <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  backup_dir <- file.path(dirname(path), "_backups")
  dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)
  ext <- tools::file_ext(path)
  backup_name <- paste0(tools::file_path_sans_ext(basename(path)), "_",
                        timestamp_tag(), if (nzchar(ext)) paste0(".", ext) else "")
  file.copy(path, file.path(backup_dir, backup_name), overwrite = TRUE)
  cat("[INFO] Saved backup: _backups/", backup_name, "\n", sep = "")
  invisible(file.path(backup_dir, backup_name))
}


#### Folders ###################################################################
# The "working directory" for the DataCleanup layer (called Original_WD) is the
# folder that holds your FCS files. It is kept in memory between scripts; if it
# is not there (new R session), the user is asked.
get_working_dir <- function(what = "fcs") {
  # what = "fcs"     DataCleanup layer: the folder that holds the FCS files
  # what = "project" Analysis / Visualization layers: the project folder that is searched (including its
  #                  sub-folders) for the output of the previous step
  is_fcs  <- identical(what, "fcs")
  q_text  <- if (is_fcs) "Is this the folder that contains your FCS files?" else
             paste0("Is this the correct working folder (project folder) for this analysis?\n",
                    "The script will look for the output of the previous step inside it and its sub-folders.")
  cap     <- if (is_fcs) "Select the folder that contains your FCS files" else
             "Select the working folder (project folder) for this analysis"
  cur <- get0("Original_WD", envir = .GlobalEnv, inherits = FALSE)
  if (!is.null(cur) && length(cur) == 1 && is.character(cur) && dir.exists(cur)) {
    cat("[INFO] Using the folder from the previous script:\n       ", cur, "\n", sep = "")
  } else {
    cat("\n[INFO] 'Original_WD' not found in memory. Requesting from user...\n")
    if (ask_yes_no(paste0(q_text, "\n\n", getwd()))) {
      cur <- getwd()
    } else {
      cat("[INFO] ", cap, "\n", sep = "")
      cur <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                 rstudioapi::isAvailable()) {
        rstudioapi::selectDirectory(caption = cap)
      } else {
        res <- svDialogs::dlg_dir(default = getwd(), title = cap)$res
        if (length(res) == 0 || !nzchar(res)) NULL else res
      }
      if (is.null(cur)) {
        stop("No folder was selected. Please run the script again and choose a folder.",
             call. = FALSE)
      }
    }
  }
  cur <- normalizePath(cur, winslash = "/", mustWork = TRUE)
  setwd(cur)
  assign("Original_WD", cur, envir = .GlobalEnv)
  cur
}

# Path of `path` relative to `base` ("." if identical; absolute if not inside it)
rel_path <- function(path, base) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  base <- sub("/+$", "", normalizePath(base, winslash = "/", mustWork = FALSE))
  if (identical(path, base)) return(".")
  prefix <- paste0(base, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1) else path
}

# Create a new output folder called <timestamp>_<label> inside parent_dir.
# Returns its full path. If the name is already taken, adds _2, _3, ...
make_stage_dir <- function(parent_dir, label) {
  clean_label <- gsub("[^A-Za-z0-9._+-]+", "-", label)
  base_name <- paste0(timestamp_tag(), "_", clean_label)
  path <- file.path(parent_dir, base_name)
  k <- 1
  while (dir.exists(path)) {
    k <- k + 1
    path <- file.path(parent_dir, paste0(base_name, "_", k))
  }
  dir.create(path, recursive = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# Find a layer folder by (part of) its name, ignoring number prefixes.
# Example: find_layer_dir("DataCleanup") finds ScriptsToUse/00_DataCleanup.
# Layer folders are searched for next to the Helpers folder.
find_layer_dir <- function(layer_name) {
  root <- dirname(HELPERS_DIR)
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  hit  <- dirs[grepl(layer_name, basename(dirs), fixed = TRUE)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# Find support files (templates, Python scripts, ...). [v2] Searches the layer
# folder first, then everything inside the folder that holds the Helpers folder
# (normally ScriptsToUse), so the file is found even if the layer folder has a
# different name than expected. Returns a (possibly empty) vector of full paths,
# best matches first, with duplicates removed. Backup and hidden folders
# (names starting with "_" or ".") are ignored.
find_support_file <- function(file_pattern, layer_name) {
  layer_dir <- find_layer_dir(layer_name)
  roots <- unique(c(if (!is.na(layer_dir)) layer_dir, dirname(HELPERS_DIR)))
  hits  <- character(0)
  for (r in roots) {
    found <- list.files(r, pattern = file_pattern, recursive = TRUE, full.names = TRUE)
    found <- found[!grepl("(^|/)[_.][^/]*/", sub(r, "", found, fixed = TRUE))]
    hits  <- c(hits, found)
  }
  unique(normalizePath(hits, winslash = "/", mustWork = FALSE))
}
# Copy the clustering scripts for the chosen environment into dest_dir.  [v6]
#   locally : 03_Leiden.py, 04_UMAP.py, 04_UMAP_3D.py (same Python files for HPC)
#   HPC     : the same three .py files + Pipeline2_03_Leiden.slurm, Pipeline2_04_UMAP.slurm,
#             Pipeline2_04_UMAP_3D.slurm (older *_UpdateUMAPImport files if those are absent)
# Files are found by name inside ScriptsToUse (a folder with "LeidenClustering"
# in its name is searched first). Returns list(copied, not_found).
copy_clustering_scripts <- function(dest_dir, environment = "locally") {
  # Current files: one set of Python scripts for local AND HPC; HPC also gets
  # the three Pipeline2_*.slurm files.
  py_files     <- c("03_Leiden.py", "04_UMAP.py", "04_UMAP_3D.py")
  slurm_files  <- c("Pipeline2_03_Leiden.slurm", "Pipeline2_04_UMAP.slurm",
                    "Pipeline2_04_UMAP_3D.slurm")
  # Older HPC-specific files (used only if the current slurm files are absent)
  legacy_hpc   <- c("03_Leiden_UpdateUMAPImport.py", "04_UMAP_UpdateUMAPImport.py",
                    "04_UMAP_3D_UpdateUMAPImport.py")
  copy_names <- function(nms) {
    copied <- character(0); missing <- character(0)
    for (nm in nms) {
      pat <- paste0("^", gsub(".", "\\.", nm, fixed = TRUE), "$")
      src <- find_support_file(pat, "LeidenClustering")
      if (length(src) == 0) { missing <- c(missing, nm); next }
      file.copy(src[1], file.path(dest_dir, nm), overwrite = TRUE)
      copied <- c(copied, nm)
    }
    list(copied = copied, not_found = missing)
  }
  if (environment == "HPC") {
    have_new_slurm <- all(vapply(slurm_files, function(nm)
      length(find_support_file(paste0("^", gsub(".", "\\.", nm, fixed = TRUE), "$"),
                               "LeidenClustering")) > 0, logical(1)))
    if (have_new_slurm) {
      res <- copy_names(py_files)
      sres <- copy_names(slurm_files)
      res$copied <- c(res$copied, sres$copied)
    } else {
      res <- copy_names(legacy_hpc)
      slurm <- unique(basename(find_support_file("UpdateUMAPImport.*\\.slurm$", "LeidenClustering")))
      if (length(slurm) > 0) {
        sres <- copy_names(slurm)
        res$copied <- c(res$copied, sres$copied)
      } else {
        res$not_found <- c(res$not_found, slurm_files)
      }
    }
  } else {
    res <- copy_names(py_files)
  }
  if (length(res$not_found) > 0) {
    cat("[INFO] Not found, so not copied: ", paste(res$not_found, collapse = ", "),
        "\n       (they should be inside ScriptsToUse, in the LeidenClustering folder)\n", sep = "")
  }
  res
}


#### Run record (run_parameters.json) ##########################################
# One JSON file per folder. Each script adds (or replaces) its own section, so
# the file shows exactly which settings produced the data in that folder.
# New stage folders start from a copy of their parent's record (parent_dir), so
# the full history travels with the data.

read_run_record <- function(dir) {
  p <- file.path(dir, RUN_RECORD_NAME)
  if (!file.exists(p)) return(NULL)
  tryCatch(jsonlite::read_json(p, simplifyVector = FALSE),
           error = function(e) {
             warning("Could not read ", p, " (", conditionMessage(e), ")",
                     call. = FALSE)
             NULL
           })
}

# module_name: e.g. "01_Get_Metadata";  Inputs: a named list of what to record
# write_dir:   folder that gets the record
# reset:       TRUE = start a new record and back up the old one (script 01 only)
# parent_dir:  folder this one was made from (copies its record the first time)
record_params <- function(module_name, Inputs, write_dir, reset = FALSE,
                          parent_dir = NULL) {
  json_path <- file.path(write_dir, RUN_RECORD_NAME)
  params <- list()
  if (file.exists(json_path)) {
    if (reset) {
      backup_file(json_path)
    } else {
      existing <- read_run_record(write_dir)
      if (!is.null(existing)) params <- existing
    }
  } else if (!is.null(parent_dir)) {
    parent_record <- read_run_record(parent_dir)
    if (!is.null(parent_record)) params <- parent_record
  }

  params[["Folder_Info"]] <- NULL
  params[[module_name]]   <- Inputs
  folder_info <- list(
    Stage_Folder         = basename(normalizePath(write_dir, winslash = "/")),
    Parent_Folder        = if (is.null(parent_dir)) "none" else basename(parent_dir),
    # [v3] sibling folders are written as "../<folder>", so the record still
    # works if the whole project folder is moved
    Parent_Path_Relative = if (is.null(parent_dir)) "none" else
      if (identical(normalizePath(dirname(write_dir), winslash = "/"),
                    normalizePath(parent_dir, winslash = "/"))) ".." else
        if (identical(normalizePath(dirname(write_dir), winslash = "/"),
                      normalizePath(dirname(parent_dir), winslash = "/")))
          paste0("../", basename(parent_dir)) else
          normalizePath(parent_dir, winslash = "/"),
    Last_Updated         = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  params <- c(list(Folder_Info = folder_info), params)

  assign("PipelineParams", params, envir = .GlobalEnv)  # kept for older scripts
  jsonlite::write_json(params, json_path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null", digits = NA)
  invisible(params)
}

# Sequence check. Stops with a plain-language message unless every module in
# `modules` has a section in this folder's run record.
#   require_modules(Original_WD, c("01_Get_Metadata", "01.2_Update_Metadata"))
require_modules <- function(dir, modules, how_to_fix = "") {
  rec <- read_run_record(dir)
  if (is.null(rec)) {
    stop("No ", RUN_RECORD_NAME, " was found in:\n  ", dir, "\n",
         "The earlier scripts have not been run in this folder. ",
         "Start with 01_Get_Metadata.", call. = FALSE)
  }
  not_run <- setdiff(modules, names(rec))
  if (length(not_run) > 0) {
    stop("These earlier steps are not recorded as run in this folder:\n  ",
         paste(not_run, collapse = "\n  "), "\nFolder: ", dir,
         "\nPlease run them first, in order. ", how_to_fix, call. = FALSE)
  }
  invisible(rec)
}

# rec[[module]][[key]], or `default` if it is not there
get_record_value <- function(rec, module, key, default = NULL) {
  v <- rec[[module]][[key]]
  if (is.null(v)) default else v
}

# Checksums (MD5) of files, as a named list for the run record
file_md5 <- function(paths) {
  stats::setNames(as.list(unname(tools::md5sum(paths))), basename(paths))
}

# R version, OS, locale and package versions for the run record
session_info_list <- function(pkgs = character()) {
  os <- utils::sessionInfo()$running
  if (is.null(os)) os <- Sys.info()[["sysname"]]
  list(R_Version        = R.version.string,
       Platform         = R.version$platform,
       OS               = os,
       Collation_Locale = Sys.getlocale("LC_COLLATE"),
       Package_Versions = as.list(vapply(pkgs, function(p)
         as.character(utils::packageVersion(p)), character(1))))
}


#### Reading the pipeline's CSV files ##########################################
# These make sure columns have the same type every time, no matter which
# script (or R session) reads them.

read_metadata <- function(dir) {
  p <- file.path(dir, METADATA_NAME)
  if (!file.exists(p)) stop(METADATA_NAME, " not found in:\n  ", dir,
                            "\nRun 01_Get_Metadata first.", call. = FALSE)
  md <- read.csv(p, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("File_Name", "File_Num")
  if (!all(need %in% colnames(md))) {
    stop(METADATA_NAME, " is missing column(s): ",
         paste(setdiff(need, colnames(md)), collapse = ", "), call. = FALSE)
  }
  md$File_Num <- as.integer(md$File_Num)
  md
}

read_panel <- function(dir) {
  p <- file.path(dir, PANEL_NAME)
  if (!file.exists(p)) stop(PANEL_NAME, " not found in:\n  ", dir,
                            "\nRun 01_Get_Metadata first.", call. = FALSE)
  pn <- read.csv(p, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("Name", "Desc", "Fixed_Param", "Analysis", "Plot", "Transform", "Factor")
  if (!all(need %in% colnames(pn))) {
    stop(PANEL_NAME, " is missing column(s): ",
         paste(setdiff(need, colnames(pn)), collapse = ", "), call. = FALSE)
  }
  pn$Analysis <- as.integer(pn$Analysis)
  pn$Plot     <- as.integer(pn$Plot)
  pn$Factor   <- as.numeric(pn$Factor)
  pn
}


#### Random numbers ############################################################
# Use for any step that samples cells. Leave `seed` as NULL to get a new random
# seed each time; type a number to repeat an earlier run exactly. Returns a list
# for the run record. (Set the seed BEFORE the sampling happens.)
setup_seed <- function(seed = NULL) {
  if (is.null(seed) || length(seed) != 1 || is.na(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
    seed_source <- "auto (random)"
  } else {
    seed <- as.integer(seed)
    seed_source <- "user-defined"
  }
  set.seed(seed)
  list(Random_Seed = seed,
       Seed_Source = seed_source,
       RNG_Kind    = paste(RNGkind(), collapse = " / "))
}
