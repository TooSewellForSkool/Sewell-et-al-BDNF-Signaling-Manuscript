================================================================================
README: DataCleanup layer
Scripts: 01_Get_Metadata, 01.2_Update_Metadata, 01.3_Update_ChannelsForBatchCorr,
         BatchCorrection/02_RunBatch (+ 03_BatchAdjust), 02_FCS_to_CSV,
         02.1_SubsetAndExport
================================================================================

WHAT THIS LAYER DOES
--------------------
Turns a folder of CyTOF FCS files into the CSV files used for Leiden clustering
(next layer, in Python), with a record of every choice made along the way.

    FCS files
       |
    01_Get_Metadata ............ reads the files; writes metadata.csv, panel.csv
       |
    01.2_Update_Metadata ....... you check/edit metadata and choose channels
       |
    (optional, only for batch correction)
    01.3_Update_ChannelsForBatchCorr ... choose channels to batch-adjust
    BatchCorrection/02_RunBatch ........ writes batch-corrected FCS files
       |
    02_FCS_to_CSV .............. writes the expression CSV files for clustering
       |
    Python clustering (next layer)
       |
    (optional, between clustering rounds)
    02.1_SubsetAndExport ....... keep chosen cells / channels for the next round

Run the scripts in this order. Each one checks that the earlier ones were run
and tells you if not.


ONE-TIME SETUP (applies to every script below)
----------------------------------------------
1. Folder layout. Keep this structure (the numbers in folder names may change):

     <project folder>
        <folder with your FCS files>     <- your WORKING FOLDER; all results go in it
        ScriptsToUse
           Helpers                      (00_Helpers.R -- do not move or edit)
           00_DataCleanup               (these scripts)
              BatchCorrection           (02_RunBatch.R and 03_BatchAdjust.R together)
           01_LeidenClustering          (Python and SLURM files)
           ...

   Every script finds Helpers by searching upward from where it is saved, and
   02_FCS_to_CSV / 02.1 find the clustering files by name inside ScriptsToUse.

2. R packages. Each script checks and tells you what is missing and how to
   install it. Altogether this layer uses:
     CRAN:         svDialogs, rstudioapi, data.table, jsonlite, shiny, miniUI,
                   rhandsontable           install.packages("packagename")
     Bioconductor: flowCore                install.packages("BiocManager")
                                           BiocManager::install("flowCore")

3. Run every script in RStudio, saved on disk (not an unsaved tab), by clicking
   "Source". Works on Mac and Windows. No XQuartz is needed: questions are
   normal pop-up windows (type a number to choose from a list) and tables open
   in your web browser.

4. Starting a NEW dataset? Restart R first (Session > Restart R) so nothing
   remembered from the previous dataset is used. Within one dataset, keep R open
   between scripts; if you close it, each script asks for the working folder and
   reads what it needs from the files.


HOW YOUR RESULTS ARE ORGANISED
------------------------------
* Original FCS files are never changed.
* 01, 01.2 and 01.3 write metadata.csv, panel.csv and ChannelsToAdjust.txt in
  your working folder. Anything they replace is first copied to _backups
  (time-stamped).
* Every later step (batch correction, CSV export, subsetting) writes a NEW
  folder directly inside your working folder, named
     <date>_<time>_<what it is>        e.g. 2026-10-21_15-43-30_BatchCorrection_quantile
                                            2026-10-22_01-19-02_CompiledData_FullDataSet_BatchCorrected-quantile
                                            2026-10-22_01-23-26_Subset_Files_<label>
  so nothing is overwritten and folders sort in the order they were made.
* run_parameters.json (the run record) is a plain-text file that lists every
  setting used. A new folder starts with a copy of the record of the folder it
  was made from and adds its own section, so the record in the last folder
  describes the whole history. Its "Folder_Info" section says which folder it
  came from. Keep it with your results; it is what makes a run repeatable.
* Random choices (subsampling in 02_FCS_to_CSV) use a seed that is saved in the
  record. Other scripts use no random numbers.
* Packages and versions (R, flowCore, ...) are recorded too.


WHICH SCRIPTS DO I NEED?
------------------------
  No batch correction, no clustering yet:   01, 01.2, 02_FCS_to_CSV
  With batch correction:                    01, 01.2, 01.3, 02_RunBatch, 02_FCS_to_CSV
  Another round of clustering on part of
  the data (or different channels):         ... 02_FCS_to_CSV, clustering, 02.1


================================================================================
THE SCRIPTS
================================================================================


================================================================================
SCRIPT 1: 01_Get_Metadata.R
================================================================================

WHAT THIS SCRIPT DOES
---------------------
Looks at the FCS files (CyTOF data) in your working folder and creates two
spreadsheets for you to check and edit in the next step (01.2_Update_Metadata):

  metadata.csv   one row per FCS file: sample information read from the file
                 name, plus the number of cells (events) in the file
  panel.csv      one row per channel (marker) in your antibody panel

It also saves run_parameters.json, a record of what was used, so the analysis
can be repeated exactly.

The script only reads the FILE HEADERS (plus the first file, to define the
panel), so it is fast and does not need much memory.


BEFORE YOU START
----------------
1. Folder layout. Keep this structure (names of the numbered folders can change):

     <project folder>
        <folder with your FCS files>        <- this is your WORKING FOLDER
        ScriptsToUse
           Helpers          (contains 00_Helpers.R -- do not move or edit)
           00_DataCleanup   (this script lives here)
           ...

2. R packages. The script checks for these and tells you what is missing:
     svDialogs, rstudioapi, data.table, jsonlite   (from CRAN)
     flowCore                                      (from Bioconductor)
   CRAN packages:      install.packages("packagename")
   flowCore:           install.packages("BiocManager")
                       BiocManager::install("flowCore")

3. Run in RStudio, with this script saved on disk (not an unsaved tab).
   Works on Mac and Windows.

4. Starting a NEW dataset? Restart R first (Session > Restart R). The scripts
   remember the working folder between steps; a restart makes sure they do not
   remember the previous dataset.

5. File names. Parts are separated by underscores:

     [BarcodeID]_[Normalization]_[Condition]_[TimePoint]_[Rep]_[Gating]_[Extra].fcs

   Examples:
     Sample001_normalized_Veh_8hr_1_CellsToAnalyze_EmptyGated.fcs
     Sample004_normalized_Cont._Depr._Test_1hr_1_4_Cells_To_Analyze_157Gated.fcs
     Sample001_normalized_Univ_1_CellsToAnalyze_EmptyGated.fcs     (universal sample)

   Rules the script relies on:
     - exactly ONE part contains "norm" (e.g. "normalized")
     - everything between "normalized" and the time label is the CONDITION. It
       may contain underscores (Cont._Depr.) and symbols (BDNF+K252a_Test)
     - the time label is a whole part like 8hr, 16hr, 15min or 1.5hr
       (minutes are converted to hours: 5min = 0.083 h)
     - files without a time label (Universal_1, HEK_Test_1) are allowed: the
       first part is the Condition and TimePoint is left blank (NA)
     - samples with "Univ" in the name are labelled Condition = Univ
     - the name ends in .fcs
   If a file does not follow these rules, the script stops BEFORE reading the
   data and lists the files that need attention.

6. All FCS files should come from the SAME antibody panel. The panel is read
   from the first file; you get a warning if the others differ.


HOW TO RUN
----------
1. Open 01_Get_Metadata.R in RStudio.
2. There is nothing you must edit. (The only setting, TIMEPOINT_DIGITS, controls
   rounding of minutes to hours and normally stays at 3.)
3. Click "Source" (top right of the script window).
4. Answer the pop-up questions:
     a) "Is this the folder that contains your FCS files?"  If not, click No and
        pick the folder.
     b) "Identified metadata components ... Are these correct?"  Read the list of
        CONDITIONS carefully. Every different treatment should appear once. If
        two things that should be separate are merged (or one is split in two),
        click No and tell the script owner.
     c) (Only if files from an earlier run exist) "Overwrite them?"  Click No
        to keep your edited files. Click Yes to start over; a timestamped copy
        of each old file is saved in the _backups folder first.
     d) "Will you be performing batch correction?"  Click Yes only if you will
        correct for batch effects; a template for choosing channels (if
        one exists in ScriptsToUse) is then copied into your working folder.
5. When you see "Metadata files output complete", go on to
   01.2_Update_Metadata.


WHAT YOU GET
------------
Written to your working folder (next to the FCS files):

  metadata.csv     Columns:
                     File_Name         the FCS file
                     File_Num          position of the file in the list (do not
                                       change; later steps use it)
                     SampleCollection  "Univ" for universal samples; otherwise
                                       "User-Defined" -- YOU check/fill in (01.2)
                     Expt_ID           placeholder -- YOU fill in (01.2)
                     BarcodeSet        first part of the file name
                     Treatment         Condition + time label, from the file name
                     Condition         from the file name (see rules above)
                     TimePoint         hours; blank (NA) if there is no time label
                     Replicate         1, 2, 3... counted within each Treatment,
                                       in file order
                     SampleType        placeholder -- YOU fill in (01.2)
                     Event_Count       number of cells (events) in the file

  panel.csv        Columns:
                     Name, Desc        channel name and marker as stored in the
                                       FCS file
                     Fixed_Param       cleaned-up marker name used downstream
                     Analysis          0/1, starts at 0 -- YOU set this in 01.2
                                       (1 = channel used for clustering)
                     Plot              0/1, starts at 0 -- YOU set this in 01.2
                                       (1 = channel included in plots)
                     Transform         Asinh for mass channels, otherwise Linear
                     Factor            5 for Asinh channels, otherwise 1

  run_parameters.json   what was used: file list with event counts, R and package
                        versions, date/time, checksums of the two CSV files.
                        Keep it with your results.

  ChannelsToAdjust_Template.txt   (only if you answered Yes to batch correction)

  _backups/        timestamped copies of anything that was replaced


IMPORTANT NOTES
---------------
* Do not rename metadata.csv or panel.csv; later scripts look for these names.
* File_Num follows the order your computer sorts the file names. Do not re-run
  this script on a different computer partway through an analysis.
* This script uses no random numbers, so no seed is needed.
* Re-running this script starts a new run_parameters.json (the old one is kept
  in _backups).
* Keep your R session open between scripts; they pass information to each other.
  If you close R, each later script can find what it needs from the files in the
  working folder.


COMMON PROBLEMS
---------------
"Could not find 00_Helpers.R"
    The Helpers folder is missing or this script was moved out of ScriptsToUse.

"Missing R package(s): ..."
    Install the packages listed (see BEFORE YOU START).

"Some FCS file names do not follow the expected naming convention"
    Rename the listed files, or ask for the script to be adapted to your naming.

"No .fcs files were found in: ..."
    Your working folder is not the folder with the FCS files. Restart R, run the
    script again and choose the right folder.

"The pop-up window was closed without an answer"
    Run the script again and answer with Yes or No.

Warning: "FCS header event counts could not be confirmed"
    The script is counting events by reading every file in full. The result is
    still correct, it just takes longer.

Warning: "Not every FCS file has the same channel names / markers as the first file"
    The files were not all run with the same panel. panel.csv reflects the
    first file only. Check that this is what you expect.

Message: "ChannelsToAdjust_Template.txt was not found"
    Not a problem. 01.3 builds its list from the channels marked Plot = 1.

"Stopped so that your existing files are not overwritten"
    You clicked No at the overwrite question. Nothing was changed.


================================================================================
SCRIPT 2: 01.2_Update_Metadata.R
================================================================================

WHAT THIS SCRIPT DOES
---------------------
Opens metadata.csv and panel.csv (created by 01_Get_Metadata) in an editor in
your web browser. You check them, fill in what is missing, choose which channels
to use, and save. Nothing is changed until you click "Done".


BEFORE YOU START
----------------
1. 01_Get_Metadata must have been run on this folder. The script checks this and
   tells you if it has not.
2. R packages (the script tells you what is missing):
     svDialogs, rstudioapi, data.table, jsonlite, shiny, miniUI, rhandsontable
   All install with install.packages("packagename").
3. Run in RStudio, with the script saved on disk. A web browser (Chrome, Safari,
   Edge...) is needed for the editor.
4. Keep the same R session open from 01_Get_Metadata if you can. If you closed
   R, the script asks for your working folder (the folder with your FCS files)
   and reads everything it needs from the files in it.


HOW TO RUN
----------
1. Open 01.2_Update_Metadata.R in RStudio. Nothing needs editing.
2. Click "Source".
3. An editor opens in your browser with two tabs.

   TAB 1 - Metadata Editor (one row per FCS file)
     Grey cells cannot be changed (File_Name, File_Num, Event_Count).
     Fill in / check:
       Expt_ID           experiment number (replace "User-Defined #")
       SampleType        pick Expt or Ctrl from the drop-down
       SampleCollection  "Univ" for universal samples, otherwise anything
                         useful to you
       Condition, Treatment, TimePoint, Replicate   read from the file names -
                         correct any that are wrong

   TAB 2 - Panel Editor (one row per channel)
     Name (grey) is the channel name in the FCS file and cannot be changed.
       Analysis    1 = channel is used for CLUSTERING; 0 = not
       Plot        1 = channel is included in PLOTS; 0 = not
       Fixed_Param marker name used in the later files; must be unique
       Transform   Asinh (mass channels) or Linear (e.g. Time, Width)
       Factor      divisor applied before the transform (5 for Asinh, 1 for
                   Linear)

   You can click a column header to sort the table while you work. The saved
   files always keep their original order.

4. Click "Done" (top right). The script checks your entries:
     - Problems that MUST be fixed (for example no Analysis channel, a Factor
       of 0, a repeated Fixed_Param) bring the editor back with your changes.
     - Reminders (for example SampleType or Expt_ID not filled in) ask whether
       you want to go back (Yes) or save anyway (No).
5. When you see "Cleanup complete", continue with:
     01.3_Update_ChannelsForBatchCorr   if you will do batch correction, or
     02_FCS_to_CSV                      if not.
   Clicking "Cancel" in the editor stops the script without saving anything.

You can run this script again at any time to change your choices. It always
starts from the files currently in the folder (including changes you made in
Excel).


WHAT YOU GET
------------
Updated in your working folder:
  metadata.csv, panel.csv     same columns as before, with your edits
  run_parameters.json         gains a "01.2_Update_Metadata" section: which
                              channels you chose for Analysis and Plot, which
                              columns you edited, SampleType counts, a
                              "Reviewed" time stamp and checksums of the files
  _backups/                   a time-stamped copy of the files before this run
                              (the first copy is exactly what 01 produced)


IMPORTANT NOTES
---------------
* Do not rename metadata.csv or panel.csv.
* The order of rows in panel.csv must match the order of the channels in the FCS
  files. The script keeps this order for you even if you sort the table.
* This script uses no random numbers, so no seed is needed.
* If you edit the CSV files in Excel instead, save them as CSV (comma
  delimited) and keep the column names unchanged.


COMMON PROBLEMS
---------------
"No ... was found" / "These earlier steps are not recorded as run in this folder"
    Run 01_Get_Metadata first, in the folder that contains your FCS files.

"Could not find 00_Helpers.R"
    The Helpers folder is missing or this script was moved out of ScriptsToUse.

"Missing R package(s)"
    Install the packages listed.

The editor does not open
    Look for a new browser tab. If R prints a web address (http://127.0.0.1:...),
    paste it into your browser.

"These problems must be fixed before the files can be saved"
    Read the list; correct the entries in the editor and click Done again.

"Editing was cancelled. No files were changed."
    You clicked Cancel. Run the script again when ready.


================================================================================
SCRIPT 3 (optional): 01.3_Update_ChannelsForBatchCorr.R
================================================================================

WHAT THIS SCRIPT DOES
---------------------
Batch correction adjusts the signal of selected channels so that barcode sets
(batches) can be compared. This script lets you confirm WHICH channels are
adjusted, and saves the list as ChannelsToAdjust.txt (one channel name per
line) in your working folder. The next step (BatchCorrection/02_RunBatch) uses
this file.

Skip this script if you are NOT doing batch correction; go straight to
02_FCS_to_CSV.


BEFORE YOU START
----------------
1. 01_Get_Metadata and 01.2_Update_Metadata must have been run on this folder.
   The script checks this and tells you if not.
2. You should have answered "Yes" to batch correction in 01_Get_Metadata. That
   (optionally) copied ChannelsToAdjust_Template.txt into your working folder.
   If you answered "No" the script asks whether to continue anyway.
   The starting list is the channels marked Plot = 1 in panel.csv (set in 01.2).
   The template is only used to show you where your list differs from it; it is
   used as the list only if fewer than 2 channels have Plot = 1.
3. R packages: svDialogs, rstudioapi, data.table, jsonlite, shiny, miniUI,
   rhandsontable (the script tells you what is missing).
4. Run in RStudio with the script saved on disk. A web browser is needed if you
   use the editor.


HOW TO RUN
----------
1. Open 01.3_Update_ChannelsForBatchCorr.R in RStudio. Nothing needs editing.
2. Click "Source".
3. A pop-up lists the channels that will be adjusted, grouped by whether they
   are marked Plot = 1 or Plot = 0 in panel.csv, and any listed channels that
   are not in panel.csv.
     Yes = adjust ALL of the channels listed. The list is saved as it is.
     No  = a table opens in your browser. Set Include to Y (adjust) or N (leave
           alone) for each channel; the starting list is pre-set to Y and every other metal channel to N. Click
           "Done" (top right) to save. "Cancel" stops without saving.
4. When you see "Cleanup complete", continue with BatchCorrection/02_RunBatch.

Which channels should be adjusted? Typically the channels measuring your
antibodies (not barcoding, DNA, bead or open channels). The template contains a
suggested list. Channel names must match the FCS files exactly (for example
Y89Di) - the table shows the names from panel.csv.


WHAT YOU GET
------------
In your working folder:
  ChannelsToAdjust.txt    one channel name per line (no header, no quotes)
  run_parameters.json     gains a "01.3_Update_ChannelsForBatchCorr" section:
                          the channels saved, where the starting list came from,
                          whether you accepted or edited it, and a checksum
  _backups/               a copy of any ChannelsToAdjust.txt that was replaced


IMPORTANT NOTES
---------------
* Do not rename ChannelsToAdjust.txt; 02_RunBatch looks for this name.
* Running the script again starts again from the Plot = 1 channels (not from
  your last choice), so re-check the list.
* This script uses no random numbers, so no seed is needed.


COMMON PROBLEMS
---------------
"These earlier steps are not recorded as run in this folder"
    Run 01_Get_Metadata and 01.2_Update_Metadata first.

"Nothing to do: batch correction was not selected"
    You answered No in 01 and No here. Continue with 02_FCS_to_CSV.

"Could not find 00_Helpers.R" / "Missing R package(s)"
    See README_01_Get_Metadata.txt.

Warning: "These channels are not in panel.csv but are included in ChannelsToAdjust.txt"
    The template lists channels your panel does not have. Click No at the pop-up
    to review the list, or ignore if you know the template is generic.

"No channels were selected"
    You set every channel to N. Run the script again.


================================================================================
SCRIPT 4 (optional): BatchCorrection/02_RunBatch.R (with 03_BatchAdjust.R)
================================================================================




================================================================================
SCRIPT 5: 02_FCS_to_CSV.R
================================================================================

WHAT THIS SCRIPT DOES
---------------------
Turns your FCS files into the CSV files used for clustering. It reads every FCS
file, optionally keeps a random subset of cells from each, applies the transforms
set in panel.csv (divide by Factor, then asinh for Asinh channels), and writes
the results to a NEW folder next to your FCS files.


BEFORE YOU START
----------------
1. 01_Get_Metadata and 01.2_Update_Metadata have been run on this folder (the
   script checks). If you did batch correction, 02_RunBatch has been run.
2. R packages: svDialogs, rstudioapi, data.table, jsonlite (CRAN) and flowCore
   (Bioconductor). The script tells you what is missing and how to install it.
3. Run in RStudio with the script saved on disk inside ScriptsToUse (the
   clustering Python files are found there).


HOW TO RUN
----------
1. Open 02_FCS_to_CSV.R in RStudio. Nothing needs editing.
2. Click "Source".
3. Answer the pop-ups:
   a. Which FCS files to convert (only if batch-corrected folders exist): type the
      number of the original files or of a corrected folder (newest first). The
      suggested number is already filled in.
   b. Full dataset? Yes = every cell. No = a random subset of cells per file.
        - Number of cells per file.
        - If some files have fewer cells than that, you can change the number or
          keep going (those files then contribute all their cells).
        - Random seed: leave blank for an automatic one, or type a number to
          repeat an earlier subsample exactly.
   c. Cluster locally or on an HPC (type 1 or 2): decides which scripts are copied.
   All choices are typed answers in a normal window (no XQuartz needed).
4. When you see "Expression files output complete", continue with Python
   (03_Leiden.py, 04_UMAP.py) in the new folder.


WHAT YOU GET
------------
A NEW folder inside your FCS folder named
  <date>_<time>_CompiledData_<FullDataSet | Subsampled-N>[_FilesWithFewerCells][_BatchCorrected-<method>]
containing:
  expression_matrix_analysis.csv   one row per cell; the channels with Analysis = 1
  expression_matrix_other.csv      the remaining channels
  filenums.csv                     File_Num of the file each cell came from
  subsample_idxs.csv               (subsampled runs) which cells were kept; one
                                   column per file, no header
  metadata.csv, panel.csv          copies from the folder the FCS files came from
  clustering scripts               copied from ScriptsToUse (LeidenClustering folder):
                                   03_Leiden.py, 04_UMAP.py, 04_UMAP_3D.py (same
                                   files locally and on the HPC); HPC also gets
                                   the three Pipeline2_*.slurm files.
                                   See 01_LeidenClustering/README_LeidenClustering.txt
  run_parameters.json              history of every earlier step + this one:
                                   input folder, subsampling, seed, channels in
                                   each file, cells kept per file, checksums


IMPORTANT NOTES
---------------
* Row order of panel.csv must match the channel order in the FCS files. The
  script checks this and stops if it does not.
* Columns are named with panel.csv's Fixed_Param.
* Running again makes another new folder; nothing is overwritten.
* With the same seed and the same input, a subsampled run gives the same cells.
* Need only some cells or samples for a later clustering round? That is the
  job of 02.1_SubsetAndExport.


COMMON PROBLEMS
---------------
"The channels in <file> do not match panel.csv"
    The panel does not belong to these FCS files, or its rows were reordered.
    Rerun 01_Get_Metadata (it rebuilds the panel) and then 01.2.
"<file> has N events but metadata.csv says M"
    metadata.csv is from different files. Rerun 01_Get_Metadata.
"panel.csv has an empty or non-numeric Factor"
    Set Factor in 01.2_Update_Metadata (5 for Asinh, 1 for Linear).
"Not found, so not copied: ..."
    The clustering Python files are not in ScriptsToUse. Copy them to the new
    folder by hand.


================================================================================
SCRIPT 6 (optional): 02.1_SubsetAndExport.R
================================================================================

WHAT THIS SCRIPT DOES
---------------------
Takes the cells you want to carry into the NEXT round of clustering and writes
them, with the channels you want to cluster on, to a new data folder.

    round N data --Python clustering--> cluster_R<N> files --02.1--> round N+1 data

Cells are selected by FCS file (condition, time point, ...) and/or by cluster or
cluster group from the clustering you already did. Cells you select are KEPT.
You can also use it before any clustering, just to drop files.


BEFORE YOU START
----------------
1. 02_FCS_to_CSV has been run (the script checks).
2. To select by cluster: Python clustering has been run in that data folder, so
   it contains cluster_R<N>_assigns.csv and cluster_R<N>_groups.csv. To select
   by GROUP, fill in the Group_Assigns column of the groups file first (Excel is
   fine; save as CSV). Group 0 means "carry through without sub-clustering".
3. R packages: svDialogs, rstudioapi, data.table, jsonlite, shiny, miniUI,
   rhandsontable (the script tells you what is missing).
4. Run in RStudio with the script saved on disk inside ScriptsToUse. A web
   browser is needed for the tables.


HOW TO RUN
----------
1. Open 02.1_SubsetAndExport.R in RStudio. Nothing needs editing.
2. Click "Source" and answer the steps (typed answers; no XQuartz):
   a. Which data folder to subset (type the number; newest first).
   b. Which cluster round to select from (type e.g. 1), or 0 to ignore clusters.
   c. Files: Yes opens a table; set Keep = N for files to leave out. Click Done.
   d. Clusters: the console lists cells per cluster and your groups. Then type
        group numbers (e.g. 1,3)   keep clusters in those groups
        all                        keep every cluster
        clusters                   then type cluster numbers, e.g. 1,4,7-9
   e. A summary shows cells before/after; confirm.
   f. Analysis channels: a table opens with the current Analysis and Plot flags.
      Set Analysis = 1 for the channels the next round should cluster on. Click Done.
   g. An optional label for the folder name, and local vs HPC (type 1 or 2).
3. When you see "Subset complete", go to Python. Open the clustering scripts in the
   new folder and set CLUSTER_ROUND = N + 1 (the banner tells you the number).


WHAT YOU GET
------------
A NEW folder inside your FCS folder named
  <date>_<time>_Subset_R<N>-to-R<N+1>[_label]      (or _Subset_Files)
containing:
  expression_matrix_analysis.csv   kept cells; the NEW Analysis channels
  expression_matrix_other.csv      kept cells; all other channels
  filenums.csv                     File_Num of each kept cell (not renumbered)
  metadata.csv                     files that still have cells; Event_Count is the
                                   number of cells now in the matrix
  panel.csv                        with the new Analysis / Plot flags
  cluster_R*_assigns.csv           earlier-round assignments for the kept cells
  cluster_R*_groups.csv            only clusters that still have cells (your group
                                   numbers and labels unchanged)
  umap_xy_R*.csv, umap_xyz_3D*.csv  earlier coordinates for the kept cells (if present)
  clustering scripts               for local or HPC, as in 02_FCS_to_CSV
  run_parameters.json              all earlier history + a "02.1_SubsetAndExport_R<N>toR<N+1>"
                                   section: files/clusters kept, cells before and
                                   after per file and per cluster, channels
                                   before and after, checksums
The input folder is never changed.


IMPORTANT NOTES
---------------
* Expression values are not recomputed: a channel has the same numbers as in the
  input folder. Only the split into analysis/other changes.
* subsample_idxs.csv is not carried over (it refers to the original files).
* File_Num keeps its original value so files can be traced across rounds.
* Repeating this step on a later round adds a new section to the history; it does
  not replace the earlier one.
* This script uses no random numbers, so no seed is needed.


COMMON PROBLEMS
---------------
"No folder with expression_matrix_analysis.csv was found"
    Run 02_FCS_to_CSV first, in the same FCS folder.
"... has N rows but filenums.csv has M"
    The cluster file does not belong to this data folder (for example it was
    made from a subsampled run). Re-run the Python clustering in this folder.
"Group_Assigns is empty or missing"
    Fill in Group_Assigns in the groups file, or choose clusters individually.
"The columns of the expression matrices do not match panel.csv"
    The folder was changed by hand. Re-create it with 02_FCS_to_CSV.
"No cells match your selection"
    The kept files and kept clusters have no cells in common.
