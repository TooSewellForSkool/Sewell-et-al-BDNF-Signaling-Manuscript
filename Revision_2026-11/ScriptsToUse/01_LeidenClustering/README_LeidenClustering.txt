================================================================================
README: LeidenClustering layer
Scripts: 03_Leiden.py, 04_UMAP.py, 04_UMAP_3D.py
         Pipeline2_03_Leiden.slurm, Pipeline2_04_UMAP.slurm,
         Pipeline2_04_UMAP_3D.slurm   (HPC only)
(The 03 / 04 numbering is historical; they are the first and second steps of
 this layer.)
================================================================================

WHAT THIS LAYER DOES
--------------------
Takes expression_matrix_analysis.csv (made by 02_FCS_to_CSV in the DataCleanup
layer) and produces
   03_Leiden      a cluster number for every cell
   04_UMAP        2D x/y coordinates for every cell (for plotting)
   04_UMAP_3D     3D x/y/z coordinates for every cell (optional). It always
                  lays out all cells in the matrix; the round number only
                  names the file (umap_xyz_3D_R<N>.csv).
The three scripts are independent of each other, so on the HPC you can submit
all three at the same time.

The same Python files are used locally and on the HPC. The only HPC-only files
are the three .slurm files, which tell the cluster how to run them.

Which R script put the files here? 02_FCS_to_CSV (and 02.1_SubsetAndExport)
copy these scripts into the data folder they create, so you run everything from
that folder. Run the scripts from the folder that contains
expression_matrix_analysis.csv.


THE BIG PICTURE (what the algorithm does, in plain words)
---------------------------------------------------------
1. For each cell, find its 15 most similar cells (nearest neighbours) in
   "channel space". Done with an approximate, fast method (HNSW), because exact
   search is too slow for ~1 million cells.
2. Turn those neighbour relationships into a weighted network ("fuzzy graph":
   cells are nodes, similar cells are linked, link strength says how similar).
3. 03_Leiden: cut the network into communities of tightly connected cells
   (the Leiden method). Each community = one cluster.
   04_UMAP:   arrange the network on a 2D (or 3D) plane so connected cells sit
   close together. This is only a picture; it is not used to define clusters.
Clusters (03) and layouts (04) are therefore built from the same neighbour
network but are calculated separately.


ROUNDS (clustering, then re-clustering part of the data)
--------------------------------------------------------
Round 1   all cells are clustered / laid out together.
Round 2+  you pick some of the round-1 clusters and cluster only those again,
          to split them further. You choose them in the groups file:

   cluster_R1_groups.csv  (made by round 1; open in Excel)
     R1               cluster numbers (do not edit)
     Cluster_Labels   your names for the clusters (optional, used later)
     Group_Assigns    REQUIRED before round 2. For every cluster enter
                        0        keep this cluster as it is (carry it through)
                        1, 2, 3  clusters with the same number are pooled and
                                 re-clustered together; use a different number
                                 for each pool
     Groups_Unique, Group_Labels   optional, used in later layers

   Then run round 2. It reads cluster_R1_assigns.csv and cluster_R1_groups.csv
   and writes cluster_R2_assigns.csv / cluster_R2_groups.csv. Round 3 works the
   same way from the round-2 files, and so on.

   In round-2+ output, clusters you carried through (0) are numbered first
   (1, 2, ...), followed by the new sub-clusters of group 1, then group 2, etc.

   04_UMAP in round 2+ makes one separate layout per group. Cells in group 0
   are NOT laid out: they are written as x = 0, y = 0, Group = 0. Remove
   Group = 0 cells before plotting round-2+ layouts.

   Always use the same round number for 03 and 04.

   Alternative to re-clustering all cells: 02.1_SubsetAndExport (R) can keep only
   chosen cells/clusters and write a new folder; its "continue numbering" option
   prepares that folder so you run the next round number there.


HOW TO RUN LOCALLY
------------------
Needs Python with these packages (install once, e.g. with conda or pip):
   pandas, numpy, leidenalg, igraph (python-igraph), umap-learn, hnswlib
In VS Code (or any terminal), go to the data folder and run:

   python 03_Leiden.py                  round 1
   python 03_Leiden.py --round 2        round 2
   python 04_UMAP.py --round 2
   python 04_UMAP_3D.py --round 2      (round is only a label in the file name)

If you do not give --round, the number written at the top of the script
(CLUSTER_ROUND = 1) is used. You do not need to edit the script for each round.
Large datasets take hours and several GB of memory.


HOW TO RUN ON THE HPC
---------------------
1. Copy the whole data folder (it contains the scripts, the .slurm files and
   expression_matrix_analysis.csv) to the cluster.
2. In that folder:
      sbatch Pipeline2_03_Leiden.slurm          round 1
      sbatch Pipeline2_04_UMAP.slurm
      sbatch Pipeline2_04_UMAP_3D.slurm         (round 1 label)
      sbatch Pipeline2_03_Leiden.slurm 2        round 2 (the 2 is the round)
      sbatch Pipeline2_04_UMAP.slurm 2
      sbatch Pipeline2_04_UMAP_3D.slurm 2
3. Check on a job with:  squeue -u $USER      Cancel with:  scancel <job id>
4. Progress is written to a log named <jobname>_<jobid>.out (errors in .err).
   Each job has its own log, so they no longer overwrite each other.

The settings you may need to edit are at the top of each .slurm file:
   --mail-user      your email (or leave the file alone and use
                    sbatch --mail-user=you@virginia.edu Pipeline2_03_Leiden.slurm)
   -A               the allocation (account) to charge: zunderlab or deppmannlab
   CONDA_ENV_NAME   the name of your conda environment; the full python path is
                    built from it using your own home folder
Everything else is usually left alone.

What the resource lines mean (for a non-expert):
   --partition=standard   which group of computers the job waits in
   -n 1                   "number of tasks" = how many separate programs SLURM
                          starts. You run one Python program, so 1 is right.
                          It does NOT choose how many CPU cores the program gets.
   --cpus-per-task        (not set) how many CPU cores that one program gets.
                          If not set, SLURM gives 1 core. Some of the packages
                          used here can use several cores, so asking for more
                          could make jobs faster, but it may also change the
                          neighbour search slightly (see "Reproducibility").
                          It is left unset on purpose. Ask HPC support or test
                          on a small dataset before changing it.
   --mem=50000            memory in MB (about 50 GB). If a job dies with
                          "out of memory" / "oom-kill" in the .err log, raise it.
   --time=72:00:00        longest the job may run (72 h). The job is stopped at
                          this limit, so raise it if a large run is cut off.
   --mail-type            emails at start, end and failure.
The .slurm files only load the module and call Python; you do not run Python on
the login node. The module name (miniforge) and partition are specific to your
cluster; if HPC support renames them, change them in the three .slurm files.


FILES IN / OUT
--------------
In:   expression_matrix_analysis.csv
      cluster_R<N-1>_assigns.csv and cluster_R<N-1>_groups.csv (round 2+)
Out:  03_Leiden:   cluster_R<N>_assigns.csv   one cluster number per cell
                   cluster_R<N>_groups.csv    cluster list for you to fill in
      04_UMAP:     umap_xy_R<N>.csv           x/y per cell (round 2+ also Group)
      04_UMAP_3D:  umap_xyz_3D_R<N>.csv       x/y/z per cell (N is only a label)
      each script: run_parameters_<script>_R<N>.json   the run record
Rows of all these files are in the same order as the cells in
expression_matrix_analysis.csv. If you re-run a step, the previous output files
are moved to a _backups folder (time-stamped) instead of being overwritten.


THE RUN RECORD
--------------
Every run writes a small JSON file (open in any text editor):
   Status, Started_At, Completed_At, Runtime_Seconds
   Inputs        each input file with a checksum (MD5), number of cells and the
                 channel names; for round 2+ the groups used
   Parameters    every setting that affects the result (listed below)
   Seeds         the random-number seeds used (see "Reproducibility")
   Upstream_Record  whether run_parameters.json from the R scripts was found and
                 whether expression_matrix_analysis.csv still matches the
                 checksum the R scripts recorded (a warning is printed if not)
   Outputs       each output file with its checksum
   Results       number of cells, number of clusters, cells per cluster
   Environment   Python version, versions of umap-learn, hnswlib, leidenalg,
                 igraph, numpy, pandas ..., the computer, and the SLURM job
Keep it with the results. It is the record of how they were made.


WHAT THE SETTINGS DO
--------------------
These are the original lab values. Changing any of them changes your clusters or
layouts, so leave them unless you have a reason and re-validate. All are near the
top of each script.

Round / files
  CLUSTER_ROUND          default round if --round is not given.
  CLUSTER_BASENAME, ASSIGNS_BASENAME, GROUPS_BASENAME, INPUT_EXPRS,
  UMAP_XY_BASENAME, UMAP_XYZ_OUTPUT     file names. Other layers expect them.

Neighbour search (HNSW) -- finds each cell's nearest neighbours quickly
  HNSW_SPACE = 'l2'      how distance is measured. 'l2' is squared straight-line
                         (Euclidean) distance between cells.
  HNSW_M = 16            how many links each cell keeps in the search index.
                         More = more accurate, more memory.
  HNSW_EF_CONSTRUCTION = 200   effort spent building the index. More = more
                         accurate, slower to build.
  HNSW_SET_EF = 20       effort spent on each search. More = more accurate,
                         slower. Must be at least UMAP_N_NEIGHBORS.

Neighbour network (used by 03, 04 and 04_3D)
  UMAP_N_NEIGHBORS = 15  how many neighbours define a cell's neighbourhood. The
                         most influential setting. Lower = finer detail, usually
                         more and smaller clusters; higher = broader structure.
  UMAP_METRIC = 'euclidean'    distance used when building the network.
  UMAP_ANGULAR = False   uses angle-based shortcuts only for angular metrics;
                         irrelevant for Euclidean.
  UMAP_SET_OP_MIX_RATIO = 1.0  how "A is near B" and "B is near A" are combined:
                         1 = either one is enough (union), 0 = both needed.
  UMAP_LOCAL_CONNECTIVITY = 1.0  how many nearest neighbours of each cell are
                         treated as certainly connected.
  UMAP_RANDOM_SEED       starting number for random steps (42 in 03 and 04,
                         99 in 04_3D).
  UMAP_VERBOSE = True    print UMAP's own progress messages.

Leiden clustering (03 only)
  LEIDEN_SEED = 1        starting number for random steps.
  LEIDEN_N_ITERATIONS = -1     repeat until the clusters stop improving.
  Method: modularity. There is no "resolution" setting; the number of clusters
  comes from the data and UMAP_N_NEIGHBORS.

Layout (04_UMAP and 04_UMAP_3D only)
  UMAP_N_COMPONENTS      2 (04_UMAP) or 3 (04_UMAP_3D) dimensions.
  UMAP_MIN_DIST = 0      how close points may sit. 0 lets cells pack tightly, so
                         clusters look compact (cell-to-cell distances in the
                         picture are not meaningful).
  UMAP_SPREAD = 1.0      overall scale of the picture; works together with
                         MIN_DIST.
  UMAP_N_EPOCHS = 1000   how long the layout is optimised. More = better
                         settled, slower.
  UMAP_INIT = 'spectral' starting positions come from the network itself, not
                         random placement.
  UMAP_INITIAL_ALPHA = 1.0     starting step size of the optimisation.
  UMAP_GAMMA = 1.0       strength of the push-apart force.
  UMAP_NEGATIVE_SAMPLE_RATE = 5  how many push-apart comparisons per pull-
                         together comparison.
  UMAP_PARALLEL = False  layout optimisation runs on one thread.
  UMAP_DENSMAP, UMAP_DENSMAP_KWDS, UMAP_OUTPUT_DENS, UMAP_OUTPUT_METRIC,
  UMAP_OUTPUT_METRIC_KWDS, UMAP_METRIC_KWDS, UMAP_EUCLIDEAN_OUTPUT
                         extra UMAP options left at the standard (off/default)
                         values.


REPRODUCIBILITY
---------------
Seeds are recorded in every run record: UMAP 42 (2D layout and graph), 99 (3D),
Leiden 1.
Even with the same seeds, a re-run can give slightly different results. The
nearest-neighbour search (hnswlib) uses several threads by default and, as with
most approximate-neighbour methods, the order in which threads add cells can
change the index. Different package versions can also change results. Therefore:
  * Keep the cluster_R*.csv and umap_*.csv files from any run used in a figure
    or paper; do not rely on re-running to get them back.
  * The package versions in the run record let you recreate the same software.
  * This behaviour is as in the original lab scripts and was left unchanged.


THINGS TO KNOW WHEN INTERPRETING OR REPORTING THE RESULTS
---------------------------------------------------------
These are inherited from the original scripts. They were NOT changed. They are
listed so that the methods can be described accurately. They have not been
tested here for how much they affect your results.
1. Edge weights. Leiden is run on a network whose link weight is
   (1 - fuzzy membership strength), i.e. stronger similarity gives a SMALLER
   weight (the original comment says this is "1-x instead of 1/x"). Most
   graph-clustering workflows weight links by similarity directly. A link with
   membership 1.0 gets weight 0. Describe the method as "Leiden (modularity) on
   a UMAP fuzzy graph with edge weights 1 - membership" rather than as the
   standard approach.
2. Distances. hnswlib with 'l2' returns squared Euclidean distances, and these
   are passed to UMAP's graph-building function together with
   metric='euclidean', which normally works with plain (not squared)
   distances. The neighbours found are the same, but the link strengths are
   computed from squared distances.
3. The 2D picture (04) and the clusters (03) come from separate calculations,
   so a cluster can look split or merged in the picture.
4. Round 2+ layouts exist only for the chosen groups; group 0 is at (0, 0).
5. Different seeds are used for Leiden (1), the 2D layout (42) and the 3D layout
   (99).


WHAT THE SCRIPTS CHECK BEFORE RUNNING
-------------------------------------
Each check stops the run with a plain-language message:
  - expression_matrix_analysis.csv missing, empty, non-numeric, or containing
    missing values; fewer cells than UMAP_N_NEIGHBORS
  - round 2+: previous-round files missing; a different number of cells than the
    expression matrix; a cluster that is not in the groups file; empty or
    non-numeric Group_Assigns; a group with too few cells; (04_UMAP) no group
    above 0
It only WARNS (and carries on) if the matrix no longer matches the checksum in
run_parameters.json, or if that file is missing.


COMMON PROBLEMS
---------------
"was not found in <folder>"
    You are not in the folder containing the matrix. Change into the data folder
    (locally: cd; HPC: submit with sbatch from that folder).
"has N rows but the expression matrix has M cells"
    The cluster files are from different data (e.g. before the data were
    subset or re-exported). Re-run the earlier round on this data, or use the
    folder made by 02.1 with "continue numbering".
"Group_Assigns is empty or not a number for cluster(s) ..."
    Fill in every row of the groups file (0 to keep a cluster as it is).
ModuleNotFoundError: No module named 'umap' (or hnswlib, leidenalg ...)
    The Python environment does not have that package, or the wrong Python was
    used. Locally, activate the environment; on the HPC check CONDA_ENV_NAME.
Job stops with "oom-kill" / "out of memory"
    Raise --mem in the .slurm file.
Job stops at the time limit
    Raise --time in the .slurm file.
Python path not found in the .err log
    CONDA_ENV_NAME in the .slurm file does not match an environment under
    ~/.conda/envs on the cluster (check with: ls ~/.conda/envs).
