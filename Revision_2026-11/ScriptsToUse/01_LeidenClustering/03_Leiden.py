# -*- coding: utf-8 -*-
"""
03_Leiden.py   --   Leiden clustering of the expression matrix
Original author: ezund (Zunder lab), 2021.   Script Version: 2

WHAT IT DOES
    Round 1:   clusters ALL cells (nearest neighbours -> fuzzy graph -> Leiden).
    Round >1:  re-clusters only the cluster GROUPS you chose in
               cluster_R<N-1>_groups.csv (column Group_Assigns); clusters given
               0 are carried through unchanged.

HOW TO RUN  (from the folder that contains expression_matrix_analysis.csv)
    python 03_Leiden.py               -> round set by CLUSTER_ROUND below
    python 03_Leiden.py --round 2     -> round 2 (overrides CLUSTER_ROUND)
    On the HPC:  sbatch Pipeline2_03_Leiden.slurm 2     (the 2 = round number)

INPUT   expression_matrix_analysis.csv  (+ cluster_R<N-1>_assigns.csv and
        cluster_R<N-1>_groups.csv when round > 1)
OUTPUT  cluster_R<N>_assigns.csv, cluster_R<N>_groups.csv,
        run_parameters_03_Leiden_R<N>.json (settings, seeds, versions)
See README_LeidenClustering.txt for what every setting means.

CHANGE LOG  (Version 1 = the lab's original script; Version 2 = this file)
  APPLIED -- do not change results:
    A1  Round number can be given on the command line (--round N), so the file
        does not need editing for each round (important on the HPC).
        CLUSTER_ROUND below is the default when nothing is given.
    A2  One script for local AND HPC. The HPC copy only differed in how umap
        was imported; `import umap.umap_ as umap_internal` works in both.
    A3  Removed repeated parameter definitions (same values were set 2-3 times),
        the unused INPUT_PANEL, the unused UMAP_* settings that this script never
        reads, and the `import umap` / `import umap.umap_ as umap` clash.
    A4  Leiden seed (1) and iterations (-1) are now named settings
        (LEIDEN_SEED, LEIDEN_N_ITERATIONS); the values are unchanged.
    A5  Round > 1 loops replaced by array indexing (also removes a loop whose
        variable names overwrote each other). Output verified identical.
    A6  Input checks that stop early with a clear message (see
        read_expression_matrix / load_round_inputs): missing file, non-numeric
        or missing values, assigns file with a different number of cells,
        blank Group_Assigns, clusters missing from the groups file, groups too
        small to cluster.
    A7  Existing output files are moved to _backups/ instead of overwritten.
    A8  Progress messages and a JSON run record (parameters, seeds, input and
        output MD5s, package versions, SLURM job info).
  NOT CHANGED: fuzzy_leiden() (the algorithm), all parameter values, the seeds,
        HNSW multithreading, and the 1 - fuzzy_set.data edge weights.
"""
## Import packages
import os
import sys
import json
import time
import shutil
import hashlib
import argparse
import platform
import datetime
import importlib.metadata
import pandas as pd
import numpy as np
import leidenalg
import igraph as ig
import umap.umap_ as umap_internal   # same import works on old and new umap-learn
import hnswlib

SCRIPT_NAME = "03_Leiden"
SCRIPT_VERSION = 2

# ------------------------------------------------------------------------------
# USER SETTING
# ------------------------------------------------------------------------------
# Default round, used only if --round is not given on the command line.
CLUSTER_ROUND = 1 # start with 1, then do 2 for subclustering, then do 3 for sub-subclustering, etc.

# Input/output file names
CLUSTER_BASENAME = "cluster_R" # For output files
ASSIGNS_BASENAME = "_assigns.csv" # For cluster assignment file
GROUPS_BASENAME = "_groups.csv" # For cluster grouping file
INPUT_EXPRS = "expression_matrix_analysis.csv" # Expression matrix dataset

# ------------------------------------------------------------------------------
# ALGORITHM SETTINGS (original values -- changing them changes the clusters)
# ------------------------------------------------------------------------------
# Parameters for HNSW approximate nearest neighbor determination
# Currently set to default values.  Not sure how changing these will effect results or runtime. . .
HNSW_SPACE = 'l2' #'l2' Squared L2: (sum((Ai-Bi)^2)), 'ip' Inner Product: (1.0 - sum(Ai*Bi)), or 'cosine' Cosine Similarity: 1.0 - sum(Ai*Bi) / sqrt(sum(Ai*Ai) * sum(Bi*Bi))
HNSW_EF_CONSTRUCTION = 200
HNSW_M = 16
HNSW_SET_EF = 20

# Parameters for the fuzzy graph (the same settings UMAP uses)
UMAP_N_NEIGHBORS = 15 # lower=local emphasis, higher=global emphasis
UMAP_RANDOM_SEED = 42 # See how much it changes from run to run!
UMAP_METRIC = 'euclidean' # Many other metrics possible.  Manhattan, Cosine, etc.
UMAP_ANGULAR = False
UMAP_SET_OP_MIX_RATIO = 1.0
UMAP_LOCAL_CONNECTIVITY = 1.0
UMAP_VERBOSE = True # Might want to turn this off some times?

# Parameters for Leiden clustering
LEIDEN_SEED = 1
LEIDEN_N_ITERATIONS = -1 # -1 = keep iterating until the clustering stops improving

# ==============================================================================
# SUPPORT CODE  (added in Version 2 -- NOT part of the clustering/embedding
# algorithm). Progress messages, input checks, backups and the run record.
# The same block is used in 03_Leiden.py, 04_UMAP.py and 04_UMAP_3D.py.
# ==============================================================================
START_TIME = datetime.datetime.now()
START_STAMP = START_TIME.strftime("%Y-%m-%d_%H-%M-%S")
UPSTREAM_RECORD = "run_parameters.json"   # written by the R scripts (DataCleanup layer)


def log(msg):
    # flush=True so messages show up immediately in SLURM log files
    print("[" + datetime.datetime.now().strftime("%H:%M:%S") + "] " + str(msg), flush=True)


def fail(msg):
    print("\nERROR: " + str(msg) + "\n", flush=True)
    sys.exit(1)


def file_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def pkg_version(name):
    try:
        return importlib.metadata.version(name)
    except Exception:
        return "not found"


def backup_existing(path):
    # Never overwrite earlier results: move them to _backups/ with a time stamp
    if os.path.exists(path):
        os.makedirs("_backups", exist_ok=True)
        dest = os.path.join("_backups", START_STAMP + "_" + os.path.basename(path))
        shutil.move(path, dest)
        log("Existing " + path + " moved to " + dest)


def read_upstream_record(expr_file):
    """Look for the R run record in this folder and check that the expression
    matrix is the one the R scripts wrote (MD5 comparison). Warns, never stops."""
    info = {"File": UPSTREAM_RECORD, "Found": False}
    if not os.path.exists(UPSTREAM_RECORD):
        log("NOTE: no " + UPSTREAM_RECORD + " in this folder (run record from the R scripts). "
            "Continuing, but the origin of the data cannot be traced.")
        return info
    info["Found"] = True
    info["MD5"] = file_md5(UPSTREAM_RECORD)
    try:
        with open(UPSTREAM_RECORD, "r", encoding="utf-8") as fh:
            rec = json.load(fh)
    except Exception as e:
        log("WARNING: could not read " + UPSTREAM_RECORD + " (" + str(e) + ")")
        return info
    info["Stage_Folder"] = rec.get("Folder_Info", {}).get("Stage_Folder")
    recorded = None
    for section in rec.values():          # last section that lists this file wins
        if isinstance(section, dict):
            md5s = section.get("Output_MD5")
            if isinstance(md5s, dict) and expr_file in md5s:
                recorded = md5s[expr_file]
    actual = file_md5(expr_file)
    info["Expression_MD5_In_Record"] = recorded
    info["Expression_MD5_Actual"] = actual
    if recorded is None:
        info["Expression_MD5_Matches"] = None
        log("NOTE: " + expr_file + " is not listed in " + UPSTREAM_RECORD + ".")
    else:
        info["Expression_MD5_Matches"] = (recorded == actual)
        if recorded == actual:
            log("Check passed: " + expr_file + " matches the R run record (MD5).")
        else:
            log("WARNING: " + expr_file + " does NOT match the MD5 in " + UPSTREAM_RECORD +
                ". The file was changed or replaced after the R scripts wrote it.")
    return info


def read_expression_matrix(expr_file, min_rows):
    """Read the expression matrix and stop early, with a clear message, if it
    is missing, not numeric, has missing values or is too small."""
    if not os.path.exists(expr_file):
        fail(expr_file + " was not found in " + os.getcwd() + ".\n"
             "Run this script from the folder that contains the expression matrix "
             "(the CompiledData / Subset folder made by the R scripts).")
    df = pd.read_csv(expr_file)
    if df.shape[1] == 0 or df.shape[0] == 0:
        fail(expr_file + " is empty.")
    non_numeric = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        fail(expr_file + " has non-numeric columns: " + ", ".join(map(str, non_numeric)) +
             "\nThe first row should be channel names and everything else numbers.")
    values = df.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        fail(expr_file + " contains missing (NA/NaN) or infinite values.")
    if df.shape[0] < min_rows:
        fail(expr_file + " has only " + str(df.shape[0]) + " cells; at least " + str(min_rows) +
             " are needed (UMAP_N_NEIGHBORS).")
    log("Read " + expr_file + ": " + str(df.shape[0]) + " cells x " + str(df.shape[1]) + " channels")
    return df


def load_round_inputs(cluster_round, n_cells, min_group_cells,
                      cluster_basename, assigns_basename, groups_basename):
    """Round > 1: read and check the previous round's assigns/groups files."""
    prev = cluster_round - 1
    prev_col = "R" + str(prev)
    assigns_file = cluster_basename + str(prev) + assigns_basename
    groups_file = cluster_basename + str(prev) + groups_basename
    for f in (assigns_file, groups_file):
        if not os.path.exists(f):
            fail(f + " was not found in " + os.getcwd() + ".\n"
                 "Round " + str(cluster_round) + " needs the files from round " + str(prev) +
                 ". Check the round number, or that you are in the right folder.")
    assigns_in = pd.read_csv(assigns_file)
    groups_in = pd.read_csv(groups_file)
    if prev_col not in assigns_in.columns:
        fail(assigns_file + " has no column named " + prev_col + ".")
    for col in (prev_col, "Group_Assigns"):
        if col not in groups_in.columns:
            fail(groups_file + " has no column named " + col + ".")
    if len(assigns_in) != n_cells:
        fail(assigns_file + " has " + str(len(assigns_in)) + " rows but the expression matrix has " +
             str(n_cells) + " cells.\nThey must be the same cells in the same order. This usually "
             "means the cluster files belong to a different dataset or an earlier version of the "
             "data (for example the data were subset or re-exported after clustering).")
    if groups_in[prev_col].duplicated().any():
        fail(groups_file + ": some clusters are listed more than once in column " + prev_col + ".")
    ga = pd.to_numeric(groups_in["Group_Assigns"], errors="coerce")
    blank = groups_in.loc[ga.isna(), prev_col].tolist()
    if blank:
        fail(groups_file + ": column Group_Assigns is empty or not a number for cluster(s) " +
             ", ".join(map(str, blank)) + ".\nEnter 0 for clusters that should be carried through "
             "unchanged, or a group number (1, 2, 3, ...) for clusters that should be re-clustered "
             "together. Every cluster needs an entry.")
    if ((ga < 0) | (ga != ga.round())).any():
        fail(groups_file + ": Group_Assigns must be whole numbers 0 or higher.")
    groups_in["Group_Assigns"] = ga.astype(int)
    unlisted = sorted(set(assigns_in[prev_col].unique()) - set(groups_in[prev_col]))
    if unlisted:
        fail(assigns_file + " contains cluster(s) not listed in " + groups_file + ": " +
             ", ".join(map(str, unlisted)))
    summary = {}
    for g in sorted(groups_in["Group_Assigns"].unique()):
        clusters = groups_in.loc[groups_in["Group_Assigns"] == g, prev_col].tolist()
        n = int(assigns_in[prev_col].isin(clusters).sum())
        summary[int(g)] = {"Clusters": [int(c) for c in clusters], "Cells": n}
        if g > 0 and n < min_group_cells:
            fail("Group " + str(g) + " has only " + str(n) + " cells; at least " +
                 str(min_group_cells) + " are needed (UMAP_N_NEIGHBORS). "
                 "Check Group_Assigns in " + groups_file + ".")
    log("Round " + str(cluster_round) + ": " + str(sum(1 for g in summary if g > 0)) +
        " group(s) to process; group 0 (carried through): " +
        str(summary.get(0, {}).get("Cells", 0)) + " cells")
    return assigns_in, groups_in, prev_col, summary


def environment_info():
    keys = ["SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_JOB_PARTITION", "SLURM_JOB_ACCOUNT",
            "SLURM_NTASKS", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURMD_NODENAME"]
    return {
        "Python": sys.version.split()[0],
        "Platform": platform.platform(),
        "Computer": platform.node(),
        "CPU_Count_Seen_By_Python": os.cpu_count(),
        "SLURM": {k: os.environ[k] for k in keys if k in os.environ},
        "Packages": {p: pkg_version(p) for p in
                     ["umap-learn", "hnswlib", "leidenalg", "igraph", "python-igraph",
                      "numpy", "pandas", "scipy", "numba", "scikit-learn", "pynndescent"]},
    }


def write_record(path, record):
    def _plain(o):
        return o.item() if hasattr(o, "item") else str(o)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2, default=_plain)
    log("Run record written: " + path)

# ------------------------------------------------------------------------------
# CORE ALGORITHM -- unchanged from the original (only constants named, see A4)
# ------------------------------------------------------------------------------
def fuzzy_leiden(umap_np):

    ## Use HNSW to find approximate nearest neighbors (ANN)
    p = hnswlib.Index(space = HNSW_SPACE, dim = umap_np.shape[1])
    data_labels = np.arange(umap_np.shape[0])
    p.init_index(max_elements = umap_np.shape[0],
                 ef_construction = HNSW_EF_CONSTRUCTION, M = HNSW_M)
    p.add_items(umap_np, data_labels)
    p.set_ef(HNSW_SET_EF) # ef should always be > k
    labels, distances = p.knn_query(umap_np, k = UMAP_N_NEIGHBORS)

    ## Convert ANN to fuzzy simplicial set (think of this like the knn UMAP uses)
    fuzzy_set = umap_internal.fuzzy_simplicial_set(
        X=umap_np,
        n_neighbors=UMAP_N_NEIGHBORS,
        random_state=np.random.RandomState(seed=UMAP_RANDOM_SEED),
        metric=UMAP_METRIC,
        metric_kwds={},
        knn_indices=labels,
        knn_dists=distances,
        angular=UMAP_ANGULAR,
        set_op_mix_ratio=UMAP_SET_OP_MIX_RATIO,
        local_connectivity=UMAP_LOCAL_CONNECTIVITY,
        verbose=UMAP_VERBOSE,
    )[0] # Just take the sparse matrix, not the extras

    ## Build igraph graph from fuzzy simplicial set sparse matrix
    leiden_graph = ig.Graph(np.column_stack((fuzzy_set.nonzero())), directed=False)
    weights = np.array(1-fuzzy_set.data) # do 1-x instead of 1/x because range is from 0-1
    leiden_graph.es['weight'] = weights

    ## Do Leiden clustering on fuzzy simplicial set graph
    leiden_out = leidenalg.find_partition(leiden_graph,
                                          leidenalg.ModularityVertexPartition,
                                          n_iterations=LEIDEN_N_ITERATIONS,seed=LEIDEN_SEED)

    # Add one to every cluster ID to make it 1-indexed instead of 0-indexed!
    plus_one = [x+1 for x in leiden_out.membership]
    return plus_one


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def get_cluster_round(default):
    parser = argparse.ArgumentParser(description="Leiden clustering (see header of this file).")
    parser.add_argument("--round", type=int, default=None,
                        help="cluster round (1 = all cells, 2 = subcluster round 1, ...)")
    args = parser.parse_args()
    value = args.round if args.round is not None else default
    if value < 1:
        fail("Cluster round must be 1 or higher (got " + str(value) + ").")
    return value


CLUSTER_ROUND = get_cluster_round(CLUSTER_ROUND)
RECORD_FILE = "run_parameters_" + SCRIPT_NAME + "_R" + str(CLUSTER_ROUND) + ".json"
assigns_outfile = CLUSTER_BASENAME + str(CLUSTER_ROUND) + ASSIGNS_BASENAME
groups_outfile = CLUSTER_BASENAME + str(CLUSTER_ROUND) + GROUPS_BASENAME

log(SCRIPT_NAME + " Version " + str(SCRIPT_VERSION) + " -- cluster round " + str(CLUSTER_ROUND) +
    " -- folder: " + os.getcwd())

record = {
    "Module": SCRIPT_NAME, "Script_Version": SCRIPT_VERSION, "Status": "Running",
    "Cluster_Round": CLUSTER_ROUND, "Started_At": START_TIME.strftime("%Y-%m-%d %H:%M:%S"),
    "Working_Folder": os.getcwd(), "Stage_Folder": os.path.basename(os.getcwd()),
}

try:
    expr_df = read_expression_matrix(INPUT_EXPRS, UMAP_N_NEIGHBORS)
    record["Upstream_Record"] = read_upstream_record(INPUT_EXPRS)
    record["Inputs"] = {INPUT_EXPRS: {"MD5": file_md5(INPUT_EXPRS), "Cells": expr_df.shape[0],
                                      "Channels": list(map(str, expr_df.columns))}}
    record["Parameters"] = {
        "HNSW_SPACE": HNSW_SPACE, "HNSW_EF_CONSTRUCTION": HNSW_EF_CONSTRUCTION,
        "HNSW_M": HNSW_M, "HNSW_SET_EF": HNSW_SET_EF,
        "UMAP_N_NEIGHBORS": UMAP_N_NEIGHBORS, "UMAP_METRIC": UMAP_METRIC,
        "UMAP_ANGULAR": UMAP_ANGULAR, "UMAP_SET_OP_MIX_RATIO": UMAP_SET_OP_MIX_RATIO,
        "UMAP_LOCAL_CONNECTIVITY": UMAP_LOCAL_CONNECTIVITY,
        "LEIDEN_N_ITERATIONS": LEIDEN_N_ITERATIONS,
        "Leiden_Partition_Type": "ModularityVertexPartition",
        "Edge_Weight": "1 - fuzzy membership strength",
    }
    record["Seeds"] = {
        "UMAP_RANDOM_SEED (fuzzy graph)": UMAP_RANDOM_SEED,
        "LEIDEN_SEED": LEIDEN_SEED,
        "HNSW random_seed": "hnswlib default (not set in this script; 100 in current hnswlib)",
        "Note": "hnswlib uses several threads by default, so neighbour search can differ slightly "
                "between runs/computers even with the same seeds.",
    }
    record["Environment"] = environment_info()
    backup_existing(RECORD_FILE)
    write_record(RECORD_FILE, record)

    # Convert to numpy array for input to UMAP
    umap_np = np.array(expr_df)

    if CLUSTER_ROUND > 1:
        # If round > 1, then use previous round clustering from cluster_ids file
        # and grouping from cluster metadata file to subset, and then do Leiden
        # clustering on each subset.
        assigns_in, groups_in, prev_col, group_summary = load_round_inputs(
            CLUSTER_ROUND, umap_np.shape[0], UMAP_N_NEIGHBORS,
            CLUSTER_BASENAME, ASSIGNS_BASENAME, GROUPS_BASENAME)
        record["Inputs"]["Groups_Used"] = {str(k): v for k, v in group_summary.items()}
        for f in (CLUSTER_BASENAME + str(CLUSTER_ROUND - 1) + ASSIGNS_BASENAME,
                  CLUSTER_BASENAME + str(CLUSTER_ROUND - 1) + GROUPS_BASENAME):
            record["Inputs"][f] = {"MD5": file_md5(f)}

        assigns_prev = assigns_in[prev_col].to_numpy()
        out_clusts = np.zeros(len(umap_np), dtype=int)
        add_for_unique = 0 #use this to clusters from separate subclustering rounds unique ids

        # Carry through 0-labeled clusters (don't subcluster them).
        # Rename them [0:n] to keep everything in order nicely.
        zero_clusts = np.array(groups_in.loc[groups_in['Group_Assigns'] == 0, prev_col])
        for c in range(len(zero_clusts)):
            out_clusts[assigns_prev == zero_clusts[c]] = c+1
        add_for_unique = len(zero_clusts)

        # Get all non-zero groups.  (Already took care of the zero group by carrying it through)
        all_groups = np.unique(groups_in['Group_Assigns'])
        nonzero_groups = all_groups[all_groups > 0]

        for g in nonzero_groups:
            g_clusts = np.array(groups_in.loc[groups_in['Group_Assigns'] == g, prev_col])
            in_group = np.isin(assigns_prev, g_clusts)
            log("Group " + str(g) + ": clustering " + str(int(in_group.sum())) + " cells ...")
            umap_np_sub = umap_np[in_group]
            clust_ids = fuzzy_leiden(umap_np_sub)
            clust_ids = np.array(clust_ids) + add_for_unique
            add_for_unique = np.max(clust_ids)
            out_clusts[in_group] = clust_ids
            log("Group " + str(g) + ": done (" + str(len(np.unique(clust_ids))) + " clusters)")
    else:
        # If round 1, then just do regular clustering
        log("Clustering " + str(len(umap_np)) + " cells ...")
        out_clusts = fuzzy_leiden(umap_np)

    ## Output cluster ids (existing files are moved to _backups first)
    backup_existing(assigns_outfile)
    backup_existing(groups_outfile)
    assigns_out = pd.DataFrame(out_clusts, columns=['R' + str(CLUSTER_ROUND)])
    assigns_out.to_csv(assigns_outfile, encoding='utf-8', index=False)

    ## Output cluster metadata - this file will have to be manually updated for cluster grouping
    group_vals = np.unique(out_clusts)
    groups_out = pd.DataFrame(group_vals, columns=['R' + str(CLUSTER_ROUND)])
    groups_fill_cols = pd.DataFrame(np.full([len(group_vals),4], ""),
                                   columns=['Cluster_Labels', 'Group_Assigns',
                                            'Groups_Unique', 'Group_Labels'])
    groups_out = pd.concat([groups_out, groups_fill_cols], axis=1)
    #groups_out = groups_out.replace(np.nan, '', regex=True)
    groups_out.to_csv(groups_outfile, encoding='utf-8', index=False)

    ## Run record
    ids, counts = np.unique(out_clusts, return_counts=True)
    record["Outputs"] = {assigns_outfile: {"MD5": file_md5(assigns_outfile)},
                         groups_outfile: {"MD5": file_md5(groups_outfile)}}
    record["Results"] = {"N_Cells": int(len(out_clusts)), "N_Clusters": int(len(ids)),
                         "Cells_Per_Cluster": {str(int(i)): int(n) for i, n in zip(ids, counts)}}
    record["Status"] = "Completed"
except SystemExit:
    record["Status"] = "Stopped by an input check (see the log)"
    raise
except BaseException as e:
    record["Status"] = "Failed: " + repr(e)
    raise
finally:
    record["Completed_At"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    record["Runtime_Seconds"] = round((datetime.datetime.now() - START_TIME).total_seconds(), 1)
    if "Environment" in record:        # only once the first record was written
        write_record(RECORD_FILE, record)

log("Finished. Cluster ids: " + assigns_outfile + " ; cluster list to edit: " + groups_outfile)
