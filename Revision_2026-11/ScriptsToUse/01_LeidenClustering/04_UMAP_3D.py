# -*- coding: utf-8 -*-
"""
04_UMAP_3D.py   --   3D UMAP layout (x, y, z coordinates)
Original author: ezunder (Zunder lab), 2021.   Script Version: 2

WHAT IT DOES
    One 3-component UMAP of ALL cells in expression_matrix_analysis.csv
    -> umap_xyz_3D_R<N>.csv. The round number N is only a LABEL for the file
    name (so it sits with cluster_R<N> files); it does not change the
    calculation, and there are no per-group layouts as in 04_UMAP.py. If you
    subset the data with 02.1_SubsetAndExport, run it again in the subset folder.
    Uses its own seed (99) so it is independent of the 2D layout (seed 42).

HOW TO RUN  (from the folder that contains expression_matrix_analysis.csv)
    python 04_UMAP_3D.py              -> label set by CLUSTER_ROUND below
    python 04_UMAP_3D.py --round 2    -> umap_xyz_3D_R2.csv
    On the HPC:  sbatch Pipeline2_04_UMAP_3D.slurm 2     (the 2 = round label)

INPUT   expression_matrix_analysis.csv
OUTPUT  umap_xyz_3D_R<N>.csv, run_parameters_04_UMAP_3D_R<N>.json
See README_LeidenClustering.txt for what every setting means.

CHANGE LOG  (Version 1 = the lab's original script; Version 2 = this file)
  APPLIED -- do not change results:
    A2  One script for local AND HPC (the HPC copy only differed in how umap
        was imported; `import umap.umap_ as umap_internal` works in both).
    A3  Removed repeated parameter definitions (UMAP_N_EPOCHS was set twice)
        and the unused INPUT_PANEL.
    A6  Input checks that stop early with a clear message (missing file,
        non-numeric/missing values, too few cells).
    A7  An existing output file is moved to _backups/ instead of overwritten.
    A8  Progress messages and a JSON run record (parameters, seeds, input and
        output MD5s, package versions, SLURM job info).
    A9  (your edit, kept) output file name carries the cluster round
        (umap_xyz_3D_R<N>.csv). Added: --round on the command line, as in the
        other two scripts, and the same label in the run-record name.
  NOT CHANGED: run_umap() (the algorithm), all parameter values, the seed (99)
        and HNSW multithreading.
"""
import os
import sys
import json
import shutil
import hashlib
import argparse
import platform
import datetime
import importlib.metadata
import numpy as np
import pandas as pd
import umap.umap_ as umap_internal   # same import works on old and new umap-learn
import hnswlib

SCRIPT_NAME = "04_UMAP_3D"
SCRIPT_VERSION = 2

# ------------------------------------------------------------------------------
# USER SETTING
# ------------------------------------------------------------------------------
# Default round label, used only if --round is not given on the command line.
# Need to change this for the level of subsetting/subclustering!
CLUSTER_ROUND = 1

# Input/output file names
INPUT_EXPRS = "expression_matrix_analysis.csv"
UMAP_XYZ_BASENAME = "umap_xyz_3D_R"   # output is <basename><round>.csv

# ------------------------------------------------------------------------------
# ALGORITHM SETTINGS (original values -- changing them changes the layout)
# ------------------------------------------------------------------------------
# Parameters for HNSW approximate nearest neighbor determination
# Currently set to default values.  Not sure how changing these will effect results or runtime. . .
HNSW_SPACE = 'l2' #'l2' Squared L2: (sum((Ai-Bi)^2)), 'ip' Inner Product: (1.0 - sum(Ai*Bi)), or 'cosine' Cosine Similarity: 1.0 - sum(Ai*Bi) / sqrt(sum(Ai*Ai) * sum(Bi*Bi))
HNSW_EF_CONSTRUCTION = 200
HNSW_M = 16
HNSW_SET_EF = 20

# Parameters for UMAP embedding
UMAP_N_NEIGHBORS = 15 # lower=local emphasis, higher=global emphasis
UMAP_RANDOM_SEED = 99 # See how much it changes from run to run!
UMAP_MIN_DIST = 0
UMAP_N_COMPONENTS = 3 # Can try 3D if 2D has too much overlap
UMAP_N_EPOCHS = 1000 # Does running longer give a better-defined layout?
UMAP_METRIC = 'euclidean' # Cosine distance is less scale-dependent
UMAP_ANGULAR = False
UMAP_SET_OP_MIX_RATIO = 1.0
UMAP_LOCAL_CONNECTIVITY = 1.0
UMAP_VERBOSE = True # Might want to turn this off some times?
UMAP_SPREAD = 1.0 # Lowest spread (~1.0) generally looks "better" and less blobby
UMAP_INITIAL_ALPHA = 1.0
UMAP_GAMMA = 1.0
UMAP_NEGATIVE_SAMPLE_RATE = 5
UMAP_INIT = 'spectral'
UMAP_METRIC_KWDS = {}
UMAP_DENSMAP = False
UMAP_DENSMAP_KWDS = {}
UMAP_OUTPUT_DENS = False
UMAP_OUTPUT_METRIC = 'euclidean'
UMAP_OUTPUT_METRIC_KWDS = {}
UMAP_EUCLIDEAN_OUTPUT = True
UMAP_PARALLEL = False

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
# CORE ALGORITHM -- unchanged from the original (only the module name umap_internal)
# ------------------------------------------------------------------------------
def run_umap(umap_np):

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
    )

    ## Get the optimal a and b embedding values
    ## (regular UMAP function umap_internal.UMAP() takes care of this automatically,
    ## but simplicial_set_embedding doesn't, so we do it here instead)
    a, b = umap_internal.find_ab_params(UMAP_SPREAD, UMAP_MIN_DIST)

    ## Get the simplicial set embedding (this is like running UMAP on our own graph)
    data_embed = umap_internal.simplicial_set_embedding(
        data=umap_np,
        graph=fuzzy_set[0],
        n_components=UMAP_N_COMPONENTS,
        initial_alpha=UMAP_INITIAL_ALPHA,
        a=a,
        b=b,
        gamma=UMAP_GAMMA,
        negative_sample_rate=UMAP_NEGATIVE_SAMPLE_RATE,
        n_epochs=UMAP_N_EPOCHS,
        init=UMAP_INIT,
        random_state=np.random.RandomState(seed=UMAP_RANDOM_SEED),
        metric=UMAP_METRIC,
        metric_kwds=UMAP_METRIC_KWDS,
        densmap=UMAP_DENSMAP,
        densmap_kwds=UMAP_DENSMAP_KWDS,
        output_dens=UMAP_OUTPUT_DENS,
        output_metric=UMAP_OUTPUT_METRIC,
        output_metric_kwds=UMAP_OUTPUT_METRIC_KWDS,
        euclidean_output=UMAP_EUCLIDEAN_OUTPUT,
        parallel=UMAP_PARALLEL,
        verbose=UMAP_VERBOSE)[0] # Outputs more than just sparse matrix, so just take first element (sparse matrix)

    # ## This is the standard way of running umap, but it takes longer because it
    # ## doesn't let you import your own knn values.
    # data_embed = umap.UMAP(n_neighbors=UMAP_N_NEIGHBORS, min_dist=UMAP_MIN_DIST,
    #                  n_components=UMAP_N_COMPONENTS, n_epochs=UMAP_N_EPOCHS,
    #                  metric=UMAP_METRIC, ).fit_transform(umap_np)
    return data_embed

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def get_cluster_round(default):
    parser = argparse.ArgumentParser(description="3D UMAP layout (see header of this file).")
    parser.add_argument("--round", type=int, default=None,
                        help="cluster round label used in the output file name")
    args = parser.parse_args()
    value = args.round if args.round is not None else default
    if value < 1:
        fail("Cluster round must be 1 or higher (got " + str(value) + ").")
    return value


CLUSTER_ROUND = get_cluster_round(CLUSTER_ROUND)
UMAP_XYZ_OUTPUT = UMAP_XYZ_BASENAME + str(CLUSTER_ROUND) + '.csv'
RECORD_FILE = "run_parameters_" + SCRIPT_NAME + "_R" + str(CLUSTER_ROUND) + ".json"
log(SCRIPT_NAME + " Version " + str(SCRIPT_VERSION) + " -- cluster round label " +
    str(CLUSTER_ROUND) + " -- folder: " + os.getcwd())

record = {
    "Module": SCRIPT_NAME, "Script_Version": SCRIPT_VERSION, "Status": "Running",
    "Cluster_Round_Label": CLUSTER_ROUND, "Started_At": START_TIME.strftime("%Y-%m-%d %H:%M:%S"),
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
        "UMAP_N_NEIGHBORS": UMAP_N_NEIGHBORS, "UMAP_MIN_DIST": UMAP_MIN_DIST,
        "UMAP_N_COMPONENTS": UMAP_N_COMPONENTS, "UMAP_N_EPOCHS": UMAP_N_EPOCHS,
        "UMAP_METRIC": UMAP_METRIC, "UMAP_ANGULAR": UMAP_ANGULAR,
        "UMAP_SET_OP_MIX_RATIO": UMAP_SET_OP_MIX_RATIO,
        "UMAP_LOCAL_CONNECTIVITY": UMAP_LOCAL_CONNECTIVITY, "UMAP_SPREAD": UMAP_SPREAD,
        "UMAP_INITIAL_ALPHA": UMAP_INITIAL_ALPHA, "UMAP_GAMMA": UMAP_GAMMA,
        "UMAP_NEGATIVE_SAMPLE_RATE": UMAP_NEGATIVE_SAMPLE_RATE, "UMAP_INIT": UMAP_INIT,
        "UMAP_DENSMAP": UMAP_DENSMAP, "UMAP_EUCLIDEAN_OUTPUT": UMAP_EUCLIDEAN_OUTPUT,
        "UMAP_PARALLEL": UMAP_PARALLEL,
    }
    record["Seeds"] = {
        "UMAP_RANDOM_SEED": UMAP_RANDOM_SEED,
        "HNSW random_seed": "hnswlib default (not set in this script; 100 in current hnswlib)",
        "Note": "hnswlib uses several threads by default, so neighbour search can differ slightly "
                "between runs/computers even with the same seeds.",
    }
    record["Environment"] = environment_info()
    backup_existing(RECORD_FILE)
    write_record(RECORD_FILE, record)

    # Convert to numpy array for input to UMAP
    umap_np = np.array(expr_df)

    log("3D UMAP of " + str(len(umap_np)) + " cells ...")
    data_embed = run_umap(umap_np)
    # Extract x and y coordinates from UMAP output
    out_df = pd.DataFrame({'umap_x': data_embed[:, 0],'umap_y': data_embed[:, 1], 'umap_z': data_embed[:,2]})

    # Output umap XYZ coordinates; old file moved to _backups
    backup_existing(UMAP_XYZ_OUTPUT)
    out_df.to_csv(UMAP_XYZ_OUTPUT, index=False)

    record["Outputs"] = {UMAP_XYZ_OUTPUT: {"MD5": file_md5(UMAP_XYZ_OUTPUT)}}
    record["Results"] = {"N_Cells": int(len(out_df))}
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
    if "Environment" in record:
        write_record(RECORD_FILE, record)

log("Finished. UMAP coordinates: " + UMAP_XYZ_OUTPUT)
