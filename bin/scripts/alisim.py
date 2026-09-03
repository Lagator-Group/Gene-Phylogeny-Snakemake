# alisim.py
# Simulation-based null for gene-gene phylogenetic congruence (Mantel r and normalized Robinson-Foulds).
# For a gene pair, simulate many alignments along gene2's tree under gene2's best-fit model (AliSim),
# re-infer a tree from each simulated alignment, then compare every simulated tree back to gene1's tree.
# The spread of the simulated statistics gives an empirical p-value for the observed gene1-vs-gene2 congruence.
# -
# - detect the iqtree binary the conda env provides (iqtree3 or iqtree)
# - parse gene1 + gene2 trees ONCE, prune both to their shared taxa ONCE, and reuse those
#   small trees everywhere (gene1 can have thousands of tips; re-parsing it per replicate is what
#   made different-gene pairs appear to hang)
# - read gene2's best-fit model + alignment length from its .iqtree file
# - simulate {num_bootstraps} alignments along gene2's tree with AliSim (one call)
# - re-infer a tree from every simulated alignment, each with its OWN prefix (so none overwrite)
# - observed + simulated-null Mantel r and normalized RF vs gene1's tree
# Input:  gene1/gene2 .iqtree + .treefile (from sfile_iqtree)
# Output: data/alisim/{gene1}-{gene2}.csv  (single row; sentinel NaN row if a metric can't be computed)

from ete3 import Tree
import pandas as pd
import numpy as np
from skbio.stats.distance import mantel, DistanceMatrix
import random, os, subprocess, sys, shutil, re
from pathlib import Path
import glob

sys.setrecursionlimit(1000000)


### Helpers: environment + iqtree output parsing
def find_iqtree_binary():
    # the env's binary may be 'iqtree3' (v3 suffix) or plain 'iqtree' - use whichever exists
    for _bin in ("iqtree3", "iqtree"):
        if shutil.which(_bin):
            return _bin
    # fail loudly: a missing binary was the silent killer before (capture_output hid it)
    raise RuntimeError(
        "No iqtree/iqtree3 binary on PATH. The alisim rule's conda env must include iqtree."
    )


def read_best_model(iqtree_path):
    # e.g. "Best-fit model according to BIC: CPREV+I"
    text = Path(iqtree_path).read_text()
    match = re.search(r"Best-fit model according to BIC:\s+(\S+)", text)
    if match is None:
        raise ValueError(f"Could not find a best-fit model line in {iqtree_path}")
    return match.group(1)


def read_num_sites(iqtree_path):
    # e.g. "Input data: 35 sequences with 406 amino-acid sites"
    text = Path(iqtree_path).read_text()
    match = re.search(r"with\s+(\d+)\s+amino-acid sites", text)
    if match is None:
        # not fatal: AliSim will fall back to its default length, just warn
        print(f"WARNING: could not read alignment length from {iqtree_path}; AliSim will use its default")
        return None
    return int(match.group(1))


### Tree comparison helpers (all operate on trees ALREADY pruned to the shared taxa)
def cophenetic_matrix(tree, labels):
    # pairwise patristic (branch-length) distance matrix over `labels`
    n = len(labels)
    mat = np.zeros((n, n), dtype=float)

    for i, a in enumerate(labels):

        for j, b in enumerate(labels):

            if i == j:
                mat[i, j] = 0.0

            else:
                d = tree.get_distance(a, b)

                if np.isnan(d):
                    d = 0.0

                mat[i, j] = d

    mat = (mat + mat.T) / 2      # force perfect symmetry
    np.fill_diagonal(mat, 0)     # enforce exact zero diagonal

    return mat


def prune_copy(tree, labels):
    # return an independent copy of `tree` pruned to `labels`, keeping branch lengths
    t = Tree(tree.write())
    t.prune(labels, preserve_branch_length=True)
    return t


def norm_rf(t_a, t_b):
    # normalized RF between two trees that share the same taxon set
    rf, max_rf, *_ = t_a.robinson_foulds(t_b, unrooted_trees=True)
    _nrf = rf / max_rf if max_rf > 0 else np.nan
    return _nrf, rf, max_rf


def cophenetic_dm(tree, labels):
    # unrooted cophenetic distance matrix as an skbio DistanceMatrix over `labels`
    t = Tree(tree.write())
    t.unroot()
    return DistanceMatrix(cophenetic_matrix(t, labels), ids=labels)


def randomize_leaf_labels(tree):
    # randomly permute tip labels while preserving topology (label-shuffle null)
    t = Tree(tree.write())

    leaves = t.get_leaves()
    labels = [leaf.name for leaf in leaves]

    shuffled = labels[:]
    random.shuffle(shuffled)

    for leaf, label in zip(leaves, shuffled):
        leaf.name = label

    return t


### Execution
#### snakemake imports
gene1_iqtree = snakemake.input.gene1_iqtree
gene2_iqtree = snakemake.input.gene2_iqtree
gene1_tree = snakemake.input.gene1_tree
gene2_tree = snakemake.input.gene2_tree
path_out = snakemake.output[0]

gene1 = snakemake.wildcards.gene1
gene2 = snakemake.wildcards.gene2

num_bootstraps = snakemake.params.num_bootstraps
min_shared_leaves = snakemake.params.min_shared_leaves
seed = snakemake.params.random_seed

print(f"{gene1} + {gene2}")

tmp_dir = f"data/tmp/{gene1}_{gene2}"
os.makedirs(tmp_dir, exist_ok=True)

# output columns (kept in one place so every exit path writes the same shape)
out_cols = [
    "Gene1", "Gene2",
    "Mantel r", "p", "n", "Mantel null mean", "Mantel null sd", "Mantel Sim p",
    "RF", "Max RF", "Normalized RF", "Shared taxa", "Permutation p-value",
    "nRF null mean", "nRF null sd", "nRF Sim p",
]

# initialise EVERYTHING to a sentinel up front, so building the row never NameErrors
mantel_r_obs = mantel_p = mantel_n = np.nan
mantel_null_mean = mantel_null_sd = mantel_sim_p = np.nan
rf_obs = maxrf_obs = nrf_obs = perm_p = np.nan
nrf_null_mean = nrf_null_sd = nrf_sim_p = np.nan


def write_row():
    # assemble the single-row output from whatever is currently set (sentinels for anything missing)
    row = pd.DataFrame([{
        "Gene1": gene1, "Gene2": gene2,
        "Mantel r": mantel_r_obs, "p": mantel_p, "n": mantel_n,
        "Mantel null mean": mantel_null_mean, "Mantel null sd": mantel_null_sd,
        "Mantel Sim p": mantel_sim_p,
        "RF": rf_obs, "Max RF": maxrf_obs, "Normalized RF": nrf_obs,
        "Shared taxa": n_shared, "Permutation p-value": perm_p,
        "nRF null mean": nrf_null_mean, "nRF null sd": nrf_null_sd, "nRF Sim p": nrf_sim_p,
    }], columns=out_cols)
    row.to_csv(path_out, index=False)


#### parse both trees ONCE and prune to shared taxa ONCE (the key performance fix)
t_g1 = Tree(gene1_tree)
t_g2 = Tree(gene2_tree)
shared = sorted(set(t_g1.get_leaf_names()) & set(t_g2.get_leaf_names()))
n_shared = len(shared)
print(f"gene1 tips: {len(t_g1)}  gene2 tips: {len(t_g2)}  shared: {n_shared}")

# guard: too few shared taxa -> write a sentinel row and stop (don't waste a 1000-alignment simulation)
if n_shared < min_shared_leaves:
    print(f"only {n_shared} shared taxa (< min_shared_leaves={min_shared_leaves}); writing sentinel row")
    write_row()
    shutil.rmtree(tmp_dir)
    sys.exit(0)

# small pruned trees, reused for EVERY comparison below (gene1 is parsed/pruned only once)
g1_shared = prune_copy(t_g1, shared)
g2_shared = prune_copy(t_g2, shared)
dm_g1 = cophenetic_dm(g1_shared, shared)   # gene1's cophenetic matrix - reused for every null Mantel


#### observed statistics (gene1 vs gene2), all on the small pruned trees
# Mantel r + its own permutation p (skbio)
dm_g2 = cophenetic_dm(g2_shared, shared)
mantel_r_obs, mantel_p, mantel_n = mantel(dm_g1, dm_g2, method="pearson", permutations=1000)
print(f"Mantel r = {mantel_r_obs:.4f}  p = {mantel_p:.4f}  n = {mantel_n}")

# normalized RF
nrf_obs, rf_obs, maxrf_obs = norm_rf(g1_shared, g2_shared)
print(f"RF = {rf_obs}  Max RF = {maxrf_obs}  nRF = {nrf_obs:.4f}")

# label-shuffle permutation p for RF (shuffle gene2 tips, recompute nRF vs gene1) - fast on pruned trees
random.seed(seed)
perm_null = []
for _ in range(num_bootstraps):
    g2_perm = randomize_leaf_labels(g2_shared)
    _nrf, _, _ = norm_rf(g1_shared, g2_perm)
    if not np.isnan(_nrf):
        perm_null.append(_nrf)
perm_p = (sum(x <= nrf_obs for x in perm_null) + 1) / (len(perm_null) + 1)
print(f"Permutation p-value = {perm_p:.4f}")


#### simulate alignments along gene2's tree, then re-infer a tree from each
iqtree_bin = find_iqtree_binary()
print(f"iqtree binary: {iqtree_bin}")

# simulate along gene2's tree, so use gene2's model + alignment length
model = read_best_model(gene2_iqtree)
n_sites = read_num_sites(gene2_iqtree)
print(f"model: {model}  length: {n_sites}")

sim_prefix = f"{tmp_dir}/sim"

# simulate along gene2 PRUNED TO THE SHARED TAXA, not the full gene2 tree. We only ever compare on
# the shared taxa, so simulating/re-inferring the full tree (which can be thousands of tips - e.g. a
# 9240-tip gene) is both wasteful and intractable: a single 9240-taxon ML re-inference can take hours.
# The pruned subtree keeps every re-inference small (n_shared tips) regardless of how big gene2 is.
sim_tree_path = f"{tmp_dir}/gene2_shared.nwk"
g2_shared.write(outfile=sim_tree_path, format=5)  # 29-taxon newick, branch lengths, no support labels

# 1) one AliSim call generates {num_bootstraps} alignments: sim_1.phy, sim_2.phy, ...
#    NOTE: pass a LIST with shell=False. subprocess.run(list, shell=True) drops every
#    argument after the program name on Linux, which is why AliSim produced nothing before.
alisim_cmd = [
    iqtree_bin,
    "--alisim", sim_prefix,
    "-t", sim_tree_path,
    "-m", model,
    "--num-alignments", str(num_bootstraps),
    "--seed", str(seed),
]
if n_sites is not None:
    alisim_cmd += ["--length", str(n_sites)]  # match the real alignment length

# print(" ".join(alisim_cmd))
_sim = subprocess.run(alisim_cmd, capture_output=True, text=True)  # no shell=True with a list!
if _sim.returncode != 0:
    print("=== AliSim FAILED ===")
    print(_sim.stdout[-1500:])
    print(_sim.stderr[-1500:])
    raise RuntimeError(f"AliSim failed for {gene2} (return code {_sim.returncode})")

# AliSim writes PHYLIP (.phy) by default; fall back across common extensions just in case
phy_files = []
for _ext in ("phy", "fa", "fasta", "phylip"):
    phy_files = sorted(glob.glob(f"{sim_prefix}*.{_ext}"))  # captured BEFORE any re-inference runs
    if phy_files:
        break
print(f"simulated {len(phy_files)} alignments")
if len(phy_files) == 0:
    print("tmp_dir contents:", os.listdir(tmp_dir))
    raise RuntimeError("AliSim wrote no alignment files - check the model string and flags")

# 2) re-infer a tree from each simulated alignment, each with a UNIQUE prefix so none overwrite
sim_treefiles = []
for _i, phy in enumerate(phy_files, 1):  # `phy` is already a full path - do NOT prepend tmp_dir again
    rep = Path(phy).stem                 # e.g. "sim_1" -> unique prefix per replicate
    rep_prefix = f"{tmp_dir}/{rep}"

    if _i == 1 or _i % 50 == 0:
        print(f"  re-inferring replicate {_i}/{len(phy_files)} for {gene1}-{gene2}")

    build_cmd = [
        iqtree_bin,
        "-s", phy,
        "-st", "AA",
        "-m", model,                     # fix the generating model (no ModelFinder per replicate = much faster)
        "-ntmax", str(snakemake.threads),
        "-fast",                         # fast search; 1000 full searches would otherwise be very slow
        "--prefix", rep_prefix,
        "--seed", str(seed),
        "-redo",
    ]
    # print(" ".join(build_cmd))
    _b = subprocess.run(build_cmd, capture_output=True, text=True)

    tf = f"{rep_prefix}.treefile"
    if _b.returncode == 0 and os.path.exists(tf):
        sim_treefiles.append(tf)
    else:
        print(f"  tree build failed for {phy} (rc={_b.returncode}) - skipping this replicate")
        # print(_b.stderr[-500:])

print(f"rebuilt {len(sim_treefiles)} simulated trees")


#### simulated-null distributions (each sim tree carries gene2's taxa, so prune to `shared` and
#### compare against the pre-pruned gene1 / its cached cophenetic matrix - no gene1 re-parsing)
mantel_null = []
nrf_null = []
for _j, tf in enumerate(sim_treefiles, 1):
    if _j == 1 or _j % 50 == 0:
        print(f"  scoring replicate {_j}/{len(sim_treefiles)} for {gene1}-{gene2}")

    try:                                   # one bad replicate must not sink the whole null
        t_sim = Tree(tf)
        sim_leaves = set(t_sim.get_leaf_names())

        keep = [x for x in shared if x in sim_leaves]  # normally == shared (sim has all gene2 tips)
        if len(keep) < 4:
            continue

        sim_shared = prune_copy(t_sim, keep)

        # if a tip is somehow missing from this sim tree, re-prune gene1 to match (rare); else reuse cache
        if len(keep) == n_shared:
            g1_ref = g1_shared
            dm_ref = dm_g1
            labels = shared
        else:
            g1_ref = prune_copy(g1_shared, keep)
            dm_ref = cophenetic_dm(g1_ref, keep)
            labels = keep

        # nRF vs gene1
        _nrf, _, _ = norm_rf(g1_ref, sim_shared)
        if not np.isnan(_nrf):
            nrf_null.append(_nrf)

        # Mantel r vs gene1 (permutations=0 -> point estimate only, fast)
        dm_sim = cophenetic_dm(sim_shared, labels)
        _r, _, _ = mantel(dm_ref, dm_sim, method="pearson", permutations=0)
        if not np.isnan(_r):
            mantel_null.append(_r)

    except Exception as e:
        print(f"  null replicate failed ({tf}): {e} - skipping")
        continue

# Mantel Sim p: congruence = HIGH r, so upper tail (how often the null reaches the observed r)
if len(mantel_null) > 0:
    _m = np.array(mantel_null)
    mantel_null_mean = _m.mean()
    mantel_null_sd = _m.std()
    mantel_sim_p = (np.sum(_m >= mantel_r_obs) + 1) / (len(_m) + 1)

# nRF Sim p: congruence = LOW nRF, so lower tail (how often the null is as congruent as observed)
if len(nrf_null) > 0:
    _r = np.array(nrf_null)
    nrf_null_mean = _r.mean()
    nrf_null_sd = _r.std()
    nrf_sim_p = (np.sum(_r <= nrf_obs) + 1) / (len(_r) + 1)

print(f"null sizes -> Mantel: {len(mantel_null)}  nRF: {len(nrf_null)}")
print(f"Mantel Sim p = {mantel_sim_p}   nRF Sim p = {nrf_sim_p}")


#### save + clean up
write_row()
shutil.rmtree(tmp_dir)


### Standalone test block (run outside Snakemake)
"""
class _S:
    pass
snakemake = _S()
snakemake.input = _S()
snakemake.input.gene1_iqtree = "data/iqtree_protein/P03066/P03066.iqtree"
snakemake.input.gene2_iqtree = "data/iqtree_protein/P0AFY6/P0AFY6.iqtree"
snakemake.input.gene1_tree = "data/iqtree_protein/P03066/P03066.treefile"
snakemake.input.gene2_tree = "data/iqtree_protein/P0AFY6/P0AFY6.treefile"
snakemake.output = ["data/alisim/P03066-P0AFY6.csv"]
snakemake.wildcards = _S(); snakemake.wildcards.gene1 = "P03066"; snakemake.wildcards.gene2 = "P0AFY6"
snakemake.params = _S(); snakemake.params.num_bootstraps = 20
snakemake.params.min_shared_leaves = 20; snakemake.params.random_seed = 1
snakemake.threads = 8
"""
