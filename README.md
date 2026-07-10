# Plasmid Transcription Factor Phylogeny
The available Snakemake pipelines align sequences with Muscle, build trees with IQ-TREE, compare the phylogenies of two genes, and test whether a gene tree's structure is associated with sample metadata traits.

Pipelines:
- **sfile_muscle**: Aligns sequences from a single file using [Muscle](https://github.com/rcedgar/muscle).
- **sfile_iqtree**: Constructs gene phylogenetic trees with [IQ-TREE](https://github.com/Cibiv/IQ-TREE).
- **sfile_geneGene**: Compares the phylogenies of two genes from their `.treefile`s (Mantel + Robinson-Foulds).
- **sfile_geneTrait**: Tests whether a gene tree's structure is associated with metadata traits. Handles both categorical and continuous traits (see below).

## Gene–Trait analysis (`sfile_geneTrait`)
A single run of `sfile_geneTrait` produces two independent summaries, one per trait type. Which columns feed which branch is controlled by two keys in `config.yml`.

### Categorical / discrete traits — `trait_col`
For each unique value of the trait, `bin/scripts/gene_trait.py` binarises the tree tips (value vs. not-value) and runs three tests: a tree-purity permutation test, a Fitch parsimony-score permutation test, and a binary Mantel test. Results are merged and Benjamini-Hochberg corrected within each (gene, trait).
- Output: **`data/gene_trait_summary.csv`** (one row per gene × trait × value).
- Best for nominal labels of usable cardinality (e.g. `predicted_host_range_overall_name`, `predicted_mobility`, `Country`).
- Note: cells with comma-separated multi-values (e.g. `AMR`, `rep_type(s)`, `Inc_group`) are tested as exact-combination categories — explode them first if you want per-item results.

### Continuous / numeric traits — `trait_num_col`
`bin/scripts/gene_trait_numeric.py` treats the column as a single continuous variable and measures phylogenetic signal with four complementary statistics: **Blomberg's K**, **Pagel's λ**, **Moran's I**, and a **continuous Mantel test** (cophenetic distance vs. pairwise trait difference). Results are merged and BH-corrected within each trait.
- Output: **`data/gene_trait_numeric_summary.csv`** (one row per gene × trait).
- Best for interval/ratio variables where numeric distance is meaningful (e.g. `Last_update_year`, `predicted_host_range_count`, `LOCATION_lat`/`LOCATION_lng`).
- Non-numeric cells are dropped; duplicate accessions are averaged to one value; each accession is a single tip (no pseudo-replication).
- Caveats: numeric-looking ID columns (NCBI taxids) are nominal, not continuous — keep them out of `trait_num_col`. Cyclic numerics like `Last_update_month` (Dec≈Jan) belong in `trait_col`.

## Interpreting the output
Both summaries answer one question: **does the trait cluster on the gene tree more than expected by chance?** A significant result means the gene's phylogeny carries information about the trait (tips with similar trait values tend to sit near each other). Every test is backed by a permutation/likelihood null, so the **p-value is what decides significance** — the effect-size statistics (r, K, λ, I, purity) tell you the *direction and strength*, but a large effect with a non-significant p-value (often from small N) is not evidence of anything.

**Read the p-values first, and prefer the BH-corrected ones.** Because many values/traits are tested at once, use the `... p-value BH` (q-value) columns rather than the raw p-values for your final call. Convention: **q < 0.05 = significant**. The `Significant BH (q<0.05)` boolean is TRUE when *any* test for that row passes. Raw p-values below ~0.05 that don't survive BH correction should be treated as suggestive at best.

**Always filter on `Status` and N first.** Only `Status == ok` rows have valid statistics. `error` (e.g. zero variance, or a category with too few members) and `insufficient_data` (fewer than the minimum accessions) rows carry no result. Very small `N_Accessions` makes even a large effect statistically meaningless.

### Categorical traits (`gene_trait_summary.csv`)
| Statistic | Range | Strong signal (interesting) | No signal (not interesting) |
|---|---|---|---|
| **Observed Purity** | ~0.5–1.0 | Higher than **Null Mean** (tips of one value share clades) | ≈ Null Mean |
| **Null Mean** | ~0.5–1.0 | — (baseline: mean purity of random tip shuffles) | — |
| **Purity p-value** (+ BH) | 0–1 | **< 0.05**: clustering beats chance | **> 0.05** |
| **Parsimony Score** | integer ≥1 | **Lower** than **Parsimony Null Mean** (few state changes needed = same-value tips clustered) | ≈ Parsimony Null Mean |
| **Parsimony Null Mean** | integer | — (baseline: mean parsimony score of random tip shuffles) | — |
| **Parsimony p-value** (+ BH) | 0–1 | **< 0.05**: fewer state changes than chance | **> 0.05** |
| **Mantel r** | −1 to 1 | Positive (same-value tips closer on the tree); even 0.1–0.3 can matter | ≈ 0 or negative |
| **Mantel p-value** (+ BH) | 0–1 | **< 0.05** | **> 0.05** |

Purity, parsimony, and Mantel test the same idea in different ways; `Significant BH` is TRUE if *any* of the three is significant after BH. Of the three, the parsimony-score test is the best-established statistic (Maddison & Slatkin 1991) and, unlike purity, is not inflated by trivially-pure single-tip clades — weight it accordingly. Mantel is the weakest here (see Harmon & Glor 2010 below); treat it as corroborating.

### Continuous traits (`gene_trait_numeric_summary.csv`)
| Statistic | Range | Meaning / strong signal | No signal (not interesting) |
|---|---|---|---|
| **Blomberg's K** | 0 → ∞ | K≈1 = Brownian-motion expectation; **K>1** = close relatives *more* alike than Brownian (strong conservatism); 0<K<1 = weaker than Brownian | **K≈0** (trait ~ random on the tree) |
| **Blomberg K p-value** (+ BH) | 0–1 | **< 0.05**: signal exceeds random tip shuffles | **> 0.05** |
| **Pagel's λ** | 0–1 | **λ≈1** = signal as strong as Brownian motion; intermediate = partial | **λ≈0** (trait independent of tree) |
| **Pagel λ p-value** (+ BH) | 0–1 | **< 0.05** (likelihood-ratio test vs λ=0) | **> 0.05** |
| **Moran's I** | ≈ −1 to 1 | **Positive** = neighbouring tips have similar values (expected ≈ 0 under no autocorrelation) | ≈ 0 (random) or negative (over-dispersed) |
| **Moran's I p-value** (+ BH) | 0–1 | **< 0.05** | **> 0.05** |
| **Mantel r** | −1 to 1 | **Positive** = phylogenetically distant tips also differ more in the trait | ≈ 0 or negative |
| **Mantel p-value** (+ BH) | 0–1 | **< 0.05** | **> 0.05** |

The four continuous metrics capture different facets and **can legitimately disagree** — e.g. Blomberg's K can read ≈0 for a low-variance trait while Pagel's λ and Moran's I detect structure. Treat agreement across metrics (and significant BH q-values) as the strongest evidence; a single significant metric is weaker, especially if the effect size is small. `Trait_Mean`/`Trait_SD` are descriptive only — a very small SD means the trait is nearly constant and signal is hard to detect regardless of the tree.

### Method references
Papers that introduce or apply each test in the use-case it serves here. Reviews are the best starting point if you only read one or two.

**Choosing among the continuous signal metrics (review):**
- Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T., Schiffers, K., & Thuiller, W. (2012). How to measure and test phylogenetic signal. *Methods in Ecology and Evolution*, 3(4), 743–756. doi:10.1111/j.2041-210X.2012.00196.x — Compares Blomberg's K, Pagel's λ, and Moran's I, with practical guidance on statistical power and when each is appropriate.
- Revell, L.J., Harmon, L.J., & Collar, D.C. (2008). Phylogenetic signal, evolutionary process, and rate. *Systematic Biology*, 57(4), 591–601. doi:10.1080/10635150802302427 — How to interpret K and λ in terms of the underlying evolutionary process (i.e. why signal ≠ a specific model).

**Blomberg's K:**
- Blomberg, S.P., Garland, T. Jr., & Ives, A.R. (2003). Testing for phylogenetic signal in comparative data: behavioral traits are more labile. *Evolution*, 57(4), 717–745. doi:10.1111/j.0014-3820.2003.tb00285.x — Original definition of K and the K=1 Brownian-motion benchmark; uses tip randomization for significance, as done here.

**Pagel's λ:**
- Pagel, M. (1999). Inferring the historical patterns of biological evolution. *Nature*, 401(6756), 877–884. doi:10.1038/44766 — Introduces λ as a branch-scaling parameter.
- Freckleton, R.P., Harvey, P.H., & Pagel, M. (2002). Phylogenetic analysis and comparative data: a test and review of evidence. *The American Naturalist*, 160(6), 712–726. doi:10.1086/343873 — Establishes λ as a phylogenetic-signal statistic and the likelihood-ratio test against λ=0 used here.

**Moran's I (phylogenetic autocorrelation):**
- Moran, P.A.P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37(1/2), 17–23. doi:10.2307/2332142 — The autocorrelation statistic.
- Gittleman, J.L., & Kot, M. (1990). Adaptation: statistics and a null model for estimating phylogenetic effects. *Systematic Zoology*, 39(3), 227–241. doi:10.2307/2992183 — Applies Moran's I to phylogenies using distance-based weights, as here.

**Discrete-trait / binary association (tree-purity and parsimony-score permutation tests):**
- Maddison, W.P., & Slatkin, M. (1991). Null models for the number of evolutionary steps in a character on a phylogenetic tree. *Evolution*, 45(5), 1184–1197. doi:10.2307/2409726 — The parsimony-score permutation test used here: compare the observed minimum number of state changes (Fitch) to the distribution under randomly shuffled tips.
- Parker, J., Rambaut, A., & Pybus, O.G. (2008). Correlating viral phenotypes with phylogeny: accounting for phylogenetic uncertainty. *Infection, Genetics and Evolution*, 8(3), 239–246. doi:10.1016/j.meegid.2007.08.001 — Tip-randomization tests of whether a discrete trait clusters on a tree (association index, parsimony score, monophyletic-clade size); the BaTS tool. The clade-majority "purity" test used here is in the same family.
- Fritz, S.A., & Purvis, A. (2010). Selectivity in mammalian extinction risk and threat types: a new measure of phylogenetic signal strength in binary traits. *Conservation Biology*, 24(4), 1042–1051. doi:10.1111/j.1523-1739.2010.01455.x — The *D* statistic; a principled alternative worth knowing for binary traits.

**Mantel test (used for both trait types) — including its caveats:**
- Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research*, 27(2), 209–220. PMID:6018555 — Original matrix-correlation test.
- Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd English ed.), Developments in Environmental Modelling, Vol. 24. Elsevier, Amsterdam. ISBN 978-0-444-53868-0 — Standard treatment of Mantel tests on distance matrices.
- Harmon, L.J., & Glor, R.E. (2010). Poor statistical performance of the Mantel test in phylogenetic comparative analyses. *Evolution*, 64(7), 2173–2178. doi:10.1111/j.1558-5646.2010.00973.x — **Caution:** the paper shows the Mantel test can have low power and, in some cases, inflated type-I error on phylogenetic distances; treat the Mantel columns as corroborating, not primary, evidence.

**Implementation reference:**
- Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D., Blomberg, S.P., & Webb, C.O. (2010). Picante: R tools for integrating phylogenies and ecology. *Bioinformatics*, 26(11), 1463–1464. doi:10.1093/bioinformatics/btq166 — Reference implementations of K, Moran's I, and permutation tests that the calculations here follow.

## Instructions
### Required Software
Uses the [Snakemake](https://github.com/snakemake/snakemake) workflow manager (see [installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)). Per-rule dependencies are provided as conda environments under `bin/env/` (e.g. `ete3.yml` for the gene–trait rules).

### Required Files
Most required inputs are outputs of the [Plasmid Assembler and TF Annotation pipeline](https://github.com/Lagator-Group/Plasmid-Assembly-Characterisation-and-Annotation-Snakemake). Specifically:
- prokka: Output of [Prokka](https://github.com/tseemann/prokka) run on plasmid sequences.
- `metadata.csv`: Sample metadata table. The `gene_col` column (default `NUCCORE_ACC`) must match the leaf names in the treefiles.

### Config File
Open `config.yml` and adjust the parameters for the pipeline you are running. For `sfile_geneTrait`:
- `metadata`: path to the metadata CSV.
- `gene`: list of gene IDs to test (each needs `data/iqtree_protein/{gene}.treefile`).
- `gene_col`: metadata column holding the sample/accession name that matches the treefile leaves.
- `trait_col`: categorical trait columns (see above).
- `trait_num_col`: continuous trait columns (see above).

### Running the pipeline
`snakemake -s sfile_{pipe_name} -c8 --use-conda --conda-frontend conda` — adjust `-c#` to the available cores.

With `--use-conda`, the `sfile_*` pipelines build their conda env from `bin/env/*.yml` (portable, but the build is slow on a WSL `/mnt/c` filesystem). If the required env already exists locally, use the `_local` variant instead: `sfile_geneTrait_local` names the existing `ete3` env directly (`conda: ete3`), so `--use-conda` reuses it rather than rebuilding:

`snakemake -s sfile_geneTrait_local -c8 --use-conda --conda-frontend conda`

Typical order of operations:

`snakemake -s sfile_muscle -c8 --use-conda --conda-frontend conda`

Then check the alignments and trim as needed, saving the trimmed alignments in the `data/muscle_trimmed` directory.

`snakemake -s sfile_iqtree -c8 --use-conda --conda-frontend conda`

`snakemake -s sfile_geneTrait -c8 --use-conda --conda-frontend conda`
