# Plasmid Transcription Factor Phylogeny
The available snakemake pipelines are designed to align sequences with Muscle, build trees with IQ Tree, determine relationships between 2 trees, and determine relationships between 1 tree and a trait.

Briefly, this pipeline:
sfile_muscle: Aligns sequences already in a single file using [Muscle](https://github.com/rcedgar/muscle)
sfile_iqtree: Constructs gene phylogenetic tree with [IQ-Tree](https://github.com/Cibiv/IQ-TREE).
sfile_cophenetic: Compares phylogeny of two genes using their `.treefile`.
sfile_geneTrait: Determines if there is correlation between a gene tree structure and metadata trait.

## Instructions
### Required Software
Uses [Snakemake](https://github.com/snakemake/snakemake) pipeline for sequence alignment and annotation. Needs Snakemake environment to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
### Required Files
Most of the files required are outputs of [Plasmid Assembler and TF Annotation pipeline](https://github.com/Lagator-Group/Plasmid-Assembly-TF-Annotation-Snakemake). Specifically, these items are:
- prokka: Output of [Prokka](https://github.com/tseemann/prokka) when run on plasmid sequences.
### Config File
Open `config.yml` and adjust the necessary parameters. 
### Running the pipeline
`snakemake -s sfile_{pipe_name} -c8 --use-conda --conda-frontend conda` Adjust `-c#` depending on available cores.

Then make sure to check the alignments and trim as needed. Save the trimmed alignments in `data/muscle_trimmed` directory.

`snakemake -s snakefile_iqtree -c8 --use-conda --conda-frontend conda`
