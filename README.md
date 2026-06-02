# **To be Updated**
# Plasmid Transcription Factor Phylogeny
This snakemake pipeline will generate the files necessary to build a phylogenetic tree for genes located on plasmids, and classify plasmids according to their incompatibility markers and pangenomes. Designed to be compatible with [Phandango](http://jameshadfield.github.io/phandango/#/main).
Briefly, this pipeline:
1. Merges amino acid (AA) sequences of gene of interest from all plasmids and output to a single file.
2. Aligns AA sequences with [Muscle](https://github.com/rcedgar/muscle).
3. Constructs gene phylogenetic tree with [IQ-Tree](https://github.com/Cibiv/IQ-TREE).
## Instructions
### Required Software
Uses [Snakemake](https://github.com/snakemake/snakemake) pipeline for sequence alignment and annotation. Needs Snakemake environment to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
### Required Files
Most of the files required are outputs of [Plasmid Assembler and TF Annotation pipeline](https://github.com/Lagator-Group/Plasmid-Assembly-TF-Annotation-Snakemake). Specifically, these items are:
- fasta_plasmid: Complete plasmid sequences in `.fasta` format.
- prokka: Output of [Prokka](https://github.com/tseemann/prokka) when run on plasmid sequences.
- sprot: Output of [Plasmid Assembler and TF Annotation pipeline](https://github.com/Lagator-Group/Plasmid-Assembly-TF-Annotation-Snakemake) which contains all Swissprot annotations.
### Config File
Open `config.yml` and adjust the necessary parameters. 
### Running the pipeline
`snakemake -s sfile_{pipe_name} -c8 --use-conda --conda-frontend conda` Adjust `-c#` depending on available cores.

Then make sure to check the alignments and trim as needed. Save the trimmed alignments in `data/muscle_trimmed` directory.

`snakemake -s snakefile_iqtree -c8 --use-conda --conda-frontend conda`
