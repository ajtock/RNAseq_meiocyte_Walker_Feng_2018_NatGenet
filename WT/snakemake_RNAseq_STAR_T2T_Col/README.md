# Automated workflow for RNA-seq data processing and alignment

This is a Snakemake workflow for automated processing and alignment of single-end RNA-seq data using [STAR](https://github.com/alexdobin/STAR).

### Requirements

- Installation of [Snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- Demultiplexed single-end RNA-seq reads in gzipped FASTQ format located in the `data/` directory. These should be named according to the following naming convention: `{sample}.fastq.gz`
- A SAMtools-indexed reference genome in FASTA format and a chromosome sizes file (e.g., `T2T_Col.fa`, `T2T_Col.fa.fai`, and `T2T_Col.fa.sizes`, the latter two of which generated with `samtools faidx T2T_Col.fa; cut -f1,2 T2T_Col.fa.fai > T2T_Col.fa.sizes`), each located in `data/index/`
- a STAR-indexed reference genome with annotation in GTF format, located in `data/index/`; `data/index/STAR_genome_index_T2T_Col.genes.repo.sh` is an example script for this purpose
- Lists of technical sequences to be supplied to FASTQC and Trimmomatic, provided in this repository: `adapters/contaminants_list_fastqc.txt` and `adapters/cat_all_and_TruSeq_Single_Indexes.fa`
- `Snakefile` in this repository. This contains "rules" that each execute a step in the workflow
- `config.yaml` in this repository. This contains customizable parameters including `reference`, which should be the relative path to the reference genome file name without the `.fa` extension (e.g., `data/index/T2T_Col`)
- The bash shell script `scripts/perbase_1based_coverage.sh` in this repository; this needs to be made executable first (e.g., `chmod +x scripts/perbase_1based_coverage.sh`)
- Optional: `environment.yaml` in this repository, used to create the software environment if conda is used
- If conda is not used, the tools listed in environment.yaml must be specified in the PATH variable

These files can be downloaded together by cloning the repository:

```
git clone https://github.com/ajtock/160601_Kyuha_RNAseq.git
```

### Creating the conda environment

```
conda env create --name RNAseq_mapping --file environment.yaml
```

Or, if using [mamba](https://github.com/mamba-org/mamba):

```
mamba env create --name RNAseq_mapping --file environment.yaml
```

### Usage

In a Unix shell, navigate to the base directory containing `Snakefile`, `config.yaml`, `environment.yaml`, and the `data/`, `scripts/` and `adapters/` subdirectories, which should have a directory tree structure like this:

```
.
├── adapters
│   ├── cat_all_and_TruSeq_Single_Indexes.fa
│   └── contaminants_list_fastqc.txt
├── config.yaml
├── data
│   ├── index
│   │   ├── chrLength.txt
│   │   ├── chrNameLength.txt
│   │   ├── chrName.txt
│   │   ├── chrStart.txt
│   │   ├── exonGeTrInfo.tab
│   │   ├── exonInfo.tab
│   │   ├── geneInfo.tab
│   │   ├── Genome
│   │   ├── genomeParameters.txt
│   │   ├── Log.out
│   │   ├── SA
│   │   ├── SAindex
│   │   ├── sjdbInfo.txt
│   │   ├── sjdbList.fromGTF.out.tab
│   │   ├── sjdbList.out.tab
│   │   ├── STAR_genome_index_T2T_Col.genes.repo.sh
│   │   ├── T2T_Col.fa
│   │   ├── T2T_Col.fa.fai
│   │   ├── T2T_Col.fa.sizes
│   │   ├── T2T_Col.genes.gtf
│   │   └── transcriptInfo.tab
│   ├── WT_leaf_RNAseq_Rep1.fastq.gz
│   └── WT_meiocyte_RNAseq_Rep1.fastq.gz
├── environment.yaml
├── README.md
├── scripts
    ├── genomeBin_bedgraphToTSV.R
│   └── perbase_1based_coverage.sh
└── Snakefile
```

Then run the following commands in the base directory (`--cores` should match the `THREADS` parameter in `config.yaml`):

```
conda activate RNAseq_mapping
snakemake -p --cores 48
conda deactivate
```

### Useful Snakemake parameters

- `--cores` specifies the maximum number of threads
- `-n` performs a dry run
- `-p` prints commands
- `--use-conda`
- `--conda-prefix ~/.myconda`
- `--forcerun genomeBin_bedgraphTPM` forces rerun of a given rule (e.g., `genomeBin_bedgraphTPM`)

### Updating the conda environment

```
conda env update --name RNAseq_mapping --file environment.yaml
```

Or, if using [mamba](https://github.com/mamba-org/mamba):

```
mamba env update --name RNAseq_mapping --file environment.yaml
```
