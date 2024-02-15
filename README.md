# ONT-cappable-seq2

This document provides a comprehensive overview of the updated ONT-cappable-seq data analysis pipeline, now fully automated. For detailed information on individual steps, please refer to the [old ONT-cappable-seq repository](https://github.com/LoGT-KULeuven/ONT-cappable-seq).

## Installation

**Note:** This snakemake pipeline is designed for Linux systems! Windows users are advised to install [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install) before proceeding with the installation.

1. Install Snakemake by following the [instructions provided here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

2. Clone the ONT-cappable-seq2 repository to the desired location to obtain the Snakemake pipeline:

```bash
git clone https://github.com/LoGT-KULeuven/ONT-cappable-seq2.git
```

## Usage

### I. Data Navigation

Obtain high-quality FASTA files for the reference genomes of the organisms you are studying and store them in the `input/genome_data` directory. For better organization, it is advisable to rename the FASTA files using the respective IDs of the phage or host, referred to as `sampleID`. Additionally, download the FLuc Control Plasmid sequence used for in vitro transcription of the RNA control spike-in from the manufacturer's website. Save this sequence in an appropriate location.

```
ONT-cappable-seq2/
…
├── input/
|     └── genome_data/
|          ├── phageID.fasta
|          ├── hostID.fasta
|          └── in-vitro-RNA.fasta
```


Place the raw sequencing files in the `input/fastq_data` directory and rename them to incorporate information about `sampleID` and treatment (enriched or control sample). Ensure that each individual sample has a single decompressed .fastq file with an adequate number of reads and a suitable mean read length (> 200 bp). This naming convention and file organization will contribute to the clarity and efficiency of the analysis.

```
ONT-cappable-seq2/
…
├── input/
|     ├─ fastq_data/
|        ├── sampleID_enriched.fastq
|        └── sampleID_control.fastq
|     └── genome_data/
|        ├── phageID.fasta
|        ├── hostID.fasta
|        └── in-vitro-RNA.fasta
```
### II. Setting up the Configuration File

1. Navigate to the `config` directory within the ONT-cappable-seq project folder.

2. Open the `config.yaml` file. This file contains paths to your input files and various parameters used for annotating transcriptional boundaries. A general example of the config.yaml file is displayed below:

```yaml
Sample name: sampleID
fasta file: input/genome_data/sampleID.fasta
enriched fastq: input/fastq_data/sampleID_enriched.fastq
control fastq: input/fastq_data/sampleID_control.fastq
ID: group
termseq alpha: 0.001

cluster width:
  TSS: 15
  TTS: 30

minimum coverage:
  enriched:
    TSS: 25
    TTS: 25
  control:
    TSS: 2
    TTS: 2

peak alignment error: 2
TSS Threshold: 1.5
TTS threshold: 0.25

TSS sequence extraction:
  upstream: 40
  downstream: 0

TTS sequence extraction:
  upstream: 30
  downstream: 30
```

- **Sample name**: Contains the name of your organism of interest (e.g., sampleID, LUZ19).
- **Fasta file**: Path to your reference genome (e.g., `input/genome_data/sampleID.fasta`, `input/genome_data/LUZ19.fasta`).
- **Enriched fastq**: Path to the raw sequencing file of the enriched sample (e.g., `input/fastq_data/sampleID_enriched.fastq`, `input/fastq_data/LUZ19_enriched.fastq`).
- **Control fastq**: Path to the raw sequencing file of the control sample (e.g., `input/fastq_data/sampleID_control.fastq`, `input/fastq_data/LUZ19_control.fastq`).
- **ID**: Additional identifier for sample classification (e.g., timepoint, condition, replicate).
- **Termseq alpha**: Peak calling threshold for the termseq-peak algorithm.
- **Cluster width**: Distance within which peak positions are clustered.
- **Minimum coverage**: Minimum number of reads required for a candidate TSS/TTS position.
- **Peak alignment error**: Positional difference allowed between peak positions in enriched and control datasets.
- **TSS/TTS threshold**: Enrichment ratio and read reduction thresholds for TSS and TTS annotation.
- **TSS/TTS sequence extraction**: Upstream and downstream nucleotides relative to TSS and TTS for DNA sequence extraction.



### III. Running the Pipeline

Before executing the Snakemake pipeline, perform a dry-run to ensure the workflow is correctly defined:

```bash
cd ONT-cappable-seq2
snakemake --dry-run
```

To run the Snakemake pipeline, use the following commands:
```bash
cd ONT-cappable-seq2
snakemake --use-conda --cores 8
```
During the initial run, the pipeline will install all necessary Conda packages, which may take some time. Subsequent runs will reuse this Conda environment.

