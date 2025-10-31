# RSVrecon - RSV Genome Reconstruction Pipeline

![Bioinformatics Pipeline](https://img.shields.io/badge/bioinformatics-pipeline-blue)
![Python](https://img.shields.io/badge/python-3.10-green)
![R](https://img.shields.io/badge/R-4.3-green)
![License](https://img.shields.io/badge/license-MIT-orange)

Please visit our [nextflow implementation](https://github.com/stjudecab/rsvrecon) if you're familar with ![image](https://www.nextflow.io/img/nextflow.svg).

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Output](#output)
- [Dependencies](#dependencies)
- [Troubleshooting](#troubleshooting)

## Features
- Parallel processing with configurable thread usage
- Supports BWA alignment
- Generates consensus sequences
- Creates phylogenetic trees for RSV-A/RSV-B subtypes
- Produces PDF and HTML reports
- Quality control metrics including coverage statistics

## Installation

### 1. Clone Repository
#### Option A: Clone with Git (recommended)
```bash
git clone https://github.com/yourusername/RSVrecon.git
cd RSVrecon
```
#### Option B: Download ZIP
Download the latest release from Github
Unzip the package:
```bash
gunzip RSVrecon-main.zip
cd RSVrecon-main
```
#### 1.1: Download reference database
Download the pre-built reference database and unzip it to a location that you have read/write permission
https://github.com/stjudecab/RSVreconPy/releases/download/Pre-release/Reference.zip

### 2. Set Up Environment
We use `conda` to manage all dependencies. Please install `conda` and 'mamba' 
#### A1. Install Conda/Mamba (If you are not on a HPC)
**Install Miniconda** 

Please check conda website for a comprehansive instruction: https://www.anaconda.com/docs/getting-started/miniconda/install

**Install Mamba (recommended, it's much faster than conda)** 

Installing mamba

Once conda is installed, installing mamba with conda:
```bash
conda install mamba -c conda-forge
```

#### A2. Load module Conda/Mamba (If you are on a HPC)
Most high-performance computing (HPC) systems come with Conda/Mamba preinstalled. To use them:
Using our system as an example (please contact your HPC mamager for more details):
```bash
module load conda
module load mamba
```

#### B. Setup Env for RSVrecon
```bash
bash Set_env.sh
```

## Configuration
Example `config.yaml`:

```yaml
# Required paths
DATA_DIR: /path/to/input/fastq_files         # Please put all your paired-FASTQ files under this input folder
REFERENCE_DIR: /path/to/reference/sequences  # Please download our pre-built reference, unzip it, then paste the path here. Make sure you have both read and write permission
OUTPUT_DIR: /path/to/output/directory        # please specify a output folder path

# Performance parameters
THREAD_N: 2                     # Threads per sample, for BWA-MEM
MAX_CONCURRENT_JOBS: 10         # Parallel samples to process, notice: THREAD_N * MAX_CONCURRENT_JOBS should < than your number of CPUs

# Analysis parameters
TOOL: BWA                       # Currently only BWA supported
COV_CUTOFF: 50                  # Coverage cutoff threshold

# Optional
RSV_NEXT_PIPE_RES: /path/to/additional/results  # We allow users to compare RSVrecon with RSV-NEXT-PIPE results. Please specify the "consensus" folder of RSV-NEXT-PIPE output for the same batch of data.
```

## Quick Start
### 1. Download test dataset and prebuilt reference 
Download test dataset from [here](https://github.com/stjudecab/test_datasets/tree/rsvrecon). FastQ files are under "fastqs" folder.
A larger dataset is available [here](https://github.com/stjudecab/RSVreconPy/releases/download/Release-V0.2/Data.zip)

Download the pre-built reference database from [here](https://github.com/stjudecab/RSVreconPy/releases/download/Release-V0.2/Reference.zip)
### 2. Edit `config.yaml` with your paths
Here is an example:
```yaml
# Required paths
DATA_DIR: /path/to/input/fastq_files         # Please put all your paired-FASTQ files under this input folder
REFERENCE_DIR: /path/to/reference/sequences  # Please download our pre-built reference, unzip it, then paste the path here. Make sure you have both read and write permission
OUTPUT_DIR: /path/to/output/directory        # please specify a output folder path

# Performance parameters
THREAD_N: 2                     # Threads per sample, for BWA-MEM
MAX_CONCURRENT_JOBS: 10         # Parallel samples to process, notice: THREAD_N * MAX_CONCURRENT_JOBS should < than your number of CPUs

# Analysis parameters
TOOL: BWA                       # Currently only BWA supported
COV_CUTOFF: 50                  # upper threshold in the dual-coverage cutoff system 
COV_CUTOFF_LOW: 10              # lower threshold in the dual-coverage cutoff system

# Optional
RSV_NEXT_PIPE_RES: /path/to/additional/results  # We allow users to compare RSVrecon with RSV-NEXT-PIPE results. Please specify the "consensus" folder of RSV-NEXT-PIPE output for the same batch of data. You can disable it using "#"
```
### 3. Run pipeline:

```bash
# export path to your PATH
export PATH=/path/to/your/RSVrecon/folder:$PATH
# activate conda env
conda activate RSVreconEnv
```

```bash
# if you're on your local server
python rsvrecon_pipeline.py config.yaml

# If you're on HPC (using LSF as example)
# number of CPUs requested should >= THREAD_N * MAX_CONCURRENT_JOBS
bsub -n 20 -R "rusage[mem=10001]" -P CAB -J RSV -q priority -cwd $(pwd -P) "python rsvrecon_pipeline.py config.yaml"
```

## Output
```
Report/
├── Mapping/          # Alignment results
├── log/              # Log files
├── Temp/             # Temporary files
├── Report.csv        # Summary table
├── Sequence_*.fasta  # Consensus sequences
├── Report.pdf        # PDF report
└── Report.html       # HTML report
```

## Dependencies
Managed via `RSV_env.yml`:

```yaml
dependencies:
  # R related
  - r-base=4.3
  - r-ggplot2
  - r-biocmanager
  - bioconductor-ggtree=3.10.0
  - bioconductor-treeio
  - r-tidyverse
  - r-devtools
  
  # Python related
  - python=3.10
  - pandas=2.2.2
  - biopython=1.78
  - pyhocon
  - reportlab
  - matplotlib
  - seaborn
  - Pillow
  - pyyaml

  # Bioinformatics tools
  - bioconda::fastp=0.23.4
  - bioconda::igvtools=2.3.93
  - bioconda::kma=1.4.9
  - bioconda::nextclade
  - bioconda::samtools=1.18
  - bioconda::blast=2.14.1
  - bioconda::bwa
  - bioconda::mafft=7.505
  - bioconda::iqtree
```

## Troubleshooting
### Common Issues:

- Environment creation fails → Try `conda env create -f RSV_env.yml`
- Pipeline errors → Check `log/*.err.log` files
- Memory issues → Reduce `MAX_CONCURRENT_JOBS`

## Citation
Our preprint is on-line at [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.03.657184v1)
