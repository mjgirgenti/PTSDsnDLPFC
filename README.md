# Single cell transcriptomic and epigenomic atlas of >2M nuclei from the human PTSD brain

This is the source code repository for *Single cell transcriptomic and epigenomic atlas of >2M nuclei from the human PTSD brain*. The main analyses use both Python and R, with specific packages detailed below. We include tutorials and expected output (including runtime) to show how it runs on our dataset. 

## System Requirements

Both Python and R are required to replicate the analysis pipeline. We used:

- Python 3.9.16
- R 4.2.3

in our analyses.

Please use the attached `environment.yml` for the full list of system requirements and versions. 

## Installation

To download and install the pipeline, simply clone this repository:

```
git clone https://github.com/mjgirgenti/PTSDsnDLPFC.git
```

and create a virtual environment with `conda` using the supported `environment.yml`:

```
conda env create -f environment.yml
``` 

## Usage

The repository is divided into directories: `RNA`, `CCC`, `ATAC`, `Multiome`, and `GWAS`, representing major analyses of the paper. In each folder, we split analysis steps into separate files, with detailed comments. Please also refer to the jupyter notebook tutorials listed below for instructions. 

## Quick Tutorials

We have curated tutorials specific to each modality of the analysis to show how our pipeline works. Please review the README.md files at each modality for details.

In addition, we provided the [Jupyter notebooks](https://github.com/mjgirgenti/PTSDsnDLPFC/tree/main/notebooks) that we used to generate figures, with sample input and output. 

## Cite Our Work

To be added...
