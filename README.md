# Allele-Resolved Hybrid Embryos Reveal the Fates of Regulatory Divergence

This repository contains the data-processing, analysis, and visualization code associated with the manuscript:

**"Allele-Resolved Hybrid Embryos Reveal the Fates of Regulatory Divergence"**

This study uses reciprocal interspecific hybrids between *Ciona intestinalis* and *Ciona savignyi* to investigate allele-specific gene expression and regulatory divergence during embryonic development at single-cell resolution.

## Repository structure

### `data_processing/`

Scripts used to process sequencing data and construct the inputs for downstream
allele-specific and integrative analyses.

- `01.allele.count.construction.r`  
  Construction of allelic count matrices.

- `02.allelic.rates.estimation.py`  
  Estimation of allelic expression rates from allelic counts.

- `03.scATAC.data.processing.r`  
  Processing of single-cell ATAC-seq data for downstream analyses.

- `04.multiome.data.processing.r`  
  Processing of single-cell multiome data.

- `05.peak.counts.construction.r`  
  Construction of allele-specific chromatin accessibility/peak count matrices.

### `figures/`

R scripts containing the analyses and visualization code used to generate the main and Extended Data figures in the manuscript.

- `Figure1.r` – Figure 1
- `Figure2.r` – Figure 2
- `Figure3.r` – Figure 3
- `Figure4.r` – Figure 4
- `Extended.Data.Figures.r` – Extended Data Figures

## Data availability

The single-cell RNA-seq, single-cell ATAC-seq, and single-cell multiome datasets generated in this study have been deposited in the NCBI Gene Expression Omnibus
(GEO) under accession **GSE346238**.

Processed data required for the analyses are also available through GEO.

## Software and dependencies

The analyses were performed primarily in R and Python.

Major R packages include:

- Seurat
- Signac
- ggplot2
- dplyr
- tidyverse
- clusterProfiler
- variancePartition

Allelic rate estimation was performed using Python. Additional package dependencies are specified in the individual scripts where applicable.

## Usage

The scripts in `data_processing/` generate allele-resolved expression and chromatin accessibility data used for downstream analyses.

The scripts in `figures/` contain the analysis and plotting workflows used to generate the figures presented in the manuscript.

File paths may need to be modified according to the user's local computing environment. Input and processed datasets are available through the GEO accession listed above.

## Citation

If you use the code or data from this repository, please cite:

Lei J, et al. Allele-Resolved Hybrid Embryos Reveal the Fates of Regulatory Divergence.
DOI: https://doi.org/xxxxx (To be updated)

## Contact

For questions regarding the code or data, please contact:

[Jiali Lei]  
[jiali.lei@utdallas.edu]  
