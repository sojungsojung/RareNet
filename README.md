# RareNet: Rare Variant Network-based Association Analysis

RareNet is an R package for powerful gene-level rare-variant association analysis that leverages high-confidence protein-protein interaction networks from [STRING](https://string-db.org/) to create biologically informed gene sets. It computes per-gene p-values using [SAIGE-GENE+](https://github.com/weizhouUMICH/SAIGE) on large sequencing cohorts, then applies [GAUSS](https://github.com/diptavo/GAUSS)-powered by a pre-computed, empirically derived reference panel of gene-gene correlations- for subset-based association testing of the PPI-derived sets. RareNet merges the SAIGE-GENE+ and GAUSS p-values to deliver a single, network-award gene-level statistic with improved power and controlled type 1 error rates. 


## Installation

```r
# Install system prerequisites:
# • R ≥ 4.0 with the SAIGE and GAUSS packages installed
# • data.table, dplyr, foreach, doParallel

# Install devtools if needed
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")

# Install RareNet from GitHub
devtools::install_github("sojungsojung/RareNet")
```

## Prerequisites

- **Docker** (required to run the SAIGE-GENE+ containers)  
- **R** (≥ 4.0) with the following packages installed:
  - `SAIGE` (for SAIGE-GENE+)  
  - `GAUSS`  
  - `data.table`  
  - `dplyr`  
  - `foreach`  
  - `doParallel`  

## Usage

```r
## Usage

```r
## Usage

```r
# Load RareNet
library(RareNet)

# Run RareNet without assignment; this will write
# “rareNet_results.txt” under your workDir and close connections.
rareNet(
  phenoFile     = "path/to/phenotypes.txt",    # must include IID & Phenotype columns
  plinkPrefix   = "path/to/pruned_prefix",     # for SAIGE step1 (.bed/.bim/.fam)
  bedFile       = "path/to/genotypes.bed",     # for SAIGE step2
  bimFile       = "path/to/genotypes.bim",
  famFile       = "path/to/genotypes.fam",
  geneSetFile   = system.file("data","geneset_string_v12.txt",   package="RareNet"),
  referenceFile = system.file("data","reference_panel.txt",      package="RareNet"),
  workDir       = "results/rareNet/",          # where to save rareNet_results.txt
  threads       = 4
)
closeAllConnections()
# Output
#   results/rareNet/rareNet_results.txt   (tab-delimited: Gene, p.value)
```
