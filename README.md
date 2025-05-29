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
# Load the RareNet package
library(RareNet)

# 1) Run the full RareNet pipeline:
res <- rareNet(
  phenoFile     = "path/to/phenotypes.txt",                   # your SAIGE phenotype file
  genotypeFiles = c("data/chr1.bgen", "data/chr2.bgen"),      # provide your genotype files
  geneSetFile   = system.file("extdata",                      # default STRING v12 gene set
                              "string_v12_geneset.txt",
                              package = "RareNet"),
  workDir       = "results/rareNet/",                         # directory for intermediate outputs
  p.thresh      = 2.5e-06                                      # significance threshold
)

# 2) View the top gene hits
head(res)

# 3) Export the final gene–p-value table
write.table(
  res[, .(Gene, combined.p)],
  file      = "rareNet_results.txt",
  row.names = FALSE,
  col.names = TRUE,
  quote     = FALSE
)
```
