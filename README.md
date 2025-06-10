# RareNet: Rare Variant Network-based Association Analysis
![Figure 1. RareNet workflow](rarenet_overview.pdf)

RareNet is an R package for powerful gene-level rare-variant association analysis that leverages high-confidence protein-protein interaction networks from [STRING](https://string-db.org/) to create biologically informed gene sets. It computes per-gene p-values using [SAIGE-GENE+](https://github.com/weizhouUMICH/SAIGE) on large sequencing cohorts, then applies [GAUSS](https://github.com/diptavo/GAUSS) (powered by a pre-computed, empirically derived reference panel of gene-gene correlations) for subset-based association testing of the PPI-derived sets. RareNet merges the SAIGE-GENE+ and GAUSS p-values to deliver a single, network-aware gene-level statistic with improved power and controlled type 1 error rates. 


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
library(RareNet)

# By default, RareNet uses the shipped geneset_string_v12.txt:
rareNet(
  phenoFile     = "path/to/phenotypes.txt",
  plinkPrefix   = "path/to/pruned_prefix",
  bedFile       = "path/to/genotypes.bed",
  bimFile       = "path/to/genotypes.bim",
  famFile       = "path/to/genotypes.fam",
  workDir       = "results/rareNet/",
  threads       = 4
)
closeAllConnections()
# → writes results/rareNet/rareNet_results.txt with columns Gene and p.value
```

## Custom Gene Set

By default, RareNet uses the shipped `geneset_string_v12.txt`.  
To override this, supply your own gene-set file via the `geneSetFile` argument.  

**Requirements for your custom file** (TSV with **three** columns and a header row):  
1. **GeneSet** — the name of the primary gene (e.g. `A1BG`)  
2. **DESC**    — a short description (often the same as GeneSet)  
3. **GeneSet,Genes** — a comma-separated list of **11** genes (the primary gene plus its top 10 neighbors)  

Example first two lines of `my_custom_geneset.txt`:
```txt
GeneSet DESC GeneSet,Genes
A1BG	A1BG	A1BG,FAM3C,ORM2,OSCAR,LGALS3BP,LAIR2,VWF,FCAR,CRISP3,ITIH4,VSTM1
```

**How to use**  
- **Argument:** `geneSetFile`  
- **Description:** Path to your custom TSV; if omitted, RareNet falls back to the built-in `geneset_string_v12.txt`.  

```r
rareNet(..., geneSetFile = "data/my_custom_geneset.txt")
```
