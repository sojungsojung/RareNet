# RareNet: **Rare** Variant **Net**work-based Association Analysis

RareNet is a network based rare variant gene based association analysis tool with improved power that leverages high-confidence protein–protein interaction networks from STRING-DB. It integrates [SAIGE-GENE+](https://github.com/weizhouUMICH/SAIGE) and [GAUSS](https://github.com/diptavo/GAUSS) to compute aggregate rare‐variant gene p-values.

# RareNet: **Rare** Variant **Net**work-based **A**ssociation **A**nalysis  

RareNet is an R package for powerful gene‐level rare‐variant testing that:  
1. Leverages high‐confidence protein–protein interaction networks from STRING-DB to define biologically informed gene sets :contentReference[oaicite:0]{index=0}  
2. Runs [SAIGE-GENE+](https://github.com/weizhouUMICH/SAIGE) to compute per-gene p-values on large sequencing cohorts :contentReference[oaicite:1]{index=1}  
3. Applies [GAUSS](https://github.com/diptavo/GAUSS) for summary-statistics–based gene-set association on those PPI-derived sets :contentReference[oaicite:2]{index=2}  
4. Combines SAIGE and GAUSS p-values via a weighted Cauchy (ACAT) approach for a unified, network-aware gene test  



## Installation

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install from GitHub
devtools::install_github("sojungsojung/RareNet")
```

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
