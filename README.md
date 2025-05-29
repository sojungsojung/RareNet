# RareNet

**RareNet** integrates [SAIGE-GENE+](https://github.com/weizhouUMICH/SAIGE) and [GAUSS](https://github.com/diptavo/GAUSS) to compute aggregate rare‐variant gene p-values and then combine them via a weighted Cauchy method.  
**It leverages protein–protein interaction networks from STRING-DB to group biologically related genes, boosting the power of gene-based tests.**

## Installation

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install from GitHub
devtools::install_github("sojungsojung/RareNet")

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
