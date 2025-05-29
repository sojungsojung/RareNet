# RareNet

**RareNet** integrates [SAIGE-GENE+](https://github.com/weizhouUMICH/SAIGE) and [GAUSS](https://github.com/diptavo/GAUSS) to compute aggregate rare‐variant gene p-values and then combine them via a weighted Cauchy method.  
**It leverages protein–protein interaction networks from STRING-DB to group biologically related genes, boosting the power of gene-based tests.**

# RareNet
**R**are **V**ariant **N**etwork-based **A**ssociation **A**nalysis

RareNet leverages high-confidence protein–protein interaction networks from STRING-DB to define biologically informed gene sets :contentReference[oaicite:0]{index=0}. It pipelines SAIGE-GENE+—a scalable mixed-model, set-based test optimized for large sequencing cohorts—to compute gene-level p-values :contentReference[oaicite:1]{index=1}. Next, it runs the GAUSS R package for summary-statistics–based gene-set association on these PPI-derived sets :contentReference[oaicite:2]{index=2}. Finally, RareNet combines SAIGE-GENE+ and GAUSS p-values via a weighted Cauchy (ACAT) approach for a unified, powerful gene-level test :contentReference[oaicite:3]{index=3}.



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
