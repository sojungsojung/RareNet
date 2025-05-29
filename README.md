# RareNet: Rare Variant Network-based Association Analysis

RareNet is an R package for powerful gene-level rare-variant testing that:  
1. Defines gene sets from high-confidence STRING-DB protein–protein interaction networks 
2. Computes per-gene p-values via SAIGE-GENE+ for large-scale sequencing cohorts 
3. Leverages a pre-computed, empirically derived reference panel of gene–gene correlations 
4. Runs GAUSS on PPI-derived sets using that reference panel for subset–based association analysis
5. Integrates SAIGE-GENE+ and GAUSS p-values for a unified, network-aware gene-level test 


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
