#' RareNet: full SAIGE→GAUSS→Cauchy pipeline
#'
#' @description
#' Runs SAIGE-GENE+, then GAUSS, and for genes in the GAUSS core subset
#' combines their SAIGE and GAUSS p-values via a weighted Cauchy approach.
#' For all other genes, it retains the original SAIGE p-value.
#'
#' @param phenoFile Path to SAIGE phenotype file (tab-delimited).
#' @param genotypeFiles Character vector of genotype file paths (BGEN or PLINK).
#' @param geneSetFile Path to gene-set file (e.g., STRING v12 GMT).
#' @param workDir Directory for intermediate outputs (default: tempdir()).
#' @return data.table with columns:
#'   - Gene: gene symbol
#'   - saige.p: SAIGE-GENE+ p-value
#'   - gauss.p: GAUSS p-value (NA if gene not tested)
#'   - coreSubset: list column of genes in each core subset
#'   - combined.p: final p-value (Cauchy combined for core genes, else saige.p)
#'   - source: "combined" or "saige_only"
#' @export
rareNet <- function(
  phenoFile,
  genotypeFiles,
  geneSetFile,
  workDir = tempdir()
) {
  library(data.table)

  # 1) Run SAIGE-GENE+
  saigeDir <- file.path(workDir, "saige")
  dir.create(saigeDir, recursive = TRUE, showWarnings = FALSE)
  saigeRes <- run_saige_gene(
    phenoFile     = phenoFile,
    genotypeFiles = genotypeFiles,
    groupFile     = geneSetFile,
    outputDir     = saigeDir
  )
  setnames(saigeRes, c("Gene", "saige.p"))

  # 2) Run GAUSS
  gaussDir <- file.path(workDir, "gauss")
  dir.create(gaussDir, recursive = TRUE, showWarnings = FALSE)
  gaussRes <- run_gauss(
    saigeRes     = saigeRes[, .(Gene, p.value = saige.p)],
    geneSetFile  = geneSetFile,
    outputDir    = gaussDir
  )
  setnames(gaussRes, c("geneSet", "gauss.p", "coreSubset"))

  # 3) Merge results
  dt <- merge(
    saigeRes,
    gaussRes[, .(Gene = geneSet, gauss.p, coreSubset)],
    by = "Gene",
    all.x = TRUE
  )

  # 4) Compute combined p-value: only for genes in core subsets
  dt[, combined.p := ifelse(
    Gene %in% unlist(coreSubset),
    cauchy_combine(saige.p, gauss.p),
    saige.p
  )]
  dt[, source := ifelse(
    Gene %in% unlist(coreSubset),
    "combined",
    "saige_only"
  )]

  # 5) Return full table of all genes
  return(dt)
}
