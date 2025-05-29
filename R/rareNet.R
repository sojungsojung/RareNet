#' RareNet: full SAIGE → GAUSS → weighted‐Cauchy pipeline
#'
#' @param phenoFile     Path to phenotype file (tab-delimited)
#' @param genotypeFiles Character vector of genotype prefixes (PLINK) or BGEN
#' @param geneSetFile   Path to gene-set GMT file (default shipped geneset_string_v12.txt)
#' @param referenceFile Path to reference panel (default shipped reference_panel.txt)
#' @param workDir       Directory for intermediate outputs (default: tempdir())
#' @param threads       Number of threads for SAIGE (default: 4)
#' @return data.table with Gene, saige.p, gauss.p, coreSubset, combined.p, source
#' @export
rareNet <- function(phenoFile,
                    genotypeFiles,
                    geneSetFile    = system.file("data","geneset_string_v12.txt", package="RareNet"),
                    referenceFile  = system.file("data","reference_panel.txt",  package="RareNet"),
                    workDir        = tempdir(),
                    threads        = 4) {
  library(data.table)

  # 1) Run SAIGE-GENE+
  saigeDir <- file.path(workDir, "saige")
  saigeRes <- run_saige_gene(
    phenoFile     = phenoFile,
    genotypeFiles = genotypeFiles,
    groupFile     = geneSetFile,
    outputDir     = saigeDir,
    threads       = threads
  )
  setnames(saigeRes, c("Gene","saige.p"))

  # 2) Run GAUSS
  gaussDir <- file.path(workDir, "gauss")
  gaussRes <- run_gauss(
    saigeRes      = saigeRes[, .(Gene,p.value=saige.p)],
    geneSetFile   = geneSetFile,
    referenceFile = referenceFile,
    outputDir     = gaussDir
  )
  setnames(gaussRes, c("geneSet","gauss.p","coreSubset"))

  # 3) Merge and combine
  dt <- merge(
    saigeRes,
    gaussRes[, .(Gene=geneSet, gauss.p, coreSubset)],
    by="Gene", all.x=TRUE
  )
  dt[, combined.p := ifelse(
    Gene %in% unlist(coreSubset),
    cauchy_combine(saige.p, gauss.p),
    saige.p
  )]
  dt[, source := ifelse(
    Gene %in% unlist(coreSubset),
    "combined","saige_only"
  )]

  return(dt)
}
