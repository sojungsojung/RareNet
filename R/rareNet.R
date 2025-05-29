#' RareNet: full SAIGE → GAUSS → weighted‐Cauchy pipeline
#'
#' @description
#' 1) Runs SAIGE-GENE+ (step1 & step2) and concatenates per‐chromosome results  
#' 2) Runs GAUSS on that SAIGE summary + gene‐set + reference panel
#' 3) For genes within the GAUSS core subset, computes weighted Cauchy combine;  
#'    otherwise retains SAIGE p-value  
#'
#' @param phenoFile     Path to SAIGE phenotype file
#' @param plinkFile     Prefix for PLINK files (no extension)
#' @param bedFile       Path to .bed file
#' @param bimFile       Path to .bim file
#' @param famFile       Path to .fam file
#' @param geneSetFile   Path to gene‐set (GMT) file
#' @param referenceFile Path to reference null weight file
#' @param workDir       Directory for intermediate outputs (default: tempdir())
#' @param threads       Number of threads for SAIGE (default: 4)
#' @return data.table with columns:
#'   - Gene, saige.p, gauss.p, coreSubset, combined.p, source
#' @export
rareNet <- function(phenoFile,
                    plinkFile,
                    bedFile,
                    bimFile,
                    famFile,
                    geneSetFile,
                    referenceFile,
                    workDir = tempdir(),
                    threads = 4) {
  library(data.table)

  # 1) SAIGE summary
  saigeDir <- file.path(workDir, "saige")
  saigeRes <- run_saige_gene(
    phenoFile, plinkFile, bedFile, bimFile, famFile,
    geneSetFile, saigeDir, threads
  )
  setnames(saigeRes, c("Gene", "saige.p"))

  # 2) GAUSS
  gaussDir <- file.path(workDir, "gauss")
  gaussRes <- run_gauss(
    saigeRes     = saigeRes[, .(Gene, p.value = saige.p)],
    geneSetFile  = geneSetFile,
    referenceFile= referenceFile,
    outputDir    = gaussDir
  )
  setnames(gaussRes, c("geneSet", "p.value", "coreSubset"))

  # 3) Merge
  dt <- merge(
    saigeRes,
    gaussRes[, .(Gene = geneSet, gauss.p = p.value, coreSubset)],
    by = "Gene", all.x = TRUE
  )

  # 4) Combine p-values
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

  return(dt)
}
