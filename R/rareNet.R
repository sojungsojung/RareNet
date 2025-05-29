# R/rareNet.R

#' RareNet: full SAIGE → GAUSS → weighted‐Cauchy pipeline (side‐effecting)
#'
#' @param phenoFile     Path to phenotype file (tab-delimited; must include IID & Phenotype)
#' @param plinkPrefix   Prefix for pruned PLINK files (SAIGE step1; no “.bed/.bim/.fam”)
#' @param bedFile       Path to .bed file (SAIGE step2)
#' @param bimFile       Path to .bim file (SAIGE step2)
#' @param famFile       Path to .fam file (SAIGE step2)
#' @param geneSetFile   Path to gene‐set file (default shipped geneset_string_v12.txt)
#' @param referenceFile Path to reference panel (default shipped reference_panel.txt)
#' @param workDir       Directory for intermediate outputs & results (default: tempdir())
#' @param threads       Number of threads for SAIGE (default: 4)
#' @export
# R/rareNet.R
rareNet <- function(phenoFile,
                    plinkPrefix,
                    bedFile,
                    bimFile,
                    famFile,
                    geneSetFile    = system.file("data", "geneset_string_v12.txt", package="RareNet"),
                    referenceFile  = system.file("data", "reference_panel.txt",   package="RareNet"),
                    workDir        = tempdir(),
                    threads        = 4) {
  library(data.table)

  # 1) SAIGE‐GENE+
  saigeRes <- run_saige_gene(
    phenoFile   = phenoFile,
    plinkPrefix = plinkPrefix,
    bedFile     = bedFile,
    bimFile     = bimFile,
    famFile     = famFile,
    groupFile   = geneSetFile,   # ← uses default or user‐supplied
    outputDir   = file.path(workDir, "saige"),
    threads     = threads
  )
  setnames(saigeRes, c("Gene","saige.p"))

  # 2) GAUSS
  gaussRes <- run_gauss(
    saigeRes      = saigeRes[, .(Gene, p.value = saige.p)],
    geneSetFile   = geneSetFile,   # ← same
    referenceFile = referenceFile,
    outputDir     = file.path(workDir, "gauss")
  )
  setnames(gaussRes, c("geneSet","gauss.p","coreSubset"))

  # 3) Merge, combine, and write
  dt <- merge(
    saigeRes,
    gaussRes[, .(Gene = geneSet, gauss.p, coreSubset)],
    by = "Gene", all.x = TRUE
  )
  dt[, combined := ifelse(Gene %in% unlist(coreSubset), cauchy_combine(saige.p, gauss.p), saige.p)]
  dt[, p.value := combined][, combined := NULL]
  dt[, source := fifelse(Gene %in% unlist(coreSubset), "combined","saige_only")]

  # Write the final two-column result
  out_file <- file.path(workDir, "rareNet_results.txt")
  fwrite(dt[, .(Gene, p.value)], out_file, sep = "\t", col.names = TRUE)

  # Clean up
  closeAllConnections()
}
