#' RareNet: full SAIGE→GAUSS→Cauchy pipeline
#'
#' @param phenoFile Path to phenotype file
#' @param genotypeFiles Character vector of genotype files
#' @param geneSetFile Path to STRING‐DB gene‐set file
#' @param p.thresh Numeric, significance threshold
#' @param workDir Character, directory for intermediate files
#' @return data.table of Gene, saige.p, gauss.p, combined.p, source, coreSubset
#' @export
rareNet <- function(phenoFile, genotypeFiles, geneSetFile,
                    p.thresh = 2.5e-06,
                    workDir  = tempdir()) {
  library(data.table)
  
  saigeDir <- file.path(workDir, "saige"); gaussDir <- file.path(workDir, "gauss")
  saigeRes <- run_saige_gene(phenoFile, genotypeFiles, geneSetFile, saigeDir)
  gaussRes <- run_gauss(saigeRes, geneSetFile, gaussDir)
  
  dt <- merge(
    saigeRes, 
    gaussRes[, .(Gene=geneSet, gauss.p=p.value, coreSubset)], 
    by = "Gene", all.x = TRUE
  )
  dt[, combined.p := ifelse(Gene %in% unlist(coreSubset),
                            cauchy_combine(p.value, gauss.p),
                            p.value)]
  dt[, source := ifelse(Gene %in% unlist(coreSubset),
                        "combined","saige_only")]
  return(dt[combined.p < p.thresh])
}
