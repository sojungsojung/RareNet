# R/rareNet.R

#' RareNet: full SAIGE → GAUSS → weighted‐Cauchy pipeline
#'
#' @param phenoFile     Path to phenotype file (tab-delimited)
#' @param plinkPrefix   Prefix for pruned PLINK files (step1; no “.bed/.bim/.fam”)
#' @param bedFile       Path to .bed file (step2)
#' @param bimFile       Path to .bim file (step2)
#' @param famFile       Path to .fam file (step2)
#' @param geneSetFile   Path to gene-set GMT (default shipped geneset_string_v12.txt)
#' @param referenceFile Path to reference panel (default shipped reference_panel.txt)
#' @param workDir       Directory for intermediate outputs (default: tempdir())
#' @param threads       Number of threads for SAIGE (default: 4)
#' @return data.table with columns:
#'   - Gene  
#'   - saige.p  
#'   - gauss.p  
#'   - coreSubset  
#'   - p.value       (final RareNet p-value)  
#'   - source        ("combined" or "saige_only")
#' @export
rareNet <- function(phenoFile,
                    plinkPrefix,
                    bedFile,
                    bimFile,
                    famFile,
                    geneSetFile    = system.file("data","geneset_string_v12.txt", package="RareNet"),
                    referenceFile  = system.file("data","reference_panel.txt",  package="RareNet"),
                    workDir        = tempdir(),
                    threads        = 4) {
  library(data.table)

  # 1) Run SAIGE‐GENE+
  saigeDir <- file.path(workDir, "saige")
  saigeRes <- run_saige_gene(
    phenoFile     = phenoFile,
    plinkPrefix   = plinkPrefix,
    bedFile       = bedFile,
    bimFile       = bimFile,
    famFile       = famFile,
    groupFile     = geneSetFile,
    outputDir     = saigeDir,
    threads       = threads
  )
  setnames(saigeRes, c("Gene","saige.p"))

  # 2) Run GAUSS
  gaussDir <- file.path(workDir, "gauss")
  gaussRes <- run_gauss(
    saigeRes      = saigeRes[, .(Gene, p.value = saige.p)],
    geneSetFile   = geneSetFile,
    referenceFile = referenceFile,
    outputDir     = gaussDir
  )
  setnames(gaussRes, c("geneSet","gauss.p","coreSubset"))

  # 3) Merge results
  dt <- merge(
    saigeRes,
    gaussRes[, .(Gene = geneSet, gauss.p, coreSubset)],
    by = "Gene", all.x = TRUE
  )

  # 4) Compute final p-value (rename combined.p → p.value)
  dt[, combined.p := ifelse(
    Gene %in% unlist(coreSubset),
    cauchy_combine(saige.p, gauss.p),
    saige.p
  )]
  dt[, p.value := combined.p              ]  # new column
  dt[, source  := ifelse(
    Gene %in% unlist(coreSubset),
    "combined", "saige_only"
  )]

  # drop the temp combined.p
  dt[, combined.p := NULL]

  # 5) Write two-column output file
  out_file <- file.path(workDir, "rareNet_results.txt")
  fwrite(
    dt[, .(Gene, p.value)],
    file      = out_file,
    sep       = "\t",
    col.names = TRUE
  )

  # 6) Return full table
  return(dt)
}
