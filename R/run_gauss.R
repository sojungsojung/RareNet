# R/run_gauss.R

#' Run GAUSS on SAIGE summary statistics
#'
#' @param saigeRes      data.table with columns Gene and p.value
#' @param geneSetFile   Path to gene-set GMT (default shipped geneset_string_v12.txt)
#' @param referenceFile Path to reference panel (default shipped reference_panel.txt)
#' @param outputDir     Directory to write GAUSS outputs
#' @return data.table with columns geneSet, p.value, coreSubset
#' @export
run_gauss <- function(saigeRes,
                      geneSetFile   = system.file("data","geneset_string_v12.txt", package="RareNet"),
                      referenceFile = system.file("data","reference_panel.txt",  package="RareNet"),
                      outputDir) {
  library(GAUSS)
  library(data.table)
  library(dplyr)
  library(foreach)
  library(doParallel)

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  # write summary
  summary_file <- file.path(outputDir, "saige_summary.txt")
  fwrite(saigeRes[, .(Gene, p.value)], summary_file,
         sep = "\t", col.names = FALSE)

  # load reference panel
  ref_df <- fread(referenceFile, header = TRUE)

  # run GAUSS_All
  gauss_out <- file.path(outputDir, "gauss_output.txt")
  GAUSS::GAUSS_All(
    summary_file = summary_file,
    gene_name    = 1,
    pv_name      = 2,
    output_file  = gauss_out,
    gmt          = geneSetFile,
    ags          = "def",
    verbose      = FALSE,
    parallel     = FALSE,
    start        = 1,
    stop         = nrow(saigeRes),
    is.appx      = TRUE,
    pv.null.wt1  = ref_df
  )
  stopifnot(file.exists(gauss_out))
  res <- fread(gauss_out)
  setnames(res, c("geneSet","p.value","coreSubset"))
  return(res)
}
