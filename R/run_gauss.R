#' Run GAUSS on SAIGE summary statistics
#'
#' @param saigeRes   data.table or data.frame with columns Gene and p.value
#' @param geneSetFile Path to GMT gene-set file (e.g., geneset_string_v12.txt)
#' @param referenceFile Path to reference null weights file (reference.txt)
#' @param outputDir   Directory to write GAUSS output files
#' @return data.table with columns geneSet, p.value, coreSubset (list column)
#' @export
run_gauss <- function(saigeRes,
                      geneSetFile,
                      referenceFile,
                      outputDir) {
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  # write SAIGE summary to disk in the format GAUSS expects
  summary_file <- file.path(outputDir, "saige_summary.txt")
  data.table::fwrite(
    saigeRes[, .(Gene, p.value)],
    file = summary_file,
    sep = "\t",
    col.names = FALSE
  )

  # load reference null weights
  ref_df <- data.table::fread(referenceFile, header = TRUE)

  # call GAUSS_All on all genes (start=1, stop=n)
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

  if (!file.exists(gauss_out)) {
    stop("GAUSS did not produce output at ", gauss_out)
  }

  # read and return
  res <- data.table::fread(gauss_out)
  data.table::setnames(res, c("geneSet", "p.value", "coreSubset"))
  return(res)
}
