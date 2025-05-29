#' Run GAUSS on SAIGE summary stats
#'
#' @param saigeRes data.table with Gene, p.value
#' @param geneSetFile Path to GMT or gene‐set file
#' @param outputDir Directory for GAUSS results
#' @return data.table with columns geneSet, p.value, coreSubset
#' @export
run_gauss <- function(saigeRes, geneSetFile, outputDir) {
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  res <- GAUSS::gauss_main(
    summary_stats = saigeRes,
    geneset_file   = geneSetFile,
    out_dir        = outputDir
  )
  return(res)
}
