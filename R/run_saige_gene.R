#' Run SAIGE-GENE+ on a phenotype
#'
#' @param phenoFile Path to phenotype file
#' @param genotypeFiles Character vector of genotype file paths (BGEN/PLINK)
#' @param groupFile Path to group file
#' @param outputDir Directory to store SAIGE output
#' @return data.table with columns Gene, p.value
#' @export
run_saige_gene <- function(phenoFile, genotypeFiles, groupFile, outputDir) {
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  cmd <- sprintf(
    "Rscript %s/step2_SPAtests.R \
      --phenoFile=%s \
      --plinkFile=%s \
      --groupFile=%s \
      --outputPrefix=%s/saige_out \
      --nThreads=4 \
      --LOCO=FALSE",
    system.file("extdata", package="SAIGE"),
    phenoFile,
    paste(genotypeFiles, collapse=","),
    groupFile,
    outputDir
  )
  status <- system(cmd)
  if (status != 0) stop("SAIGE-GENE+ failed with status ", status)
  out_file <- file.path(outputDir, "saige_out.html")  # adjust to actual merged.out path
  res <- data.table::fread(out_file)
  setnames(res, c("Gene","p.value"))
  return(res)
}
