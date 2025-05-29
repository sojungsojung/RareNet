#' Run SAIGE-GENE+ on a phenotype
#'
#' @param phenoFile Path to phenotype file (tab-delimited)
#' @param genotypeFiles Character vector of genotype file paths (BGEN or PLINK)
#' @param groupFile Path to gene-set (group) file
#' @param outputDir Directory to store SAIGE output
#' @return data.table with columns Gene, p.value
#' @export
run_saige_gene <- function(phenoFile,
                           genotypeFiles,
                           groupFile,
                           outputDir) {
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  # build and run the SAIGE step2 command
  cmd <- sprintf(
    "docker run -v /data:/data wzhou88/saige:1.4.4 step2_SPAtests.R \
     --phenoFile=%s \
     --plinkFile=%s \
     --groupFile=%s \
     --outputPrefix=%s/saige_out \
     --nThreads=4 \
     --LOCO=FALSE",
    phenoFile,
    paste(genotypeFiles, collapse = ","),
    groupFile,
    outputDir
  )
  status <- system(cmd)
  if (status != 0) {
    stop("SAIGE-GENE+ failed with exit status ", status)
  }

  # SAIGE step2 produces <outputPrefix>.SAIGE_GENE_output.txt
  out_file <- file.path(outputDir, "saige_out.SAIGE_GENE_output.txt")
  if (!file.exists(out_file)) {
    stop("Cannot find SAIGE output at ", out_file)
  }

  # read and return
  res <- data.table::fread(out_file)
  data.table::setnames(res, c("Gene", "p.value"))
  return(res)
}
