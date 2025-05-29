# R/run_saige_gene.R

#' Run SAIGE-GENE+ on a phenotype (step1 → step2 → concatenate)
#'
#' @param phenoFile     Path to phenotype file (must include IID & Phenotype columns)
#' @param plinkPrefix   Prefix for pruned PLINK files (no “.bed/.bim/.fam”)
#' @param bedFile       Path to the genotype .bed file for step2
#' @param bimFile       Path to the genotype .bim file for step2
#' @param famFile       Path to the genotype .fam file for step2
#' @param groupFile     Path to gene-set file (default shipped geneset_string_v12.txt)
#' @param outputDir     Directory to hold SAIGE outputs
#' @param threads       Number of threads for SAIGE (default 4)
#' @return data.table with columns Gene, p.value
#' @export
run_saige_gene <- function(phenoFile,
                           plinkPrefix,
                           bedFile,
                           bimFile,
                           famFile,
                           groupFile    = system.file("data","geneset_string_v12.txt", package="RareNet"),
                           outputDir,
                           threads      = 4) {
  library(data.table)

  # ensure inputs exist
  stopifnot(file.exists(phenoFile),
            file.exists(paste0(plinkPrefix, ".bed")),
            file.exists(bedFile), file.exists(bimFile), file.exists(famFile),
            file.exists(groupFile))

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  # 1) Fit null GLMM (step1)
  step1_pref <- file.path(outputDir, "step1")
  cmd1 <- sprintf(
    "Rscript %s/step1_fitNULLGLMM.R \
     --plinkFile=%s \
     --phenoFile=%s \
     --phenoCol=Phenotype \
     --covarColList=Sex \
     --sampleIDColinphenoFile=IID \
     --traitType=binary \
     --invNormalize=TRUE \
     --outputPrefix=%s \
     --nThreads=%d \
     --LOCO=FALSE \
     --IsOverwriteVarianceRatioFile=TRUE",
    system.file("extdata", package="SAIGE"),
    plinkPrefix,
    phenoFile,
    step1_pref,
    threads
  )
  message("→ Running SAIGE step1 (null model)…")
  if (system(cmd1) != 0) stop("SAIGE step1 failed")

  # 2) SPAtests per chromosome (step2)
  cmd2 <- sprintf(
    "Rscript %s/step2_SPAtests.R \
     --bedFile=%s \
     --bimFile=%s \
     --famFile=%s \
     --varianceRatioFile=%s.varianceRatio.txt \
     --GMMATmodelFile=%s.rda \
     --groupFile=%s \
     --outputPrefix=%s/step2_out \
     --nThreads=%d \
     --LOCO=FALSE",
    system.file("extdata", package="SAIGE"),
    bedFile,
    bimFile,
    famFile,
    step1_pref,
    step1_pref,
    groupFile,
    outputDir,
    threads
  )
  message("→ Running SAIGE step2 (SPAtests)…")
  if (system(cmd2) != 0) stop("SAIGE step2 failed")

  # 3) Concatenate & filter “Cauchy” rows
  concat_cmd <- sprintf(
    "cd %s && \
     cat step2_out_*_chr{1..22}.txt > all_chr.txt && \
     awk -F '\\t' '$2==\"Cauchy\" {print $1, $4}' all_chr.txt | sort -k2,2g > saige_summary.txt",
    outputDir
  )
  message("→ Concatenating and filtering SAIGE results…")
  if (system(concat_cmd) != 0) stop("Concatenation/filtering failed")

  # 4) Read & return summary
  sumfile <- file.path(outputDir, "saige_summary.txt")
  stopifnot(file.exists(sumfile))
  res <- fread(sumfile, col.names = c("Gene","p.value"))
  return(res)
}
