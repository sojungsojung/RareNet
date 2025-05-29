#' Run SAIGE-GENE+ on a phenotype (step1 → step2 → concatenate)
#'
#' @param phenoFile    Path to phenotype file (tab-delimited)
#' @param genotypeFiles Character vector of genotype file prefixes (PLINK .bed/.bim/.fam or BGEN)
#' @param groupFile    Path to gene-set file (default: shipped geneset_string_v12.txt)
#' @param outputDir    Directory to store all SAIGE outputs
#' @param threads      Number of threads for SAIGE (default: 4)
#' @return data.table with columns Gene, p.value
#' @export
run_saige_gene <- function(phenoFile,
                           genotypeFiles,
                           groupFile    = system.file("data","geneset_string_v12.txt", package="RareNet"),
                           outputDir,
                           threads      = 4) {
  # load required lib for reading results
  library(data.table)

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  # 1) Fit null GLMM (step1)
  step1_prefix <- file.path(outputDir, "step1")
  cmd1 <- sprintf(
    "docker run -v /data:/data wzhou88/saige:1.4.4 step1_fitNULLGLMM.R \
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
    genotypeFiles[1], phenoFile, step1_prefix, threads
  )
  if (system(cmd1) != 0) stop("SAIGE step1 failed")  # :contentReference[oaicite:3]{index=3}

  # 2) Run SPAtests per chromosome (step2)
  cmd2 <- sprintf(
    "docker run -v /data:/data wzhou88/saige:1.4.4 step2_SPAtests.R \
     --bedFile=%s \
     --bimFile=%s \
     --famFile=%s \
     --varianceRatioFile=%s.varianceRatio.txt \
     --GMMATmodelFile=%s.rda \
     --groupFile=%s \
     --outputPrefix=%s/step2_out \
     --nThreads=%d \
     --LOCO=FALSE",
    sub("\\.bed$","", genotypeFiles[1]),                                 # :contentReference[oaicite:4]{index=4}
    paste0(sub("\\.bed$","", genotypeFiles[1]), ".bim"),
    paste0(sub("\\.bed$","", genotypeFiles[1]), ".fam"),
    step1_prefix, step1_prefix,
    groupFile,
    outputDir,
    threads
  )
  if (system(cmd2) != 0) stop("SAIGE step2 failed")

  # 3) Concatenate & filter only “Cauchy” rows
  concat_cmd <- sprintf(
    "cd %s && \
     cat step2_out_*_chr{1..22}.txt > all_chr.txt && \
     awk -F '\\t' '$2==\"Cauchy\" {print $1, $4}' all_chr.txt | sort -k2,2g > saige_summary.txt",
    outputDir
  )
  if (system(concat_cmd) != 0) stop("Concatenation/filtering failed")

  # 4) Read and return
  summary_file <- file.path(outputDir, "saige_summary.txt")
  if (!file.exists(summary_file)) stop("SAIGE summary missing: ", summary_file)
  res <- fread(summary_file, col.names = c("Gene","p.value"))
  return(res)
}
