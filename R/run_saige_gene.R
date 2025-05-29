#' Run SAIGE-GENE+ on a phenotype (step1 → step2 → concat)
#'
#' @param phenoFile    Path to phenotype file (tab-delimited)
#' @param plinkFile    Prefix for PLINK files (no “.bed/.bim/.fam” extension)
#' @param bedFile      Path to .bed file for step2
#' @param bimFile      Path to .bim file for step2
#' @param famFile      Path to .fam file for step2
#' @param groupFile    Path to gene‐set (group) file
#' @param outputDir    Directory to hold all SAIGE intermediates & summary
#' @param threads      Number of threads for SAIGE (default 4)
#' @return data.table with columns Gene, p.value (SAIGE summary)
#' @export
run_saige_gene <- function(phenoFile,
                           plinkFile,
                           bedFile,
                           bimFile,
                           famFile,
                           groupFile,
                           outputDir,
                           threads = 4) {
  stopifnot(file.exists(phenoFile),
            file.exists(paste0(plinkFile, ".bed")),
            file.exists(bedFile), file.exists(bimFile), file.exists(famFile),
            file.exists(groupFile))
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  # ── Step 1: Fit null GLMM ────────────────────────────────────────────────
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
    plinkFile, phenoFile, step1_prefix, threads
  )
  message("Running SAIGE step1 (null model)...")
  if (system(cmd1) != 0) {
    stop("SAIGE step1 failed")
  }

  # ── Step 2: SPA tests per‐chromosome ──────────────────────────────────────
  message("Running SAIGE step2 (SPAtests)...")
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
    bedFile, bimFile, famFile,
    step1_prefix, step1_prefix,
    groupFile,
    outputDir, threads
  )
  if (system(cmd2) != 0) {
    stop("SAIGE step2 failed")
  }

  # ── Step 3: Concatenate & filter Cauchy results ─────────────────────────
  message("Concatenating per-chromosome SPAtest outputs...")
  concat_sh <- sprintf(
    "cd %s && \\
     cat step2_out_*_chr{1..22}.txt > all_chr.txt && \\
     awk -F '\\t' '$2==\"Cauchy\" {print $1, $4}' all_chr.txt | sort -k2,2g > saige_summary.txt",
    outputDir
  )
  if (system(concat_sh) != 0) {
    stop("Concatenation/filtering failed")
  }

  # ── Step 4: Read & return summary ────────────────────────────────────────
  summary_file <- file.path(outputDir, "saige_summary.txt")
  if (!file.exists(summary_file)) {
    stop("Expected summary not found: ", summary_file)
  }
  res <- data.table::fread(summary_file,
                           col.names = c("Gene", "p.value"))
  return(res)
}
