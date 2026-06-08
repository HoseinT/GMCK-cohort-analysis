#!/usr/bin/env Rscript
# Run FACETS v0.6.2 CNV analysis on snp-pileup output.
# Outputs:
#   <sample>_cnv.tsv   — segment table with CNV class labels
#   <sample>_purity.txt — purity / ploidy on two lines
#   <sample>_facets_qc.txt — QC metrics

suppressPackageStartupMessages({
  library(facets)
  library(optparse)
})

option_list <- list(
  make_option("--pileup", type = "character", help = "snp-pileup CSV.gz file"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--sample", type = "character", help = "Sample name prefix")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Read and pre-process pileup ----------------------------
rcmat <- readSnpMatrix(opt$pileup)
xx    <- preProcSample(rcmat, gbuild = "hg19")

# ---- Run FACETS ---------------------------------------------
# Use increased max_depth (handled upstream by snp-pileup --max-depth).
# cval = 150 (stringent) for exome/panel data; tune if needed.
oo <- procSample(xx, cval = 150)

# ---- Fit and extract purity/ploidy -------------------------
fit <- emcncf(oo)

purity <- fit$purity
ploidy <- fit$ploidy

# ---- Label CNV classes -------------------------------------
seg <- fit$cncf

label_cnv <- function(tcn, lcn) {
  dplyr::case_when(
    tcn == 0           ~ "deep_loss",
    tcn < 2            ~ "loss",
    tcn == 2 & lcn == 0 ~ "CNLOH",
    tcn == 2 & lcn == 1 ~ "neutral",
    tcn >= 3 & tcn <= 5 ~ "gain",
    tcn >= 6            ~ "amplification",
    TRUE                ~ "unknown"
  )
}

seg$cnv_class <- label_cnv(seg$tcn.em, seg$lcn.em)

# ---- Write outputs -----------------------------------------
out_seg    <- file.path(opt$outdir, paste0(opt$sample, "_cnv.tsv"))
out_purity <- file.path(opt$outdir, paste0(opt$sample, "_purity.txt"))
out_qc     <- file.path(opt$outdir, paste0(opt$sample, "_facets_qc.txt"))

write.table(seg, out_seg, sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c("purity\tploidy", sprintf("%.4f\t%.4f", purity, ploidy)), out_purity)

qc_lines <- c(
  paste("dipLogR",    fit$dipLogR,   sep = "\t"),
  paste("loglik",     fit$loglik,    sep = "\t"),
  paste("flags",      paste(oo$flags, collapse = ";"), sep = "\t")
)
writeLines(qc_lines, out_qc)

message("FACETS done: purity=", round(purity, 3), " ploidy=", round(ploidy, 3))
