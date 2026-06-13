#!/usr/bin/env Rscript
# SuSiE-RSS fine-mapping + coloc.susie for ONE GWAS x eQTL locus.
# Shipped as a resource and invoked by finemap.py via subprocess.
#
# Usage: Rscript susie_coloc.R <work_dir> [r_libs_path]
#   <work_dir>/merged.tsv : snp,pos,beta_gwas,se_gwas,beta_eqtl,se_eqtl
#   <work_dir>/ld.tsv     : signed LD correlation matrix (ALT-dosage; same order)
#   <work_dir>/meta.json  : {"gwas_N":..., "eqtl_N":...}
# Emits machine-parseable "KEY value" lines (CS_GWAS, CS_EQTL,
# COLOC_SUSIE_PP4, LD_S_GWAS, LD_S_EQTL, DONE) on stdout.

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
ul <- if (length(args) >= 2 && nzchar(args[2])) args[2] else Sys.getenv("R_LIBS_USER")
if (nzchar(ul)) .libPaths(c(ul, .libPaths()))

ok <- requireNamespace("susieR", quietly = TRUE) &&
      requireNamespace("coloc", quietly = TRUE) &&
      requireNamespace("jsonlite", quietly = TRUE)
if (!ok) { cat("ERROR missing_r_packages\n"); quit(status = 2) }
suppressMessages({ library(susieR); library(coloc); library(jsonlite) })

m  <- read.delim(file.path(wd, "merged.tsv"))
R  <- as.matrix(read.table(file.path(wd, "ld.tsv")))
mt <- fromJSON(file.path(wd, "meta.json"))
rownames(R) <- colnames(R) <- m$snp
zg <- m$beta_gwas / m$se_gwas
ze <- m$beta_eqtl / m$se_eqtl

sg <- tryCatch(suppressWarnings(susie_rss(z = zg, R = R, n = mt$gwas_N, L = 10)),
               error = function(e) NULL)
se <- tryCatch(suppressWarnings(susie_rss(z = ze, R = R, n = mt$eqtl_N, L = 10)),
               error = function(e) NULL)
csg <- if (is.null(sg) || is.null(sg$sets$cs)) 0L else length(sg$sets$cs)
cse <- if (is.null(se) || is.null(se$sets$cs)) 0L else length(se$sets$cs)
cat(sprintf("CS_GWAS %d\n", csg))
cat(sprintf("CS_EQTL %d\n", cse))

pp4 <- 0.0
if (csg > 0 && cse > 0) {
  cr <- tryCatch(coloc.susie(sg, se), error = function(e) NULL)
  if (!is.null(cr) && !is.null(cr$summary) && nrow(cr$summary) > 0)
    pp4 <- max(cr$summary$PP.H4.abf)
}
cat(sprintf("COLOC_SUSIE_PP4 %.4f\n", pp4))

sgd <- tryCatch(estimate_s_rss(zg, R, n = mt$gwas_N), error = function(e) NA)
sed <- tryCatch(estimate_s_rss(ze, R, n = mt$eqtl_N), error = function(e) NA)
cat(sprintf("LD_S_GWAS %s\n", ifelse(is.na(sgd), "NA", sprintf("%.3f", sgd))))
cat(sprintf("LD_S_EQTL %s\n", ifelse(is.na(sed), "NA", sprintf("%.3f", sed))))
cat("DONE\n")
