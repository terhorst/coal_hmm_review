#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
mu <- args[[1]]
gen <- args[[2]],
size_hist<-read.table(args[[3]], header=TRUE)
size_hist$method = "msmc"
readr::write_csv(size_hist, args[[4]])
