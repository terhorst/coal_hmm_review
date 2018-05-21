#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
size_hist<-read.table(args[[1]], header=FALSE)[, 1:2]
names(size_hist) = c('t', 'Ne')
size_hist$method = "PSMC"
size_hist$Ne = size_hist$Ne * 1e4
size_hist$demo = args[[2]]
readr::write_csv(size_hist, args[[3]])
