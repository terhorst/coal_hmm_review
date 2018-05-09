#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
size_hist<-read.table(args[[1]], header=FALSE)[, 1:2]
names(size_hist) = c('t', 'Ne')
size_hist$method = "psmc"
readr::write_csv(size_hist, args[[2]])
