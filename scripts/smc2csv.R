#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
size_hist<-readr::read_csv(args[[1]])[, 2:3]
names(size_hist) = c('t', 'Ne')
size_hist$method = "SMC++"
size_hist$demo = args[[2]]
print(size_hist)
readr::write_csv(size_hist, args[[3]])
