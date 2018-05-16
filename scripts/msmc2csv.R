#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
mu <- as.numeric(args[[1]])
gen <- as.numeric(args[[2]])
demo <- args[[3]]
size_hist<-read.table(args[[4]], header=TRUE)
out = dplyr::transmute(size_hist, 
                    t = left_time_boundary / mu * gen,
                    Ne = (1 / lambda_00) / mu,
                    method = "MSMC",
                    demo = demo)
readr::write_csv(out, args[[5]])
