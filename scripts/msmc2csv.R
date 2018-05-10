#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
mu <- as.numeric(args[[1]])
gen <- as.numeric(args[[2]])
size_hist<-read.table(args[[3]], header=TRUE)
out = dplyr::transmute(size_hist, 
                    t = left_time_boundary / mu * gen,
                    Ne = (1 / lambda_00) / mu,
                    method = "MSMC")
readr::write_csv(out, args[[4]])
