#!/usr/bin/env Rscript
library(tidyverse, quietly = T, warn.conflicts = F)

txt = paste(readLines('stdin'), collapse="\n")
df = read_delim(txt, " ", col_names = F, n_max = 2)
theta = df[[1,2]]
N0 = theta / 2
df = read_delim(txt, "\t", col_names = F, skip = 3) %>% 
    tidyr::separate(X1, into=c("a", "t"), sep=" ") %>% 
    dplyr::mutate(t = as.numeric(t)) %>%
    dplyr::mutate(left_time_boundary = t * 2 * N0, lambda = (1 / X2 / N0))
df %>% select(left_time_boundary, lambda) %>% 
    format_delim(delim = " ", col_names = T) %>% cat
