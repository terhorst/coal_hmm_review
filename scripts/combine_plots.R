#!/usr/bin/env Rscript
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
out = args[[1]]
i = 0
df = map(read_csv, args[-1]) %>% 
     map2(seq_along(args[1]), ~ mutate(.x, i=.y)) %>%
     reduce(rbind)
ggplot(df, aes(t, Ne, group=i, color=method)) + geom_step()
ggsave(out)
