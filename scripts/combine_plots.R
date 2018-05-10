#!/usr/bin/env Rscript
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
out = args[[1]]
i = 0
df = args[-1] %>% 
     map(read_csv) %>%
     map2(seq_along(.), ~ mutate(.x, i = .y)) %>%
     reduce(rbind)
write_csv(df, str_c(out, ".csv"))
df_truth = df %>% filter(method == "truth")
df_est = df %>% filter(method != "truth")
p = ggplot(df_est, aes(t, Ne, group=interaction(i, method), color=method)) + geom_step(alpha=1/2) + 
      scale_x_log10() + scale_y_log10() + theme_minimal() + geom_step(data=df_truth, color="black", 
                                                                      linetype="dashed")
ggsave(out, p)
