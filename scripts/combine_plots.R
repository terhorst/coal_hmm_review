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
df_truth = df %>% filter(str_sub(method, 0, 5) == "truth")
df_est = df %>% filter(str_sub(method, 0, 5) != "truth")
extended = df_est %>% group_by(method) %>% top_n(1, t) %>% mutate(t=max(t, 1e7)) %>% bind_rows(df_est)
# extended = mutate(extended, t = pmax(t, 10)) %>% print(n=Inf)
df_truth = mutate(df_truth, t = pmax(t, 10)) %>% print
p = ggplot(extended, aes(t, Ne, group=interaction(i, method), color=method)) + geom_step(alpha=1/2) + 
      scale_x_log10() + scale_y_log10() + theme_minimal() + geom_step(data=df_truth, color="black", 
                                                                      linetype="dashed")
ggsave(out, p)
