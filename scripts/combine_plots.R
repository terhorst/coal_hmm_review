#!/usr/bin/env Rscript
library(scales)
library(tidyverse)
library(ggthemes)
library(lettercase)
args = commandArgs(trailingOnly=TRUE)
out = args[[1]]
i = 0
df = args[-1] %>% 
     map(read_csv) %>%
     map2(seq_along(.), ~ mutate(.x, i = .y)) %>%
     reduce(rbind) %>%
     mutate(Method = ifelse(stringr::str_to_lower(method) == "dical", "diCal", method),
            Demography = str_title_case(demo))
write_csv(df, stringr::str_c(out, ".csv"))
df_truth = df %>% filter(stringr::str_sub(method, 0, 5) == "truth")
df_est = df %>% filter(stringr::str_sub(method, 0, 5) != "truth")
extended = df_est %>% group_by(method) %>% top_n(1, t) %>% mutate(t=max(t, 1e7)) %>% bind_rows(df_est)
# extended = mutate(extended, t = pmax(t, 10)) %>% print(n=Inf)
df_truth = mutate(df_truth, t = pmax(t, 10))
p = ggplot(extended, aes(t, Ne, group=interaction(i, Method), color=Method)) + geom_step(alpha=1/3) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
              labels = trans_format("log10", math_format(10^.x))) +
      theme_minimal() + geom_step(data=df_truth, color="black", linetype="dashed") +
      scale_color_brewer(palette="Set1") + ylab(expression(N[e])) +
      xlab("Time (years; g=29)") + facet_wrap(~ Demography, nrow = 2) +
      guides(colour = guide_legend(override.aes = list(alpha = 1)))  # do not make lines illegible in legend
ggsave(out, p, width=4, height=3, units="in")
