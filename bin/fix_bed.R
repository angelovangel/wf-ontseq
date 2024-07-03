#!/usr/bin/env Rscript

# fix bed file where start > end for plasmids
# arg[1] is bed, arg[2] is len
library(dplyr)

arg <- commandArgs(trailingOnly = T)

df <- read.delim(arg[1], header = F, sep = '\t')
df1 <- df %>%
  filter(V2 < V3)
# these have to be regenerated in 2 rows - 1..end and start..len
df2 <- df %>%
  filter(V2 > V3)
df2a <- df2 %>%
  mutate(V2 = 1)
df2b <- df2 %>%
  mutate(V3 = as.numeric(arg[2]))

finaldf <- bind_rows(df1, df2a, df2b)
write.table(finaldf, file = "", sep = '\t', row.names = F, col.names = F, quote = F)


