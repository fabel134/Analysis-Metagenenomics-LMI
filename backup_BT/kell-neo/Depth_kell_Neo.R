library(ggplot2)
library(dplyr)


library(readr)
LP_all_sort <- read_delim("LP_all_sort.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE , show_col_types = FALSE)
ggplot(LP_all_sort, aes(x = Position , y = Depth)) +
  geom_bar(stat = "identity", fill = "green") +
  labs(x = "Position", y = "Depth", title = "Kellermania") +
  theme_minimal() +
  facet_wrap(~Reference)


LP_all_neo_sort <- read_delim("LP_all_neo_sort.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE) 

ggplot(LP_all_neo_sort, aes(x = Position , y = Depth)) +
  geom_bar(stat = "identity", fill = "green") +
  labs(x = "Position", y = "Depth", title = "Neophaeosphaeria") +
  theme_minimal() +
  facet_wrap(~Reference)

