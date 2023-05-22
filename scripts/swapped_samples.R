library(tidyverse)
raw <- read_tsv(file = "data/Master_counts.txt")
swap <- read_tsv(file =  "data/Master_counts_swap.txt")

indexFD <- match(swap, raw)

swapped_samples <- data.frame(Swap = colnames(swap),
                              Raw = colnames(raw)[indexFD])

write_tsv(swapped_samples, file = "data/swapped_samples.txt")
