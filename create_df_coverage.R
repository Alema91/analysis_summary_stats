# library -------

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(readxl, quietly = TRUE, warn.conflicts = FALSE)
library(forcats, quietly = TRUE, warn.conflicts = FALSE)
library(reshape2, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)

# opt

option_list <- list(
    make_option(c("-f", "--file"),
        type = "character", default = NULL,
        help = "dataset file name", metavar = "character"
    ),
    make_option(c("-o", "--out"),
        type = "character", default = "out.txt",
        help = "output file name [default= %default]", metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#### Coverage -----

# list_coverage <- list.files("./", pattern = "all_samples.mosdepth.coverage.tsv", full.names = T)
tabla_coverage <- read.table(opt$file, header = T)

df_coverage <- matrix(0, nrow = length(unique(tabla_coverage$sample)), ncol = 3)
samples <- unique(tabla_coverage$sample)
for (i in 1:length(samples)) {
    df_tmp <- tabla_coverage[tabla_coverage$sample == samples[i], ]
    df_tmp$coord <- abs(df_tmp[, 2] - df_tmp[, 3])
    # mayor10 <- round((length(df_tmp[df_tmp$coverage > 10, 2]) / length(df_tmp$chrom) * 100), 2)
    mayor10 <- round(sum(df_tmp[df_tmp$coverage > 10, 6]) / sum(df_tmp$coord) * 100, 2)
    mediancoverage <- round(median(df_tmp$coverage), 2)
    # rellenamos la tabla
    df_coverage[i, 1] <- unique(df_tmp$sample)
    df_coverage[i, 2] <- mediancoverage
    df_coverage[i, 3] <- mayor10
}

# creamos la tabla
df_coverage <- data.frame(df_coverage)
colnames(df_coverage) <- c("sample", "medianCoverage", "%Coverage>10x")
write.table(x = df_coverage, file = opt$out, quote = F, row.names = F, col.names = T, sep = "\t")
