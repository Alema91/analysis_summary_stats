# library -------

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(readxl, quietly = TRUE, warn.conflicts = FALSE)
library(forcats, quietly = TRUE, warn.conflicts = FALSE)
library(reshape2, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)
library(rmarkdown, quietly = TRUE, warn.conflicts = FALSE)
library(knitr, quietly = TRUE, warn.conflicts = FALSE)

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

# Ns per 100 kb
tabla_ns <- read.csv2(opt$file, header = T, sep = "\t")
mod_tabla_ns <- data.frame(samples = tabla_ns$Assembly, Ns100kb = tabla_ns$X..N.s.per.100.kbp)
mod_tabla_ns$samples_mod <- str_split(mod_tabla_ns$samples, ".consensus", n = Inf, simplify = T)[, 1]
# creamos la tabla
df_ns <- data.frame(samples = mod_tabla_ns$samples_mod, Nsper100kb = mod_tabla_ns$Ns100kb)
write.table(x = df_ns, file = opt$out, quote = F, row.names = F, col.names = T, sep = "\t")
