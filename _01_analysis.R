# library -------

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(readxl, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(forcats, quietly = TRUE, warn.conflicts = FALSE)
library(reshape2, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(Metrics, quietly = TRUE, warn.conflicts = FALSE)

### Not in o in ----

`%notin%` <- Negate(`%in%`)

# creamos un indice para la mapping_illumina_mod.xlsx
data_mapping_illumina <- read_xlsx("output/mapping_illumina_mod.xlsx")
data_summary_mqc <- read.csv2("output/summary_mqc_mod.csv", sep = ",")

data_luis <- read_xlsx("data/METADATA_LAB_2022_07_22_Luis.xlsx")

#### cambiamos el _ de las muestras repetidas por guion

mod_samples_id_illu <- gsub("_", "-", data_mapping_illumina$sample)
data_mapping_illumina$sample <- mod_samples_id_illu

mod_samples_id_mqc <- gsub("_", "-", data_summary_mqc$Sample)
data_summary_mqc$Sample <- mod_samples_id_mqc

#### filtramos las muestras de luis del mapping_illumina

id_luis <- as.character(data_luis$"Sample ID given for sequencing")

mapping_filter <- data_mapping_illumina[data_mapping_illumina$sample %in% id_luis, ]
mqc_filter <- data_summary_mqc[data_summary_mqc$Sample %in% id_luis, ]

#### Write data_filter

write.table(mapping_filter, "output/filter/mapping_illumina_filter.xlsx", sep = "\t", row.names = F, quote = F)
write.table(mqc_filter, "output/filter/summary_mqc_filter.csv", sep = ",", row.names = F, quote = F)


#### Analisis de la variacion en el coverage entre el original del multiqc y mi calculo de la coverage

# Data multiqc
data_mqc_multiqc_v2 <- read.csv2("data/mqc_multiqc/summary_variants_metrics_mqc_v2.csv", sep = ",")
data_mqc_multiqc_v4 <- read.csv2("data/mqc_multiqc/summary_variants_metrics_mqc_v4.csv", sep = ",")

# modificar tablas
data_mqc_multiqc <- rbind(data_mqc_multiqc_v2, data_mqc_multiqc_v4)
data_mqc_multiqc$Sample <- gsub("_", "-", data_mqc_multiqc$Sample)
filter_multiqc_coverage <- data_mqc_multiqc[, c(1, 8, 10)]
colnames(filter_multiqc_coverage) <- c("sample", "median_multiqc", "%10_multiqc")

data_coverage <- read.csv2("data/df_coverage.txt", sep = "\t")
filter_data_coverage <- data_coverage[data_coverage$sample %in% data_mqc_multiqc$Sample, ]
colnames(filter_data_coverage) <- c("sample", "median_bu", "%10_bu")

# diferencias
samples_dif <- setdiff(filter_multiqc_coverage$sample, filter_data_coverage$sample)
filter_2_multiqc_coverage <- filter_multiqc_coverage[filter_multiqc_coverage$sample %notin% samples_dif, ]

# Data final
df_coverage <- merge(filter_2_data_coverage, filter_2_multiqc_coverage)

##### Analisis ----

# df long
df_coverage_long <- melt(df_coverage, id.var = c("sample"))
df_coverage_long$value <- as.numeric(df_coverage_long$value)
df_coverage_long$value[is.na(df_coverage_long$value)] <- 0

# mediana
df_median <- df_coverage_long[df_coverage_long$variable == c("median_bu", "median_multiqc"), ]
df_median$variable <- factor(df_median$variable)

# test wilconxon (Mann-Whitney test)
wilcox.test(value ~ variable, data = df_median)

# RSME
rmse(as.numeric(df_coverage$"%10_bu"), as.numeric(df_coverage$"%10_multiqc"))

# Plots
ggplot(df_median, aes(value, color = variable)) +
    geom_histogram(aes(y = ..density..), color = "black", fill = "white") +
    geom_density(alpha = .2, fill = "#FF6666")

# %10X
df_10x <- df_coverage_long[df_coverage_long$variable == c("%10_bu", "%10_multiqc"), ]
df_10x$variable <- factor(df_10x$variable)

# test wilconxon (Mann-Whitney test)
wilcox.test(value ~ variable, data = df_10x)

# Plots
ggplot(df_10x, aes(value, color = variable)) +
    geom_histogram(aes(y = ..density..), color = "black", fill = "white") +
    geom_density(alpha = .2, fill = "#FF6666")
