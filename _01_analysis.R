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
