#!/bin/bash

echo -e "Sample\t# Input reads\t# Trimmed reads (fastp)\t% Non-host reads (Kraken 2)\t% Mapped reads\t# Mapped reads\t# Trimmed reads (iVar)\tCoverage median\t% Coverage > 1x\t% Coverage > 10x\t# SNPs\t# INDELs\t# Missense variants\t# Ns per 100kb consensus\tPangolin lineage\tNextclade clade" > header_mqc.csv
tail -n +2 summary_mqc_filter.csv > tmp_mqc.csv
cat header_mqc.csv tmp_mqc.csv > summary_mqc_filter.csv
rm tmp_mqc.csv header_mqc.csv

