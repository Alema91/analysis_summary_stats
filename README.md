# Creación del summary_stats de viralrecon sin multiqc

Un análisis de viralrecon con más de 100 muestras podría provocar problemas para generar el MultiQC y por consiguiente el archivo summary_mqc. Este archivo incluye varios datos necesarios para la generación de la summary_stats. Por lo tanto, el objetivo de este repositorio es documentar y programar un método para la obtención de la summary_stats sin necesidad de utilizar el MultiQC.

## Columnas que se sustituyen:

- Coverage median
- % Coverage > 10x
- Read virus o Mapped reads

# lablog

"""
# copiamos los scripts para calcular la coverage con mosdepth
rsync -rlv /data/bi/pipelines/script_summary_stats_sinmultiqc/create_df_coverage.R .
rsync -rlv /data/bi/pipelines/script_summary_stats_sinmultiqc/create_summary_mod_final.sh create_summary_nomultiqc.sh

# generamos el quast con los consensos
#module load QUAST/5.0.2-foss-2020a-Python-3.8.2

scratch_dir=$(echo $PWD | sed "s/\/data\/bi\/scratch_tmp/\/scratch/g")
CONSENSUS="*/variants/ivar/consensus/bcftools/"

echo "srun --partition short_idx --output QUAST.%j.log --chdir ${scratch_dir} quast $CONSENSUS/*.consensus.fa -o quast_results" > _05_quast.sh

# Lanzamos el script para el calculo de la coverage

#module load R/4.1.0-foss-2021a
echo "srun --partition short_idx --output COVERAGE.%j.log --chdir ${scratch_dir} Rscript create_df_coverage.R -i */variants/bowtie2/mosdepth/genome/all_samples.mosdepth.coverage.tsv -o df_coverage.csv" > _06_coverage.sh

# Summary stats

echo "bash create_summary_nomultiqc.sh" > _07_create_summary_report_nomultiqc.sh
"""



