echo -e "run\tuser\thost\tVirussequence\tsample\ttotalreads\treadshostR1\treadshost\t%readshost\treadsvirus\t%readsvirus\tunmappedreads\t%unmapedreads\tmedianDPcoveragevirus\tCoverage>10x(%)\tInputreads\tTrimmedreads\tNsper100kb\tVariantsinconsensusx10\tMissenseVariants\t%Ns10x\tLineage" > mapping_illumina_mod.tab
echo -e "Sample\t# Input reads\t# Trimmed reads (fastp)\t% Non-host reads (Kraken 2)\t% Mapped reads\t# Mapped reads\t# Trimmed reads (iVar)\tCoverage median\t% Coverage > 1x\t% Coverage > 10x\t# SNPs\t# INDELs\t# Missense variants\t# Ns per 100kb consensus\tPangolin lineage\tNextclade clade" > summary_mqc_mod.csv

USER=$(pwd | cut -d '/' -f6 | cut -d '_' -f4)
HOST=$(pwd | cut -d '/' -f8 | cut -d '_' -f4 | tr '[:upper:]' '[:lower:]' | sed 's/.*/\u&/')

# he creado un bucle temporal para generar los nombres de las carreras

for line in `cat samples_id.txt`;do
	run=`ls -lt ../../RAW/${line}_*R1*.fastq.gz | cut -d "/" -f7`
	sample=${line}
	echo -e "$sample\t$run" >> run_tmp.txt
done

# generamos la mapping_illumina & summary_mqc

cat samples_id_mod.txt | while read in
do
	echo -e "$(cat run_tmp.txt | grep -w ${in} | cut -f2)\t${USER}\t${HOST}\tNC_045512.2\t$(echo ${in} | sed 's/_/-/g')\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3)\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print ($1*2)}')\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f1 | sed 's/ //g')\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)\t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print ($2*100)/$1'})")\t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}')\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print $1-($2+$3)}')")\t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}')\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print (($1-($2+$3))/$1)*100}')")\t$(cat df_coverage.txt | grep ${in} | cut -f2)\t$(cat df_coverage.txt | grep ${in} | cut -f3)\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n1 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat df_nsper100.csv | grep ${in} | cut -f2)\t$(zcat */variants/ivar/consensus/bcftools/${in}.filtered.vcf.gz | grep -v '^#' | wc -l)\t$(LC_ALL=C awk -F, '{if($10 >= 0.75)print $0}' */variants/ivar/variants_long_table.csv | grep -w ${in} | grep 'missense' | wc -l)\t$(cat %Ns.tab | grep ${in} | cut -f2)\t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f2)" >> mapping_illumina_mod.tab
	echo -e "$(echo ${in} | sed 's/_/-/g')\t$(echo 'NA')\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n1 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(cat quast_results/quast.log | grep "${in}.consensus" | cut -d "#" -f2 | grep "N's per 100 kbp" | cut -d "=" -f2)\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')" >> summary_mqc_mod.csv
done

rm -rf run_tmp.txt quast_results

################################################################################################################################################################################
# Prueba mqc####################################################################################################################################################################
################################################################################################################################################################################

echo -e "Sample\t# Input reads\t# Trimmed reads (fastp)\t% Non-host reads (Kraken 2)\t% Mapped reads\t# Mapped reads\t# Trimmed reads (iVar)\tCoverage median\t% Coverage > 1x\t% Coverage > 10x\t# SNPs\t# INDELs\t# Missense variants\t# Ns per 100kb consensus\tPangolin lineage\tNextclade clade" > prueba_mqc.csv
cat samples_id_mod.txt | while read in
do
	echo -e "$(echo ${in} | sed 's/_/-/g')\t$(echo 'NA')\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n1 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(cat quast_results/quast.log | grep "${in}.consensus" | cut -d "#" -f2 | grep "N's per 100 kbp" | cut -d "=" -f2)\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')\t$(echo 'NA')" >> prueba_mqc.csv
done