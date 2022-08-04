echo -e "run\tuser\thost\tVirussequence\tsample\ttotalreads\treadshostR1\treadshost\t%readshost\treadsvirus\t%readsvirus\tunmappedreads\t%unmapedreads\tmedianDPcoveragevirus\tCoverage>10x(%)\tInputreads\tTrimmedreads\tNsper100kb\tVariantsinconsensusx10\tMissenseVariants\t%Ns10x\tLineage\tPangolinversion" > mapping_illumina_mod.tab

USER=$(pwd | cut -d '/' -f6 | cut -d '_' -f4)
HOST=$(pwd | cut -d '/' -f8 | cut -d '_' -f4 | tr '[:upper:]' '[:lower:]' | sed 's/.*/\u&/')

cat samples_id_v2.txt | while read in line; do run="ls -lt ../../RAW/${line}*_R1*_*.fastq.gz | cut -d'/' -f7 | sort -u | head -n1" \
sample="${line}" \
echo $sample\t$run>> run_id.txt; done

cat samples_id_mod.txt | while read in
do
	echo -e "${RUN}\t${USER}\t${HOST}\tNC_045512.2\t${in}\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3)\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print ($1*2)}')\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f1 | sed 's/ //g')\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)\t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print ($2*100)/$1'})")\t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}')\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print $1-($2+$3)}')")\t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}')\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print (($1-($2+$3))/$1)*100}')")\t$(cat df_coverage.txt | grep ${in} | cut -f2)\t$(cat df_coverage.txt | grep ${in} | cut -f3)\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n1 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')\t$(cat df_nsper100.csv | grep ${in} | cut -f2)\t$(zcat */variants/ivar/consensus/bcftools/${in}.filtered.vcf.gz | grep -v '^#' | wc -l)\t$(LC_ALL=C awk -F, '{if($10 >= 0.75)print $0}' */variants/ivar/variants_long_table.csv | grep -w ${in} | grep 'missense' | wc -l)\t$(cat %Ns.tab | grep ${in} | cut -f2)\t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f2)\t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f10)"  >> mapping_illumina_mod.tab
done


