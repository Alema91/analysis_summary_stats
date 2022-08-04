##############################################################################################################################
# modificado
##############################################################################################################################

echo -e "run\tuser\thost\tVirussequence\tsample\ttotalreads\treadshostR1\treadshost\t%readshost\treadsvirus\t%readsvirus\tunmappedreads\t%unmapedreads\tmedianDPcoveragevirus\tCoverage>10x(%)\tVariantsinconsensusx10\tMissenseVariants\t%Ns10x\tLineage\tPangolinversion" > mapping_illumina_v2.tab
USER=$(pwd | cut -d '/' -f6 | cut -d '_' -f4)
HOST=$(pwd | cut -d '/' -f8 | cut -d '_' -f4 | tr '[:upper:]' '[:lower:]' | sed 's/.*/\u&/')

cat samples_id_v2_mod.txt | while read in
do
	# run
    echo -e "$(cat samples_id_v2.txt | xargs -I % ls -lt ../../RAW/%* | bash | cut -d'/' -f7 | sort -u | head -n1) \
    
    # usuario
    \t${USER}
    
    # host
    \t${HOST}
    
    # Virussequence
    \tNC_045512.2
    
    # samples
    \t${in}
    
    # totalreads
    \t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \

    # readshostR1
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3) \

    # readshost
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print ($1*2)}') \

    # %readshost
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f1 | sed 's/ //g') \

    # readsvirus
    \t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1) \

    # %readsvirus
    \t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \
    \t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print ($2*100)/$1'})") \

    # unmappedreads
    \t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}') \
    \t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print $1-($2+$3)}')") \

    # %unmapped
    \t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}') \
    \t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print (($1-($2+$3))/$1)*100}')") \

    # mediancoverage
    \t$(cat df_coverage.txt | grep ${in} | cut -f2)) \

    # %coverage10x
    \t$(cat df_coverage.txt | grep ${in} | cut -f3)) \

    # inputreads
    \t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n1 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \

    # trimmedreads
    \t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \

    # Nsper100kb
    \t$(cat df_nsper100.csv | grep ${in} | cut -f2)) \

    # variants
    \t$(zcat */variants/ivar/consensus/bcftools/${in}.filtered.vcf.gz | grep -v '^#' | wc -l) \

    # missense
    \t$(LC_ALL=C awk -F, '{if($10 >= 0.75)print $0}' */variants/ivar/variants_long_table.csv | grep -w ${in} | grep 'missense' | wc -l) \

   # Ns
    \t$(cat %Ns.tab | grep ${in} | cut -f2)\t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f2)

    # lineage
    \t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f2)  
    
    # version pangolin

    \t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f10)
    
    >> mapping_illumina_v2.tab

done

##############################################################################################################################
# original
##############################################################################################################################

# user
echo -e "\t${USER}

# host
\t${HOST}

# virus sequence
\tNC_045512.2

# samples
\t${in}

# totalreads
\t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')

#readhostR1
\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3)

# readhost
\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print ($1*2)}')

# %readhost
\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f1 | sed 's/ //g')

# readvirus
\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)

# %readsvirus
\t$(cat */multiqc/summary_variants_metrics_mqc.csv | grep ^${in}, | cut -d ',' -f5)

# unmapped
\t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')
\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}')
\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print $1-($2+$3)}')")

# %unmapped
\t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g')
\t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}')
\t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print (($1-($2+$3))/$1)*100}')")

# medianDP
\t$(cat */multiqc/summary_variants_metrics_mqc.csv | grep ^${in}, | cut -d ',' -f8)

# %coverage10x
\t$(cat */multiqc/summary_variants_metrics_mqc.csv | grep ^${in}, | cut -d ',' -f10)

# variants
\t$(zcat */variants/ivar/consensus/bcftools/${in}.filtered.vcf.gz | grep -v '^#' | wc -l)

# missense
\t$(LC_ALL=C awk -F, '{if($10 >= 0.75)print $0}' */variants/ivar/variants_long_table.csv | grep -w ${in} | grep 'missense' | wc -l)

# Ns
\t$(cat %Ns.tab | grep ${in} | cut -f2)\t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f2)

# lineage
\t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f10)"  >> mapping_illumina.tab


##############################################################################################################################
# final
##############################################################################################################################


echo -e "run\tuser\thost\tVirussequence\tsample\ttotalreads\treadshostR1\treadshost\t%readshost\treadsvirus\t%readsvirus\tunmappedreads\t%unmapedreads\tmedianDPcoveragevirus\tCoverage>10x(%)\tInputreads\tTrimmedreads\tNsper100kb\tVariantsinconsensusx10\tMissenseVariants\t%Ns10x\tLineage\tPangolinversion" > mapping_illumina_v2.tab
USER=$(pwd | cut -d '/' -f6 | cut -d '_' -f4)
HOST=$(pwd | cut -d '/' -f8 | cut -d '_' -f4 | tr '[:upper:]' '[:lower:]' | sed 's/.*/\u&/')

cat samples_id_v2_mod.txt | while read in
do
    echo -e "$(cat samples_id_v2.txt | xargs -I % ls -lt ../../RAW/%* | cut -d'/' -f7 | sort -u | head -n1) \
    
    # usuario
    \t${USER}
    
    # host
    \t${HOST}
    
    # Virussequence
    \tNC_045512.2
    
    # samples
    \t${in}
    
    # totalreads
    \t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \

    # readshostR1
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3) \

    # readshost
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print ($1*2)}') \

    # %readshost
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f1 | sed 's/ //g') \

    # readsvirus
    \t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1) \

    # %readsvirus
    \t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \
    \t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print ($2*100)/$1'})") \

    # unmappedreads
    \t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}') \
    \t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print $1-($2+$3)}')") \

    # %unmapped
    \t$(echo "$(echo -e "$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \
    \t$(cat */kraken2/${in}.kraken2.report.txt | tail -n1 | cut -f3| awk '{print $1*2}') \
    \t$(cat */variants/bowtie2/samtools_stats/${in}.sorted.bam.flagstat | grep '+ 0 mapped' | cut -d ' ' -f1)" | awk '{print (($1-($2+$3))/$1)*100}')") \

    # mediancoverage
    \t$(cat df_coverage.txt | grep ${in} | cut -f2)) \

    # %coverage10x
    \t$(cat df_coverage.txt | grep ${in} | cut -f3)) \

    # Nsper100kb
    \t$(cat df_nsper100.csv | grep ${in} | cut -f2)) \

    # inputreads
    \t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n1 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \

    # totalreads
    \t$(grep 'total_reads' */fastp/${in}.fastp.json | head -n2 | tail -n1 | cut -d ':' -f2 | sed 's/,//g') \

    # variants
    \t$(zcat */variants/ivar/consensus/bcftools/${in}.filtered.vcf.gz | grep -v '^#' | wc -l) \

    # missense
    \t$(LC_ALL=C awk -F, '{if($10 >= 0.75)print $0}' */variants/ivar/variants_long_table.csv | grep -w ${in} | grep 'missense' | wc -l) \

   # Ns
    \t$(cat %Ns.tab | grep ${in} | cut -f2)\t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f2)

    # lineage
    \t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f2)  
    
    # version pangolin

    \t$(cat */variants/ivar/consensus/bcftools/pangolin/${in}.pangolin.csv | tail -n1 | cut -d ',' -f10)
    
    >> mapping_illumina_v2.tab

done
