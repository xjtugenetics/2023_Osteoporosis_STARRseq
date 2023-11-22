ls *fastq|sed 's/^.*\///'|while read m;do bowtie2 -q $m -x $pwd/bowtie2.index/6197.SNPs.out -p 8 --local --xeq > $m\.out;done
