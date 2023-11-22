#osteoporosis-relevant functional analysis for input SNPs and (if control SNPs provided) epigenetic and allelic predicted motif enrichment analysis
#/bin/bash
 
usage(){
  echo "
Usage:
  -i, --input   input SNPs <vcf format>
  -c, --control control/backgound SNPs for enrichment analysis <vcf format>
  -r, --region expanded region[bp] for SNP epigenetic annotation <numeric eg:60>
  -h, --help    display this help and exit
 
  example1: sh SNP.osteoporosis.functional.analysis.sh -i baaSNPs.vcf -c inactive.SNPs.vcf -r 60
  example2: sh SNP.osteoporosis.functional.analysis.sh --input baaSNPs.vcf --control inactive.SNPs.vcf --region 60
  example3<no control SNPs>: sh SNP.osteoporosis.functional.analysis.sh -i input.SNPs <vcf format>
"
}

main(){
while true
do
  case "$1" in
  -i|--input)
      input="$2"
      shift
      ;;
  -c|--control)
      control="$2"
      shift
      ;;
  -r|--region)
      region="$2"
      shift
      ;;
  -h|--help)
      usage
      exit
      ;;
  --)
    shift
    break
    ;;
  *)
    echo "$1 is not option"
    ;;
  esac
  shift
done
}

set -- $(getopt -o i:c:r:h --long input:,control:,region:,help -- "$@") 
main $@

path=$(dirname "$PWD")
output=`echo $input|sed 's/.vcf//'|sed "s/'//g"`
snp1=`echo $input|sed "s/'//g"`
if [ ${#} == "7" ]; then
	snp2=`echo $control|sed "s/'//g"`
fi
region=`echo $region|sed "s/'//g"`

#epigenetic enrichment analysis
if [ ${#} == "7" ]; then
	python $path/HiC.genes.functional.analysis/SNP.osteoblast.epigenetic.anno.enrichment.py --SNP1 $snp1 --SNP2 $snp2 --region $region --out $output
	mkdir epigenetic.enrichment
	mv ENCODE.osteoblast.plot.pdf ./epigenetic.enrichment/
	mv hMSC.induced.osteoblast.plot.pdf ./epigenetic.enrichment/
	mv osteoblast.HMM15.plot.pdf ./epigenetic.enrichment/
	mv *epigentic.enrichment.Fisher.txt ./epigenetic.enrichment/
	echo "epigenetic enrichment analysis finished!"
else
	echo "no eogentic enrichment analysis becuase control SNPs not provided!"
fi

#extract Hi-C interacted genes & functional annatation
python $path/HiC.genes.functional.analysis/SNP.HiC.interacted.genes.functional.anno.py --SNP $snp1
mkdir HiC.genes.annotation
mv *HiC.genes.anno.txt ./HiC.genes.annotation/
mv top10.GO.enriched.pathways.pdf ./HiC.genes.annotation/
mv HiC.interacted.genes.GO.significant.enrichment.txt ./HiC.genes.annotation/
mv SNPs.osteoblast.HiC.interacted.genes ./HiC.genes.annotation/
mv baaSNPs.HiC.genes.function.summary ./HiC.genes.annotation/

#alleleic motif prediction
if [ ${#} == "7" ]; then
	cat $snp1 $snp2|sort|uniq > meme.suit.SNPs.vcf
	awk '{print $3"\t""'$snp1'"}' $snp1|sed 's/.vcf//' > SNPs.types.summary
	awk '{print $3"\t""'$snp2'"}' $snp2|sed 's/.vcf//' >> SNPs.types.summary
else
	cat $snp1 > meme.suit.SNPs.vcf
fi
python $path/TF.network.analysis/1.motif.enrichment/1.allelic.motif.prediction/meme.suite.motif.prediction.py --vcf meme.suit.SNPs.vcf
echo "allelic motif prediction finished!"

#motif enrichment analysis
#enriched TFs pathway analysis
if [ ${#} == "7" ]; then
	python $path/TF.network.analysis/1.motif.enrichment/2.motif.enrichment/allelic.motif.enrichment.analysis.py --SNP SNPs.types.summary --motif meme-suit-motif-unique.txt
	rm -f -r SNPs.types.summary
	mkdir motif.enrichment
	mv Result.motif.fimo.txt ./motif.enrichment/
	mv meme-suit-motif-unique.txt ./motif.enrichment/
	mv *.allelic.motif ./motif.enrichment/
	mv *.allelic.motif.Fisher.significant.txt ./motif.enrichment/
	mv top10*enriched.pathways.pdf ./motif.enrichment/
	mv significant.enriched.TFs* ./motif.enrichment/
	rm meme.suit.SNPs.fa meme.suit.SNPs.vcf Rplots.pdf
	echo "motif enrichment analysis finished!"
else
	mkdir motif.enrichment
	mv Result.motif.fimo.txt ./motif.enrichment/
	mv meme-suit-motif-unique.txt ./motif.enrichment/
	rm meme.suit.SNPs.fa meme.suit.SNPs.vcf
	echo "No motif enrichment analysis becuase control SNPs not provided!"
fi
