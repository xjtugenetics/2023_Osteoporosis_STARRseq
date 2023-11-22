#comparision of epigenetic annotation between input positive/negative SNPs
import re,sys,os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--SNP1",type=str,required=True, help="input positive SNPs <vcf format>")
parser.add_argument("--SNP2",type=str,required=True, help="input backgound SNPs <vcf format,necessary for epigenetic & motif enrichment analysis>")
parser.add_argument("--region",type=str,required=False, default="60",help="expanded region for epigenetic annotation <xxx bp>")
parser.add_argument("--out",type=str,required=True,help="output file")
args = parser.parse_args()

os.environ["SNP1"]=args.SNP1
os.environ["SNP2"]=args.SNP2
os.environ["region"]=args.region
os.environ["path"]=re.sub("\/[^/]*$","",os.getcwd())

#make input bed file
os.system("awk 'BEGIN{OFS=\"\t\"}{print \"chr\"$1,$2,$2,$3\";\"\"'$SNP1'\"}' $SNP1|sed 's/.vcf$//' > compared.SNPs.bed")
os.system("awk 'BEGIN{OFS=\"\t\"}{print \"chr\"$1,$2,$2,$3\";\"\"'$SNP2'\"}' $SNP2|sed 's/.vcf$//'>> compared.SNPs.bed")
w1={} #positive/negative SNPs<eg:HiC_SNPs VS no.HiC_SNPs>
f1=open("compared.SNPs.bed")
for line in f1:
	lines=line.strip().split("\t")
	SNP=lines[3].split(";")[0]
	group=lines[3].split(";")[1]
	if group not in w1:
		w1[group]={}
	w1[group][SNP]=""
f1.close()

w2={} #epigentic:group:SNPs
os.system("awk 'BEGIN{OFS=\"\t\"}{print $1,$2-\"'$region'\",$3+\"'$region'\",$4}' compared.SNPs.bed> compared.SNPs.expaned.region.bed")
os.system("rm -f compared.SNPs.bed")
f2=os.popen("zcat $path/HiC.genes.functional.analysis/data/Combined.osteoblast.epigenetic.anno.bed.gz|bedtools intersect -a compared.SNPs.expaned.region.bed -b - -wo")
for line in f2:
	lines=line.strip().split("\t")
	SNP=lines[3].split(";")[0]
	group=lines[3].split(";")[1]
	anno=lines[7]
	if anno not in w2:
		w2[anno]={}
	if group not in w2[anno]:
		w2[anno][group]={}
	w2[anno][group][SNP]=""
f2.close()
os.system("rm -f compared.SNPs.expaned.region.bed")

#output
f3=open(args.out,"w")
group=[re.sub(".vcf$","",args.SNP1),re.sub(".vcf$","",args.SNP2)]
N1=str(len(w1.get(group[0]).keys()))
N2=str(len(w1.get(group[1]).keys()))
f3.write("epigenetic\tanno\t"+group[0]+".anno\t"+group[0]+"\t"+group[1]+".anno\t"+group[1]+"\n")
annos=w2.keys()
annos.sort()
for anno in annos:
	if group[0] in w2[anno]:
		n1=str(len(w2.get(anno).get(group[0]).keys()))
	elif group[0] not in w2[anno]:
		n1="0"
	if group[1] in w2[anno]:
		n2=str(len(w2.get(anno).get(group[1]).keys()))
	elif group[1] not in w2[anno]:
		n2="0"
	if int(n1)>=5 or int(n2)>=5: #keep annotation with at least 5 SNPs
		f3.write(re.sub(";","\t",anno)+"\t"+n1+"\t"+N1+"\t"+n2+"\t"+N2+"\n")
f3.close()

os.environ["out"]=args.out
os.system("Rscript $path/HiC.genes.functional.analysis/script/epigenetic.enrichment.test.plot.R $out")
os.system("rm -f $out")
