#Extract Hi-C interacted genes for input SNPs && annotated known bone-relevant functions or potential causal regulatory effect for osteoporosis on Hi-C genes
import re,sys,os,gzip
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--SNP",type=str,required=True, help="input SNPs <vcf format>")
args = parser.parse_args()
path=re.sub("\/[^/]*$","",os.getcwd())
os.environ["path"]=path
os.environ["SNP"]=re.sub(".vcf$","",args.SNP)

#anno.HiC.interacted.genes
os.system("awk 'BEGIN{OFS=\"\t\"}{print \"chr\"$1,$2,$2,$3}' $SNP\.vcf > input.SNPs.bed")
os.system("awk 'BEGIN{OFS=\"\t\"}{if($2>1000000)print $1,$2-1000000,$3+1000000,$4;if($2<=1000000)print $1,\"1\",$3+1000000,$4}' input.SNPs.bed|bedtools intersect -a - -b $path/HiC.genes.functional.analysis/data/gencode.v19.transcript_TSS1KB_promoter.bed -wo|awk 'BEGIN{OFS=\"\t\"}{print $1,$3-1000000-1,$3-1000000+1,$5,$6,$7,$4\";\"$8}'|awk 'BEGIN{OFS=\"\t\"}{if($2>0)print $0;if($2<0)print $1,\"1\",$3,$4,$5,$6,$7}'|awk 'BEGIN{OFS=\"\t\"}{if($5>0)print $0;if($5<0)print $1,$2,$3,$4,\"1\",$6,$7}'|bedtools pairtopair -a - -b $path/HiC.genes.functional.analysis/data/hMSC.induced.osteoblast.significant.HiC.bedpe -type both|awk 'BEGIN{OFS=\"\t\"}{split($14,m,\";\");if($9<$12)print $7,$8,$9,$10,$11,$12,$13,m[1],m[2],m[3];if($9>=$12)print $7,$11,$12,$13,$8,$9,$10,m[1],m[2],m[3]}'|sort|uniq|sed '1i SNP.gene.pairs\tchr(1)\tstart(1)\tend(1)\tchr(2)\tstart(2)\tend(2)\tInteraction.Reads\t-Log10P\tFDR' > SNPs.osteoblast.HiC.interacted.genes")
os.system("rm -f input.SNPs.bed")

#gene.anno
f2=open("SNPs.osteoblast.HiC.interacted.genes")
f2.readline() #ignore head
w2={} #HiC interacted gene
for line in f2:
	lines=line.strip().split("\t")
	gene=lines[0].split(";")[1]
	w2[gene]=""
f2.close()

m1=["All.genome.genes.smultixcan.bonferroni.significant","All.osteoporosis.GWAS.GTEx.V8.fastenloc.RCP.0.01","IMPC.release17.significant.skeletal.related.txt"]
m2=["smultixcan","fastenloc","IMPC"]
w1={}
for m in m2:
	w1[m]={}
for i in range(0,len(m1)):
	f1=gzip.open(path+"/HiC.genes.functional.analysis/data/"+m1[i]+".gz")
	f2=open(re.sub(".vcf$","",args.SNP)+"."+m2[i]+".HiC.genes.anno.txt","w")
	f2.write(f1.readline())
	for line in f1:
		lines=line.strip().split("\t")
		if i==1:
			if lines[0] in w2:
				w1[m2[i]][lines[0]]=""
				f2.write(line)
		if i==0:
			if lines[2] in w2:
				w1[m2[i]][lines[2]]=""
				f2.write(line)
		else:
			if lines[2] in w2:
				if lines[2] not in w1[m2[i]]:
					w1[m2[i]][lines[2]]={}
				w1[m2[i]][lines[2]][lines[-1]]=""
				f2.write(line)
	f1.close()
	f2.close()

#summary annotation
f3=open(re.sub(".vcf$","",args.SNP)+".HiC.genes.function.summary","w")
f3.write("gene\tskeletal phenotypic abnormalities in gene knockout mouse (IMPC)\tOsteoporosis-associated genes\n")
f4=open("HiC.interacted.genes","w")
for gene in w2:
	if gene in w1["IMPC"]:
		mm=[", ".join(w1.get("IMPC").get(gene).keys())]
	elif gene not in w1["IMPC"]:
		mm=["-"]
	mm1=[]
	if gene in w1["smultixcan"]:
		mm1.append("smultixcan")
	if gene in w1["fastenloc"]:
		mm1.append("fastenloc")
	if not mm1:
		mm.append("-")
	if mm1:
		mm.append(", ".join(mm1))
	if mm!=["-","-"]:
		f3.write(gene+"\t"+"\t".join(mm)+"\n")
	f4.write(gene+"\n")
f3.close()
f4.close()
os.system("Rscript $path/HiC.genes.functional.analysis/script/Gene.pathway.enrich.plot.R HiC.interacted.genes")
os.system("rm -f HiC.interacted.genes Rplots.pdf")
