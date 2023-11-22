#extract alleleic motif prediction comparasion
import re,sys,os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--SNP",type=str,required=True, help="list of baaSNPs & inactive.SNPs <SNP\tbaaSNPs/inactive.SNPs>")
parser.add_argument("--motif",type=str,required=True, help="allelic motif result for --SNP<meme-suit-motif-unique.txt>")
args = parser.parse_args()
pw=re.sub("2023_Osteoporosis_STARRseq.*$","2023_Osteoporosis_STARRseq",os.getcwd())+"/TF.network.analysis/1.motif.enrichment/2.motif.enrichment"
os.environ["pw"]=pw

f1=open(args.SNP)
w1={}
i=0
for line in f1:
	lines=line.strip().split("\t")
	if lines[1] not in w1:
		w1[lines[1]]={}
	w1[lines[1]][lines[0]]=""
	i+=1
	if i==1:
		compare=[lines[1]]
compare.append(lines[1])
f1.close()

f2=open(args.motif)
w2={} #motif:SNP
f2.readline()
for line in f2:
	lines=line.strip().split("\t")
	for TF1 in lines[1:]:
		if TF1!=".":
			TF2=TF1.split(", ")
			for TF3 in TF2:
				TF3=TF3.split("(")[0]
				if TF3 not in w2:
					w2[TF3]={}
				w2[TF3][lines[0]]=""
f2.close()

f3=open("_vs_".join(compare)+".allelic.motif","w")
N1=str(len(w1.get(compare[0]).keys()))
N2=str(len(w1.get(compare[1]).keys()))
for TF in w2:
	m1=0;m2=0
	for snp in w2[TF]:
		if snp in w1[compare[0]]:
			m1+=1
		if snp in w1[compare[1]]:
			m2+=1
	if m1>=5:  #keep TFs with at least 5 predicted binding baaSNPs
		f3.write(TF+"\t"+"_vs_".join(compare)+"\t"+str(m1)+"\t"+N1+"\t"+str(m2)+"\t"+N2+"\n")
f3.close()

os.environ["motif"]="_vs_".join(compare)+".allelic.motif"
os.system("Rscript $pw/script/motif.enrichment.test.R $motif")
os.system("awk '{print $1}' $motif\.Fisher.significant.txt|awk 'NR>1' > significant.enriched.TFs")

#functional anno for enriched TFs
os.system("Rscript $pw/script/TFs.pathway.enrich.plot.R significant.enriched.TFs")
f1=open("significant.enriched.TFs")
w={} #enriched.TFs
for line in f1:
	w[line.strip()]=""
f1.close()

files=os.listdir(pw+"/data")
for mm in files:
	f2=open(pw+"/data/"+mm)
	f3=open(re.sub("^TFs","significant.enriched.TFs",mm),"w")
	f3.write(f2.readline())
	for line in f2:
		lines=line.strip().split("\t")
		if lines[0] in w or lines[1] in w or lines[2] in w:
			f3.write(line)
	f2.close()
	f3.close()
print "motif enrichment and functional annotation finished!"
