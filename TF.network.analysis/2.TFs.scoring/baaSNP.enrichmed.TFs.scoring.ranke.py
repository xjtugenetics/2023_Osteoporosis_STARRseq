#weightiness score for each baaSNP-TF-target gene connection pair,including Log2FC(STARR-seq),Log2(chromatin interaction reads) and allele-specific motif prediction score on baaSNP
import numpy as np
import math
import re,sys,os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--baaSNP",type=str,required=True, help="SNP function score <eg,STARR-seq assayed logFC>")
parser.add_argument("--motifscore",type=str,required=False, default="Result.motif.fimo.txt",help="motif prediction score output <default:Result.motif.fimo.txt>")
parser.add_argument("--motif",type=str,required=False, default="meme-suit-motif-unique.txt",help="allelic motif prediction summary output <default:meme-suit-motif-unique.txt>")
parser.add_argument("--TF",type=str,required=True,help="significant TFs list <enriched TFs for scoring>")
parser.add_argument("--HiC",type=str,required=True, help="SNP-gene regulation <eg,Hi-C interaction reads or eQTL beta>")
parser.add_argument("--output",type=str,required=True, help="scores for baaSNP-TF-target gene connection pairs <output>")
args = parser.parse_args()
args = parser.parse_args()
pw=re.sub("2023_Osteoporosis_STARRseq.*$","2023_Osteoporosis_STARRseq",os.getcwd())+"/TF.network.analysis/2.TFs.scoring"
os.environ["pw"]=pw

f1=open(args.baaSNP)
f1.readline()
w1={} #SNP:functional score
for line in f1:
	lines=line.strip().split("\t")
	if lines[1]!="NA":
		w1[lines[0]]=np.abs(float(lines[1]))
	else:
		w1[lines[0]]=1000000
f1.close()

f2=open(args.motif)
f2.readline()
w2={} #baaSNP:allelic motif
for line in f2:
	lines=line.strip().split("\t")
	if lines[0] in w1:
		w2[lines[0]]={}
		for motif in lines[1:]:
			if motif!=".":
				for TF1 in motif.split(", "):
					TF2=TF1.split("(")[0]
					allele=re.sub("\)","",TF1.split("(")[1])
					w2[lines[0]][TF2]=allele
f2.close()

f3=open(args.motifscore)
f3.readline()
w3={} #baaSNP:max motif score
for line in f3:
	lines=line.strip().split("\t")
	SNP=lines[3].split("-")[0]
	allele=lines[3].split("-")[1]
	TF=lines[1]
	score=lines[7]
	if SNP in w2 and TF in w2[SNP]:
		if allele==w2[SNP][TF]:
			if SNP+";"+TF not in w3:
				w3[SNP+";"+TF]=score
			else:
				if float(score)>float(w3[SNP+";"+TF]):
					w3[SNP+";"+TF]=score
f3.close()

f4=open(args.HiC)
f4.readline()
w4={} #baaSNP-gene:HiC.interaction
for line in f4:
	lines=line.strip().split("\t")
	if lines[0] not in w4:
		w4[lines[0]]={}
	w4[lines[0]][lines[1]]=lines[2]
f4.close()

f5=open(args.TF)
w5={} #selected TFs
for line in f5:
	w5[line.strip()]=""
f5.close()

#output
f5=open(args.output,"w")
f5.write("baaSNP\tTF\tgene\tSTARR.seq\tHiC.interaction\tMotif.score\n")
for snp in w2:
	for TF in w2[snp]:
		if TF in w5:
			for gene in w4[snp]:
				f5.write(snp+"\t"+TF+"\t"+gene+"\t"+str(w1[snp])+"\t"+str(w4[snp][gene])+"\t"+str(w3[snp+";"+TF])+"\n")
f5.close()

os.environ["motif"]=args.output
os.system("Rscript $pw/scr/TF.scoring.R $motif")
