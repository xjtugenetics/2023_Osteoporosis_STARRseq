#output MPRAnalyze input files from ../1.blast/output(UMI reads:SNP\tInput1\tInput2\tInput3\tOutput1\tOutput2\tOutput3)
import re,sys,os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--eSNP",type=str,required=True, help="eSNPs list <vcf format>")
parser.add_argument("--reads",type=str,required=True, help="All SNPs UMI reads")
args = parser.parse_args()

f1=open(args.eSNP)
w1={}
for line in f1:
	lines=line.strip().split("\t")
	w1[lines[2]]=[lines[3],lines[4]]
f1.close()

f2=open(args.reads)
n=(len(f2.readline().strip().split("\t"))-1)/2 #how many replicates
w2={} #inpt
w3={} #output
for line in f2:
	lines=line.strip().split("\t")
	if "-" in lines[0]:
		snp=lines[0].split("-")[0]
		allele=lines[0].split("-")[1]
		if snp in w1:
			if snp not in w2:
				w2[snp]={}
				w3[snp]={}
			w2[snp][allele]=lines[1:(n+1)]
			w3[snp][allele]=lines[(n+1):]
f2.close()

for snp in w2:
	alleles=w2.get(snp).keys()
	if len(alleles)==1:
		for allele in w1.get(snp):
			if allele!=alleles[0]:
				w2[snp][allele]=["0","0","0"]
				w3[snp][allele]=["0","0","0"]
mm=["SNP"] #head
for i in range(1,n+1):
	mm.append("allele1_sample"+str(i))
for i in range(1,n+1):
	mm.append("allele2_sample"+str(i))

f3=open("input.DNA.eSNPs.reads.txt","w")
f4=open("output.RNA.eSNPs.reads.txt","w")
f3.write("\t".join(mm)+"\n")
f4.write("\t".join(mm)+"\n")
for snp in w1:
	alleles=w1.get(snp)
	if len(alleles)==2:
		f3.write(snp+"|"+alleles[0]+"/"+alleles[1]+"\t"+"\t".join(w2.get(snp).get(alleles[0]))+"\t"+"\t".join(w2.get(snp).get(alleles[1]))+"\n")
		f4.write(snp+"|"+alleles[0]+"/"+alleles[1]+"\t"+"\t".join(w3.get(snp).get(alleles[0]))+"\t"+"\t".join(w3.get(snp).get(alleles[1]))+"\n")
f3.close()
f4.close()
