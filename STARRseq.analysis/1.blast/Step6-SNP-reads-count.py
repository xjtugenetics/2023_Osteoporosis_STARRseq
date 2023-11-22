#extract blast read counts
#output file:Final.megred.SNPs.UMIs.count
import re,sys,os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--pair",type=str,required=True, help="mutation count files with 2 colums, sample.name\\tmerged1.mutation.count;merge2.mutation.count<xxx.mutation.count>")
args = parser.parse_args()

f1=open(args.pair)
w1={}
for line in f1:
	lines=line.strip().split("\t")
	w1[lines[0]]=lines[1]
f1.close()

w2={}
for sample in w1:
	w3={}
	for mm in w1.get(sample).split(";"):
		f2=open(mm)
		for line in f2:
			lines=line.strip().split("\t")
			if lines[3]=="120=":
				if lines[2] not in w3:
					w3[lines[2]]={}
				w3[lines[2]][lines[0]]=""
		f2.close()
	for m1 in w3:
		if m1 not in w2:
			w2[m1]={}
		w2[m1][sample]=len(w3.get(m1).keys())
	w3.clear()
	print sample+"...read..done"

sample=w1.keys()
sample.sort()

f3=open("Final.megred.SNPs.UMIs.count","w")
f3.write("SNP\t"+"\t".join(sample)+"\n")
for m in w2:
	mm1=[m]
	for ss in sample:
		mm1.append(str(w2.get(m).get(ss,"0")))
	f3.write("\t".join(mm1)+"\n")
f3.close()
print "Ste6..summary..read..counts..done"
