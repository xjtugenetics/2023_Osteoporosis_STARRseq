#extract RP and read1 matched sequences
#read1:ACACGACGCTCTTCCGATCT(20-bp) RP:TCGAAGCGGCCGGCCGAATTCG(22-bp)
import time,copy,re,sys,os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--fq",type=str,required=True, help="extracted RP sequencing files<fastq format>")
args = parser.parse_args()
start=time.time()

def Reverse(m):
	m1=m[::-1]
	m2=[]
	for p in m1:
		if p=="A":
			m2.append("T")
		if p=="T":
			m2.append("A")
		if p=="G":
			m2.append("C")
		if p=="C":
			m2.append("G")
		if p=="N":
			m2.append("N")
	return "".join(m2)

def matchcount(m,m1):
	count=0
	for i in range(0,len(m1)):
		if m[i]==m1[i]:
			count+=1
	return count

def bestlen(FP):
	FP=list(FP)
	n=len(FP)
	w={}
	mm=["A","T","G","C","N"]
	for i in range(0,n):
		FP1=copy.deepcopy(FP)
		FP1.pop(i)
		w["".join(FP1)]=len(FP)-1
	for i in range(0,n):
		FP1=copy.deepcopy(FP)
		FP1.pop(i)
		for s in range(0,n-1):
			FP2=copy.deepcopy(FP1)
			mm4=FP2[s]
			for mm3 in mm:
				if mm3 !=mm4:
					FP2[s]=mm3
					if "".join(FP2) not in w:
						w["".join(FP2)]=len(FP)-1
					FP2[s]=mm4
	for i in range(0,n):
		FP1=copy.deepcopy(FP)
		FP1.pop(i)
		for j in range(0,n-1):
			FP2=copy.deepcopy(FP1)
			FP2.pop(j)
			if "".join(FP2) not in w:
				w["".join(FP2)]=len(FP)-2
	return w
w1=bestlen("ACACGACGCTCTTCCGATCT")  #read1
w2=bestlen("TCGAAGCGGCCGGCCGAATTCG") #RP

f2=open(args.fq)
f3=open(re.sub(".fq$","",args.fq)+".sequence","w")
i=0
for line in f2:
	if line[0]=="@":
		m=line.strip()
		i=2
	else:
		if i<=4:
			if i==2:
				sequence=line.strip()
				N111=matchcount(sequence[:20],"ACACGACGCTCTTCCGATCT") #read1
				if N111>=18:
					n11=20
					n111=N111
				elif N111<18 and sequence[:19] in w1:
					n11=19
					n111=w1.get(sequence[:19])
				elif N111<18 and sequence[:18] in w1:
					n11=18
					n111=w1.get(sequence[:18])
				else:
					n11=20
					n111=matchcount(sequence[:20],"ACACGACGCTCTTCCGATCT") #read1
				mmm="N"
				for sss in [13,12,11]:
					N22=n11+sss
					N222=matchcount(sequence[N22:(N22+22)],"TCGAAGCGGCCGGCCGAATTCG") #RP
					if N222>=20 and mmm!="Y":
						n22=22
						umi=sss
						n222=N222
						mmm="Y"
					elif sequence[N22:(N22+21)] in w2 and mmm!="Y":
						n22=21
						umi=sss
						n222=w2.get(sequence[N22:(N22+21)])
						mmm="Y"
					elif sequence[N22:(N22+20)] in w2 and mmm!="Y":
						n22=20
						umi=sss
						n222=w2.get(sequence[N22:(N22+20)])
						mmm="Y"
				if mmm!="Y":
					n22=22
					w3={}
					for sss in [13,12,11]:
						N22=n11+sss
						w3[matchcount(sequence[N22:(N22+22)],"TCGAAGCGGCCGGCCGAATTCG")]=sss #RP
					n222=max(w3.keys())
					umi=w3.get(n222)
				oligo=Reverse(line.strip()[(n22+umi+n11):])
				umiseq=Reverse(line.strip()[n11:(n11+umi)])
				f3.write(m+"\t"+umiseq+"\t"+str(n11)+"\t"+str(n111)+"\t"+str(umi)+"\t"+str(n22)+"\t"+str(n222)+"\t"+oligo+"\t+\t")
			elif i==4:
				f3.write(line.strip()[(n11+umi+n22):][::-1]+"\n")
			i+=1
f2.close()
f3.close()
end=time.time()
print "Step3....extract..RP..sequences...."+'%s seconds' % (end - start)
