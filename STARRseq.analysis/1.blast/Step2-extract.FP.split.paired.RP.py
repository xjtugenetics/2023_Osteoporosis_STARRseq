#extracted paired sequences(both fowared and reverse) and extract FP match sequences
#FP sequences:TAGATTGATCTAGAGCATGCA(20-bp)
import time,copy,re,sys,os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--fq1",type=str,required=True, help="splited sample sequencing file<fastq format>")
parser.add_argument("--fq2",type=str,required=True, help="Merged multi-samples sequencing file<fastq format>")
args = parser.parse_args()
FP="TAGATTGATCTAGAGCATGCA"  ##FP sequences
start=time.time()

f1=open(args.fq1)
w1={}
for line in f1:
	if line[0]=="@":
		lines=re.split("[ \t][ \t]*",line.strip())
		mm1=lines[0]
		mm2=lines[1].split(":")[0]
		if mm1 not in w1:
			w1[mm1]={}
		if mm2=="1":
			w1[mm1][mm2]="2"
		if mm2=="2":
			w1[mm1][mm2]="1"
f1.close()
end=time.time()		

w2={}
for read in w1:
	n=len(w1.get(read).keys())
	if n==1:
		w2[read]=w1[read].values()[0]
w1.clear()
end=time.time()		

f2=open(args.fq2)
w3={}
f4=open(re.sub(".fq$","",args.fq1)+".RP.unique.fq","w")
for line in f2:
	if line[0]=="@":
		lines=re.split("[ \t][ \t]*",line.strip())
		mm1=lines[0]
		mm2=lines[1].split(":")[0]
		match="N"
		if mm1 in w2:
			if mm2==w2[mm1]:
				f4.write(lines[0]+" "+lines[1]+"\n")
				match="Y"
				w3[mm1]=""
	else:
		if match=="Y":
			f4.write(line)
f2.close()
f4.close()
w2.clear()
end=time.time()		

f2=open(args.fq1)
f3=open(re.sub(".fq$","",args.fq1)+".FP.unique.fq","w")
for line in f2:
	if line[0]=="@":
		lines=re.split("[ \t][ \t]*",line.strip())
		mm1=lines[0]
		match="N"
		if mm1 in w3:
			f3.write(line)
			match="Y"
	else:
		if match=="Y":
			f3.write(line)
f2.close()
f3.close()

#extract.FP.match.sequences
FP=list(FP)
n=len(FP)
w={}
mm=["A","T","G","C","N"]
for i in range(0,n):
	FP1=copy.deepcopy(FP)
	FP1.pop(i)
	w["".join(FP1)]=20
	for i in range(0,n):
		FP1=copy.deepcopy(FP)
		FP1.pop(i)
		for ii in range(0,n-1):
			FP2=copy.deepcopy(FP1)
			mm4=FP2[ii]
			for mm3 in mm:
				if mm3 !=mm4:
					FP2[ii]=mm3
					if "".join(FP2) not in w:
						w["".join(FP2)]=19
					FP2[ii]=mm4
	for i in range(0,n):
		FP1=copy.deepcopy(FP)
		FP1.pop(i)
		for j in range(0,n-1):
			FP2=copy.deepcopy(FP1)
			FP2.pop(j)
			if "".join(FP2) not in w:
				w["".join(FP2)]=19
	FP1=copy.deepcopy(FP)

def matchcount(m,m1):
	count=0
	for i in range(0,len(m1)):
		if m[i]==m1[i]:
			count+=1
	return count

f2=open(re.sub(".fq$","",args.fq1)+".FP.unique.fq")
f3=open(re.sub(".fq$","",args.fq1)+".FP.unique.sequence","w")
i=0
for line in f2:
	if line[0]=="@":
		m=line.strip()
		i=2
	else:
		if i<=4:
			sequence=line.strip()
			if i==2:
				NN2=matchcount(sequence[:21],FP)
				if NN2>=19:
					nn1=21
					nn2=NN2
				elif NN2<18:
					if sequence[:20] in w:
						nn1=20
						nn2=w.get(sequence[:20])
					elif sequence[:19] in w:
						nn1=19
						nn2=w.get(sequence[:19])
					else:
						nn1=21
						nn2=matchcount(sequence[:21],FP)
				tail=sequence[(nn1-4):nn1]
				oligo=sequence[nn1:]
				f3.write(m+"\t"+str(nn1)+"\t"+str(nn2)+"\t"+tail+"\t"+oligo+"\t+\t")
			elif i==4:
				f3.write(sequence[nn1:]+"\n")
			i+=1
f2.close()
f3.close()
end=time.time()
print "Step2....extract..FP..match..fastq...."+'%s seconds' % (end - start)
