##split merged fasteq into sample level according to barcodes sequences (less than 2 oligonucleotides lost/mismatch permmitted)
import time,re,sys,os,copy
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--fq",type=str,required=True, help="Merged multi-samples sequencing file<fastq format>")
parser.add_argument("--barcode",type=str,required=True, help="Barcode file for spliting different samples <sample\\tbarcode>")
args = parser.parse_args()
start=time.time()

#build wrong barcodes with 1 mis,1 lost,2 mis,1 lost+1 mis or 2 losts of oligonucleotides
def makewrong(m):
	w={}
	w["0.1.mis"]={}
	for i in range(0,len(m)):
		m1=copy.deepcopy(m)
		if i==0:
			w["0.1.mis"]["A"+m1[1:]]=""
			w["0.1.mis"]["T"+m1[1:]]=""
			w["0.1.mis"]["G"+m1[1:]]=""
			w["0.1.mis"]["C"+m1[1:]]=""
			w["0.1.mis"]["N"+m1[1:]]=""
		elif i<(len(m)-1):
			w["0.1.mis"][m1[:i]+"A"+m1[i+1:]]=""
			w["0.1.mis"][m1[:i]+"T"+m1[i+1:]]=""
			w["0.1.mis"][m1[:i]+"G"+m1[i+1:]]=""
			w["0.1.mis"][m1[:i]+"C"+m1[i+1:]]=""
			w["0.1.mis"][m1[:i]+"N"+m1[i+1:]]=""
		elif i==(len(m)-1):
			w["0.1.mis"][m1[:i]+"A"]=""
			w["0.1.mis"][m1[:i]+"T"]=""
			w["0.1.mis"][m1[:i]+"G"]=""
			w["0.1.mis"][m1[:i]+"C"]=""
			w["0.1.mis"][m1[:i]+"N"]=""
	w["1.lost"]={}
	for i in range(0,len(m)):
		m1=copy.deepcopy(m)
		m1=list(m1)
		m1.pop(i)
		m1="".join(m1)
		if m1 not in w["0.1.mis"]:
			w["1.lost"][m1]=""
	w["2.mis"]={}
	for p in w["0.1.mis"].keys():
		for i in range(0,len(p)):
			m1=copy.deepcopy(p)
			if i==0:
				if "A"+m1[1:] not in w["0.1.mis"]:
					w["2.mis"]["A"+m1[1:]]=""
				if "T"+m1[1:] not in w["0.1.mis"]:
					w["2.mis"]["T"+m1[1:]]=""
				if "G"+m1[1:] not in w["0.1.mis"]:
					w["2.mis"]["G"+m1[1:]]=""
				if "C"+m1[1:] not in w["0.1.mis"]:
					w["2.mis"]["C"+m1[1:]]=""
				if "N"+m1[1:] not in w["0.1.mis"]:
					w["2.mis"]["N"+m1[1:]]=""
			elif i<(len(m)-1):
				if m1[:i]+"A"+m1[i+1:] not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"A"+m1[i+1:]]=""
				if m1[:i]+"T"+m1[i+1:] not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"T"+m1[i+1:]]=""
				if m1[:i]+"G"+m1[i+1:] not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"G"+m1[i+1:]]=""
				if m1[:i]+"C"+m1[i+1:] not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"C"+m1[i+1:]]=""
				if m1[:i]+"N"+m1[i+1:] not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"N"+m1[i+1:]]=""
			elif i==(len(m)-1):
				if m1[:i]+"A" not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"A"]=""
				if m1[:i]+"T" not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"T"]=""
				if m1[:i]+"G" not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"G"]=""
				if m1[:i]+"C" not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"C"]=""
				if m1[:i]+"N" not in w["0.1.mis"]:
					w["2.mis"][m1[:i]+"N"]=""
	w["1.lost.1.mis"]={}
	for p in w["0.1.mis"].keys():
		for i in range(0,len(p)):
			m1=copy.deepcopy(p)
			m1=list(m1)
			m1.pop(i)
			m1="".join(m1)
			if m1 not in w["1.lost"] and m1 not in w["2.mis"] and m1 not in w["0.1.mis"]:
				w["1.lost.1.mis"][m1]=""
	w["2.lost"]={}
	for p in w["1.lost"].keys():
		for i in range(0,len(p)):
			m1=copy.deepcopy(p)
			m1=list(m1)
			m1.pop(i)
			m1="".join(m1)
			if m1 not in w["1.lost"] and m1 not in w["2.mis"] and m1 not in w["0.1.mis"] and m1 not in w["1.lost.1.mis"]:
				w["2.lost"][m1]=""
	return w

f1=open(args.barcode)
w={}
for line in f1:
	lines=line.strip().split("\t")
	w10=makewrong(lines[1])
	for p1 in w10:
		for p2 in w10[p1]:
			if p1 not in w:
				w[p1]={}
			w[p1][p2]=lines[0]
f1.close()

f1=open(args.fq)
w1={}
for line in f1:
	if line[0]=="@":
		i=1
		m=[line.strip()]
	else:
		i+=1
		if i==3:
			m.append(line.strip())
		if i==2:
			m1=line[:8]
			m2=line[:7]
			m3=line[:6]
			if m1 in w["0.1.mis"]:
				m.append(line.strip()[8:])
				length=8
				sample=w["0.1.mis"][m1]
			if m1 not in w["0.1.mis"] and m2 in w["1.lost"]:
				m.append(line.strip()[7:])
				length=7
				sample=w["1.lost"][m2]
			if m1 not in w["0.1.mis"] and m2 not in w["1.lost"] and m1 in w["2.mis"]:
				m.append(line.strip()[8:])
				length=8
				sample=w["2.mis"][m1]
			if m1 not in w["0.1.mis"] and m2 not in w["1.lost"] and m1 not in w["2.mis"] and m2 in w["1.lost.1.mis"]:
				m.append(line.strip()[7:])
				length=7
				sample=w["1.lost.1.mis"][m2]
			if m1 not in w["0.1.mis"] and m2 not in w["1.lost"] and m1 not in w["2.mis"] and m2 not in w["1.lost.1.mis"] and m3 in w["2.lost"]:
				m.append(line.strip()[6:])
				length=6
				sample=w["2.lost"][m3]
			if m1 not in w["0.1.mis"] and m2 not in w["1.lost"] and m1 not in w["2.mis"] and m2 not in w["1.lost.1.mis"] and m3 not in w["2.lost"]:
				length=0
		if i==4:
			if length>0:
				m.append(line.strip()[length:])
				if sample not in w1:
					w1[sample]={}
				w1[sample]["\n".join(m)]=""
f1.close()

for sample in w1.keys():
	f2=open(args.fq+"."+sample+".fq","w")
	for m in w1[sample]:
		f2.write(m+"\n")
	f2.close()
end=time.time()		
print "Step1....split..samples..fastq...."+'%s seconds' % (end - start)
