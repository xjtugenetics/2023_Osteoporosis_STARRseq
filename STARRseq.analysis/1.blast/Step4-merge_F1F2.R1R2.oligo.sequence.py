#merge marge FP and RP sequences into 120-bp SNP sequences
#Dut to dulplicated sequences between FP and RP, we adopted two merge strategies (FP-FP.dup+RP or FP-RP.dup+RP)
import time,re,sys,os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--fq1",type=str,required=True, help="Uptream FP sequencing file<xxx.FP.unique.sequence>")
parser.add_argument("--fq2",type=str,required=True, help="Downtream RP.UMI.read1 sequencing file<xxx.RP.unique.sequence>")
args = parser.parse_args()
start=time.time()

f1=open(args.fq1)
w1={}
for line in f1:
	lines=re.split("[ \t][ \t]*",line.strip())
	if int(lines[3])>=19:
		w1[lines[0]]=lines[5:]
f1.close()

f2=open(args.fq2)
f3=open(re.sub(".FP.unique.sequence$","",args.fq1)+".merged1.fq","w")
f4=open(re.sub(".FP.unique.sequence$","",args.fq1)+".merged2.fq","w")
for line in f2:
	lines=re.split("[ \t][ \t]*",line.strip())
	if int(lines[4])>=18 and int(lines[7])>=20:
		if lines[0] in w1:
			m=w1[lines[0]]
			n1=len(m[0])-0
			n2=len(lines[8])-0
			if n1+n2>120:
				f3.write(lines[0]+":"+lines[2]+" "+lines[1]+"\n")
				f4.write(lines[0]+":"+lines[2]+" "+lines[1]+"\n")
				seq1=m[0]
				qc1=m[2]
				seq2=lines[8]
				qc2=lines[10]
				f3.write(seq1[:120]+"\n+\n"+qc1[-120:]+"\n")
				f4.write(seq1[:(120-n2)]+seq2+"\n+\n"+qc1[:(120-n2)]+qc2+"\n")
f3.close()
f4.close()
end=time.time()
print "Step4...extracted..SNP..sequences.."+'%s seconds' % (end - start)
