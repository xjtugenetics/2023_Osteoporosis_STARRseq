#blast merged SNP sequences and summary blast results
import re,sys,os

#blast
os.system('ls *.merged1.fq|while read m;do bowtie2 -q $m -x ../bowtie2.index/example.SNPs.out -p 5 --local --xeq > $m\.out;done')
os.system("ls *.merged2.fq|while read m;do bowtie2 -q $m -x ../bowtie2.index/example.SNPs.out -p 5 --local --xeq > $m\.out;done")

blast=os.listdir(os.getcwd())
blast=[p.strip() for p in blast if re.match("^.*merge.*.out$",p.strip())]
for mm in blast:
	f1=open(mm)
	f2=open(mm+".mutation.count","w")
	f2.write("UMI\tID\tSNP\tBlast\t=\tM\tI\tD\tN\tS\tH\tP\tX\n")
	for line in f1:
		lines=re.split("[ \t][ \t]*",line.strip())
		if len(lines)>6 and line[0]!="@":
			if lines[5]!="*" and lines[5]!="":
				umi=lines[0].split(":")[-1]
				ID=":".join(lines[0].split(":")[:-1])
				SNPanno=lines[2].split("|")[0]
				SNP=lines[2].split("|")[0].split("_")[0]+"-"+lines[2].split("|")[-1]
				perfect=re.findall("[0-9][0-9]*=",lines[5])
				M=re.findall("[0-9][0-9]*M",lines[5])
				I=re.findall("[0-9][0-9]*I",lines[5])
				D=re.findall("[0-9][0-9]*D",lines[5])
				N=re.findall("[0-9][0-9]*N",lines[5])
				S=re.findall("[0-9][0-9]*S",lines[5])
				H=re.findall("[0-9][0-9]*H",lines[5])
				P=re.findall("[0-9][0-9]*P",lines[5])
				X=re.findall("[0-9][0-9]*X",lines[5])
				pe=0;m=0;i=0;d=0;n=0;ss=0;h=0;p=0;x=0
				if perfect:
					for s in perfect:
						pe+=int(re.sub("=","",s))
				if I:
					for s in I:
						i+=int(re.sub("I","",s))
				if D:
					for s in D:
						d+=int(re.sub("D","",s))
				if N:
					for s in N:
						n+=int(re.sub("N","",s))
				if S:
					for s in S:
						ss+=int(re.sub("S","",s))
				if M:
					for s in M:
						m+=int(re.sub("M","",s))
				if H:
					for s in H:
						h+=int(re.sub("H","",s))
				if P:
					for s in P:
						p+=int(re.sub("P","",s))
				if X:
					for s in X:
						x+=int(re.sub("X","",s))
				f2.write(umi+"\t"+ID+"\t"+SNPanno+"\t"+lines[5]+"\t"+str(pe)+"\t"+str(m)+"\t"+str(i)+"\t"+str(d)+"\t"+str(n)+"\t"+str(ss)+"\t"+str(h)+"\t"+str(p)+"\t"+str(x)+"\n")
	f1.close()
	f2.close()
print ("Step5..blast and extract results..finished")
