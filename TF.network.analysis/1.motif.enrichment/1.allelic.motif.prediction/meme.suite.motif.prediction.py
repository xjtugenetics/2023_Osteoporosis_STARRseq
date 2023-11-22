#prepare mysnps.vcf
from __future__ import division
import re,sys,os,time
import argparse
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
parser = argparse.ArgumentParser()
parser.add_argument("--vcf",type=str,required=True, help="Input SNP files for motif analysis<xxx.vcf format>")
parser.add_argument("--core", default=5,type=int,required=False, help="Prallel Threads for motif analysis. default 5")
args = parser.parse_args()
start=time.time()
pw=re.sub("2023_Osteoporosis_STARRseq.*$","2023_Osteoporosis_STARRseq",os.getcwd())+"/TF.network.analysis/1.motif.enrichment/1.allelic.motif.prediction"
os.environ["pw"]=pw

#extract region
f1=open(args.vcf)
w={}
for line in f1:
	lines=line.strip().split('\t')
	if lines[0] not in w:
		w[lines[0]]={}
	w[lines[0]][lines[2]]=[lines[1],lines[3],lines[4]]

#extract sequence
try:
	os.system("gzip -d $pw/files/motif_database/All-chromosome.fa.gz")
except:
	print ("All-chromosome.fa already uncompressed!")
record=SeqIO.to_dict(SeqIO.parse(pw+"/files/motif_database/All-chromosome.fa","fasta"))
f2=open(re.sub("\.vcf","",args.vcf)+".fa",'w')
for p in w:
	for q in w.get(p):
		n=w.get(p).get(q)
		n1=str(record["chr"+str(p)].seq[(int(n[0])-31):int(n[0])-1])
		n2=str(record["chr"+str(p)].seq[int(n[0]):(int(n[0])+30)])
		for d in n[1:]:
			f2.write('>'+q+'-'+d+'\n')
			f2.write(n1+d+n2+'\n')
	print "chr"+p+" is done"
f1.close()
f2.close()
print "Sequence extracted: 1/3" 

#motif analysis(Prallel)
#extract sequence info
f=open(re.sub("\.vcf",".fa",args.vcf))
w6={};linecount=0
for line in f:
	if line[0]==">":
		motif=line.strip()
	else:
		w6[motif]=line.strip()
	linecount+=1
f.close()

#adjust suitable Prallel core counts
if int(linecount)<=500:
	core=1
elif int(args.core)>=(int(linecount)/2):
	core=int(linecount/500)+1
else:
	core=int(args.core)
if core>=20:
	core=20

#split sequence file
motif=w6.keys()
w5={}
size=int((linecount/2)/core)+1
for i in range(0,core):
	w5[str(i)]={}
	for j in range(i*size,i*size+size):
		if j<len(motif):
			w5[str(i)][motif[j]]=""

#Read motif name annotation
#w9:{name:TF:URL}
f=open(pw+"/files/motif_database/fimo.motifs.anno")
w9={}
for line in f:
	lines=line.strip().split("\t")
	kk=[]
	for p in re.split("::|;",lines[1]):
		kk.append(p)
	w9[lines[0]]={}
	for p in kk:
		w9[lines[0]][p]=lines[2]
f.close()

#output_split_sequence
for p1 in w5:
	f=open("split."+str(p1),"w")
	for p2 in w5.get(p1):
		f.write(p2+"\n"+w6.get(p2)+"\n")
	f.close()
print "Sequence split finished: 2/3"

def task(n):
	db="split"+str(n)
	os.environ["db"]=db
	os.environ["sequence"]="split."+str(n)
	os.system("fimo --thresh 1e-4 --oc $db $pw/files/motif_database/fimo.meme $sequence")
	f2=open("./"+db+'/fimo.tsv')
	f3=open('fimo-'+db+'.txt.anno','w')
	f3.write('meme-db\tTF\tmotif_id\tSNP\t'+"\t".join(f2.readline().strip().split("\t")[3:])+'\tURL\n')
	for line in f2:
		if line.strip() and line[0]!="#":
			lines=line.strip().split('\t')
			del lines[1]
			line1="\t".join(lines)
			if lines[0] in w9:
				f3.write(re.sub("[$].*$","",lines[0])+'\t'+";".join(w9.get(lines[0]).keys())+"\t"+line1+"\t"+w9.get(lines[0]).values()[0].strip()+'\n')
	f2.close()
	f3.close()
	os.system("rm -f -r $db")

if  __name__ == '__main__':
	pool=Pool(core)
	pool.map(task,w5.keys())
	pool.close()
	pool.join()

os.system("head -1 fimo-split0.txt.anno>Result.motif.fimo.txt")
os.system("awk '$1!~/meme-db/' fimo-split*.txt.anno|awk '$5<=31 && $6>=31' >> Result.motif.fimo.txt")
os.system("rm -f split.* fimo-split*.txt.anno")
print "Motif analysis finished: 3/3"

#anno TF & make w1:rs.GC-db-TF
#w:motif.db
#w1:rs.GC-db-TF
w={}
w1={}
f2=open("Result.motif.fimo.txt")
f2.readline()
for line in f2:
	lines=line.strip().split('\t')
	snp=lines[3].split("-")[0]
	allele=lines[3].split("-")[1]
	database=lines[0]
	if snp not in w1:
		w1[snp]={}
	if database not in w1[snp]:
		w1[snp][database]={}
	if lines[1] not in w1[snp][database]:
		w1[snp][database][lines[1]]={}
	w1[snp][database][lines[1]][allele]=""
	w[database]=""
f2.close()

#filtered no-unique TFs
p1=w1.keys()
for pp1 in p1:
	p2=w1.get(pp1).keys()
	for pp2 in p2:
		p3=w1.get(pp1).get(pp2).keys()
		for pp3 in p3:
			if len(w1.get(pp1).get(pp2).get(pp3).keys())>=2:
				del w1[pp1][pp2][pp3]
		if not w1[pp1][pp2]:
			del w1[pp1][pp2]
	if not w1[pp1]:
		del w1[pp1]

#write unique results
f5=open("meme-suit-motif-unique.txt","w")
f5.write("SNP\t"+"\t".join(w.keys())+"\n")
for p in w1:
	mm=[]
	for pp in w:
		if pp not in w1[p]:
			mm.append(".")
		else:
			mmm=[]
			for k in w1.get(p).get(pp).keys():
				for kk in w1.get(p).get(pp).get(k).keys():
					mmm.append(k+"("+kk+")")
			mm.append(", ".join(sorted(mmm)))
	f5.write(p+"\t"+"\t".join(mm)+"\n")
f5.close()
end=time.time()
print '%s seconds' % (end - start)
