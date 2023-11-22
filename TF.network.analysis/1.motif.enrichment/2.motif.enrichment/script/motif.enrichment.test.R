#Rscript Step2.enrichment.test.R baaSNPs_vs_inactive.SNPs.allelic.motif
library("exact2x2")
argv <- commandArgs(TRUE)
data<-read.table(argv[1],header=F,sep="\t")
m1=strsplit(as.character(data[1,2]),split="_vs_")[[1]][1] #input SNPs
m2=strsplit(as.character(data[1,2]),split="_vs_")[[1]][2] #backgound SNPs
result<-data.frame(Anno=1,input.anno=1,control.anno=1,OR=1,low=1,high=1,Pvalue=1)
for(j in c(1:nrow(data))){
	x=matrix(c(data[j,3],data[j,5],(data[j,4]-data[j,3]),(data[j,6]-data[j,5])),nr=2,dimnames=list(c("positive","negative"),c("anno","Not.anno")))
	T1=fisher.test(x,alternative="greater")
	T=exact2x2(x,tsmethod="minlike")
	result[j,]=c(as.character(data[j,1]),paste(data[j,3],data[j,4],sep="/"),paste(data[j,5],data[j,6],sep="/"),as.numeric(T[3]),as.numeric(T[2]$conf.int)[1],as.numeric(T[2]$conf.int)[2],T1$p.value)
}
colnames(result)=c("TF",paste(m1,".anno",sep=""),paste(m2,".anno",sep=""),"OR","confidence.interval.low","confidence.interval.high","Pvalue")

#keep significant results
write.table(result[result$"Pvalue"<0.05,],file=paste(argv[1],".Fisher.significant.txt",sep=""),quote=F,sep="\t",row.names=F)
warnings()
