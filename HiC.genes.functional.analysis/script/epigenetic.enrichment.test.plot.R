#epigenetic enrichment test and plot
library("exact2x2")
library("ggplot2")
library("patchwork")
argv <- commandArgs(TRUE)
data<-read.table(argv[1],header=T,sep="\t")
result<-data.frame(1,2,3,4,5,6,7,8)
for(j in c(1:nrow(data))){
	x=matrix(c(data[j,3],data[j,5],(data[j,4]-data[j,3]),(data[j,6]-data[j,5])),nr=2,dimnames=list(c("positive","negative"),c("anno","No.anno")))
	T1=fisher.test(x,alternative="two.sided")
	T=exact2x2(x,tsmethod="minlike")
	result[j,]=c(as.character(data[j,1]),as.character(data[j,2]),paste(data[j,3],data[j,4],sep="/"),paste(data[j,5],data[j,6],sep="/"),as.numeric(T[3]),as.numeric(T[2]$conf.int)[1],as.numeric(T[2]$conf.int)[2],T1$p.value)
}
colnames(result)=c(colnames(data)[1],colnames(data)[2],colnames(data)[3],colnames(data)[5],"OR","confidence.interval.low","confidence.interval.high","Pvalue")

#Pvalue adjust
result$"FDR"=p.adjust(result$"Pvalue",n=nrow(result),method="fdr")
write.table(result,file=paste(argv[1],".epigentic.enrichment.Fisher.txt",sep=""),quote=F,sep="\t",row.names=F)

#adjust extremely large FDR
result$"Pvalue"=as.numeric(result$"Pvalue")
result$"log10.Pvalue"=-1*log(result$"Pvalue",10)
for (i in seq(1:nrow(result))){
	if (result$"log10.Pvalue"[i]>10){
		result$"log10.Pvalue"[i]=10
	}
}

#plot
for (m in levels(factor(result$"epigenetic"))){
	data1=result[result$"epigenetic"==m,]
	data1=droplevels(data1)
	data1=data1[order(data1$"OR",decreasing=FALSE),]
	data1$"anno"=factor(data1$"anno",levels=data1$"anno")
	print(paste(m,paste(nrow(data1[data1$"Pvalue"<0.05,]),nrow(data1),sep="/"),"at Pvalue<0.05",sep=" "))
	data1$"OR"=as.numeric(data1$"OR")
	data1$"confidence.interval.low"=as.numeric(data1$"confidence.interval.low")
	data1$"confidence.interval.high"=as.numeric(data1$"confidence.interval.high")
	data1$"log10.Pvalue"=as.numeric(data1$"log10.Pvalue")
	ggplot(data1,aes(y=round(OR,4),x=anno,color=log10.Pvalue))+geom_point(size=0.7)+geom_errorbar(aes(ymin=round(confidence.interval.low,4), ymax=round(confidence.interval.high,4),width=0.25))+labs(x="",y="OR")+scale_color_gradient(low="#CDC673",high="#A52A2A")+coord_flip()+ylim(0,round(ceiling(max(data1$"confidence.interval.high")),1))+theme(panel.background=element_blank(),panel.border=element_blank(),panel.grid=element_blank(),legend.title = element_text(size=6,colour="black"),legend.text=element_text(size=6,colour="black"),legend.key.height=unit(0.4,"cm"),legend.key.width=unit(0.3,"cm"),plot.title=element_text(size=6),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),axis.text.y=element_text(size=8),axis.line.y=element_line(),axis.line.x=element_line())+geom_hline(yintercept=1,linetype="dotted",color ="#C0C0C0")
	m2=nrow(data1)/15*3
	if(m2<=1.5){
		m2=1.5
	}
	ggsave(file=paste(m,".plot.pdf",sep=""),width=3.5,height=m2)
}
print ("Epigenetic.enrichment & plot task done!")
