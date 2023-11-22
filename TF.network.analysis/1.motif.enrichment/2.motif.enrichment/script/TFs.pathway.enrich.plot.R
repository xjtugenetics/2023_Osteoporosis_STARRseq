#TF pathway enrichment analysis
#http://www.bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#introduction
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(Rgraphviz)
library(ggplot2)
argv <- commandArgs(TRUE)
data <- read.table(argv[1],header=F)
data2 = bitr(data[,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#GO enrichment
goMF <- enrichGO(gene=data2$"ENTREZID",OrgDb="org.Hs.eg.db",ont= "MF",pAdjustMethod = "fdr",pvalueCutoff=1,qvalueCutoff=1,minGSSize = 1,readable=TRUE)
goBP <- enrichGO(gene=data2$"ENTREZID",OrgDb="org.Hs.eg.db",ont= "BP",pAdjustMethod = "fdr",pvalueCutoff=1,qvalueCutoff=1,minGSSize = 1,readable=TRUE)

BP=as.data.frame(goBP)[,c(1:4,8,9,5)]
BP=BP[BP$"Count">=5,]
BP$"FDR"=p.adjust(BP$"pvalue",n=nrow(BP),method="fdr")
write.table(BP[BP$"FDR"<0.05,],file=paste(argv[1],".GO.BP.significant.enrichment.txt",sep=""),quote=F,sep="\t",row.names=F)

MF=as.data.frame(goMF)[,c(1:4,8,9,5)]
MF=MF[MF$"Count">=1,]
MF$"FDR"=p.adjust(MF$"pvalue",n=nrow(MF),method="fdr")
write.table(MF[MF$"FDR"<0.05,],file=paste(argv[1],".GO.MF.significant.enrichment.txt",sep=""),quote=F,sep="\t",row.names=F)

#plot GO.BP top 10 pathways(FDR<0.05)
BP=BP[BP$"FDR"<0.05,]
BP$"log10.FDR"=-1*log(BP$"FDR",10)
BP=droplevels(BP)
BP=BP[order(BP$"Count",decreasing=TRUE),]
if(nrow(BP)>=10){
	BP=BP[c(1:10),]
}
BP=BP[order(BP$"Count",decreasing=FALSE),]
BP=droplevels(BP)
BP$"Ratio"=1
for(i in 1:nrow(BP)){
	BP$"Ratio"[i]=as.numeric(strsplit(BP$"GeneRatio"[i],"/")[[1]][1])/as.numeric(strsplit(BP$"GeneRatio"[i],"/")[[1]][2])
}
BP$"Description"=factor(BP$"Description",levels=BP$"Description")
ggplot(BP,aes(y=round(Ratio,2),x=Description,color=log10.FDR))+geom_point(size=3*BP$"Count"/max(BP$"Count"))+labs(x="",y="GeneRatio")+scale_color_gradient(low="#CDC673",high="#A52A2A")+coord_flip()+ylim(0,round(max(BP$"Ratio")*1.2,2))+theme(panel.background=element_blank(),legend.title = element_text(size=6,colour="black"),legend.text=element_text(size=6,colour="black"),legend.key.height=unit(0.4,"cm"),legend.key.width=unit(0.3,"cm"),plot.title=element_text(size=6),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),axis.text.y=element_text(size=8),axis.text.x=element_text(size=7),axis.line.y=element_line(),axis.line.x=element_line())
m2=nrow(BP)/10*2.5
ggsave(file=paste("top10.GO.BP.enriched.pathways.pdf",sep=""),width=4.2,height=m2)

#plot GO.MF top 10 pathways(FDR<0.05)
MF=MF[MF$"FDR"<0.05,]
MF$"log10.FDR"=-1*log(MF$"FDR",10)
MF=droplevels(MF)
MF=MF[order(MF$"Count",decreasing=TRUE),]
if(nrow(MF)>=10){
	MF=MF[c(1:10),]
}
MF=MF[order(MF$"Count",decreasing=FALSE),]
MF=droplevels(MF)
MF$"Ratio"=1
for(i in 1:nrow(MF)){
	MF$"Ratio"[i]=as.numeric(strsplit(MF$"GeneRatio"[i],"/")[[1]][1])/as.numeric(strsplit(MF$"GeneRatio"[i],"/")[[1]][2])
}
MF$"Description"=factor(MF$"Description",levels=MF$"Description")
ggplot(MF,aes(y=round(Ratio,2),x=Description,color=log10.FDR))+geom_point(size=3*MF$"Count"/max(MF$"Count"))+labs(x="",y="GeneRatio")+scale_color_gradient(low="#CDC673",high="#A52A2A")+coord_flip()+ylim(0,round(max(MF$"Ratio")*1.2,2))+theme(panel.background=element_blank(),legend.title = element_text(size=6,colour="black"),legend.text=element_text(size=6,colour="black"),legend.key.height=unit(0.4,"cm"),legend.key.width=unit(0.3,"cm"),plot.title=element_text(size=6),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),axis.text.y=element_text(size=8),axis.text.x=element_text(size=7),axis.line.y=element_line(),axis.line.x=element_line())
m2=nrow(MF)/10*2.8
ggsave(file=paste("top10.GO.MF.enriched.pathways.pdf",sep=""),width=8,height=m2)
