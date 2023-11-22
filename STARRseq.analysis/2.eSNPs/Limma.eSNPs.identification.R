#Rscript Limma.eSNPs.identification.R U2OS.filtered.expressed.10.reads
library("limma")
argv <- commandArgs(TRUE)
data<-read.table(argv[1],sep="\t",header=T)
rownames(data)<-data[,1]
data=data[,-1]

group_list <- factor(c(rep("Input",3),rep("Output",3)), levels = c("Input","Output")) #3 replications
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(data)
print(design)

v <- voom(data,design,plot=TRUE)
fit <- lmFit(v,design)
cont.matrix <- makeContrasts('Output-Input',levels = design)  #output (RNA) compared with input (DNA)
fit <- contrasts.fit(fit,cont.matrix)

fit2 <- eBayes(fit)
summary(decideTests(fit2))
plotSA(fit2)

tempOutput = topTable(fit2,coef=1,n=Inf)
DEG_voom = na.omit(tempOutput)

#output
write.table(DEG_voom,file="Limma.voom.adjusted.test.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(v$E,"Limma.voom.adjusted.expression.txt",sep="\t",row.names=T,quote=F,col.names=T)
