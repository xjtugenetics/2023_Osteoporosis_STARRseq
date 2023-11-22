#Rscipt Step2.TF.scoring.R Bone.baaSNPs.enriched.TFs.gene.paris.score
argv <- commandArgs(TRUE)
data=read.table(argv[1],header=T)
data$"STARR.seq.rank"=rank(data[,4])/nrow(data)
data$"HiC.interaction.rank"=rank(data[,5])/nrow(data)
data$"Motif.score.rank"=rank(data[,6])/nrow(data)
data$"pair.weightiness.score"=data$"STARR.seq.rank"*data$"HiC.interaction.rank"*data$"Motif.score.rank"
data$"cumulative.TF.weightiness.score"="NA"
data1=aggregate(data[,7]*data[,8]*data[,9],list(data[,2]),sum)
data1=data1[order(data1[,2],decreasing=T),]
colnames(data1)=c("TF","score")
for(i in c(1:nrow(data))){
	data[i,11]=data1[data1[,1]==data[i,2],2]
}
write.table(data[order(data[,11],data[,1],decreasing=T),c(1:3,7:11)],file=paste(argv[1],".normalized.rank",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
