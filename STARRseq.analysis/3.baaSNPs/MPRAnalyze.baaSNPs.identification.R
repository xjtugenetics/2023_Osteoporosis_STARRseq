#prepare three files: output(RNA)/input(DNA)/group.compare
library("MPRAnalyze")
argv <- commandArgs(TRUE)

agg_rna=read.table(argv[1],header=T,sep="\t") #read RNA(output) reads
agg_dna=read.table(argv[2],header=T,sep="\t") #read DNA(input) reads
rownames(agg_rna)=agg_rna[,1]
rownames(agg_dna)=agg_dna[,1]
agg_rna=agg_rna[,-1]
agg_dna=agg_dna[,-1]
agg_rna=as.matrix(agg_rna)
agg_dna=as.matrix(agg_dna)

group=read.table("group.compare.txt",header=T,sep="\t") #read groups comparing info, allele2 compared with allele1
rownames(group)=group[,1]
group=group[,-1]
obj <- MpraObject(dnaCounts=agg_dna,rnaCounts=agg_rna,colAnnot=group,controls = NA)
obj <- estimateDepthFactors(obj,lib.factor = "batch", which.lib = "both")
obj <- analyzeComparative(obj, dnaDesign = ~ batch + condition,rnaDesign = ~ condition,correctControls=F,reducedDesign=~1)
results <- testLrt(obj)
write.table(results,file="MPRAnalyze.baaSNPs.analysis.txt",quote=F,row.names = T,col.names = T,sep="\t")
