

library(limma)             
geneFile="gene.txt"        
setwd("D:\\GeneExp")    


gene=read.table(geneFile,sep="\t",header=F,check.names=F)
genelist=as.vector(gene[,1])
genelist=gsub(" ","",genelist)


files=dir()
files=grep("^symbol.",files,value=T)

outTab=data.frame()
for(i in files){
	
	CancerType=gsub("symbol\\.|\\.txt","",i)
	rt=read.table(i,sep="\t",header=T,check.names=F)
	
	
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	

	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	Type=ifelse(group==0,"Tumor","Normal")
	geneExp=t(data[genelist,])
	#geneExp=geneExp/data["TBP",]      ?
	outTab=rbind(outTab,cbind(geneExp,Type,CancerType))
}


out=cbind(ID=row.names(outTab),outTab)
write.table(out,file="GeneExp.txt",sep="\t",quote=F,row.names=F)
