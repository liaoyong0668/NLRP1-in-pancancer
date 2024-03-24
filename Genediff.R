


library(ggpubr)              
expFile="GeneExp.txt"      
setwd("D:\\diff")                    
data=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)    


Normal=data[data[,"Type"]=="Normal",]
NormalNum=table(Normal[,"CancerType"])
NormalNum=NormalNum[NormalNum>=5]
NormalCacner=names(NormalNum)
data=data[which(data[,"CancerType"] %in% NormalCacner),]


genelist=colnames(data[,(1:(ncol(data)-2))])
for(gene in genelist){
	rt1=data[,c(gene,"Type","CancerType")]
	colnames(rt1)[1]="expression"
	
	p=ggboxplot(rt1, x="CancerType", y="expression", color = "Type", 
	     ylab=paste0(gene," expression"),
	     xlab="Cancer type",
	     palette = c("blue","red") )
	p=p+rotate_x_text(60)
	p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")
	pdf(file=paste0(gene,".diff.pdf"),width=8,height=5)    
	print(p1)
	dev.off()
}
