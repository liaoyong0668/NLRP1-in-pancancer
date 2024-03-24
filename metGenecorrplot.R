


library(corrplot)                
inputFile="panGeneExp.txt"       
setwd("D:\\met_corrplot")    

rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)    
rt=rt[(rt[,"Type"]=="Tumor"),]       
data=rt[,(1:(ncol(rt)-2))]           
geneNum=ncol(data)                   


pdf("corrplot.pdf",height=6,width=6)             
par(oma=c(0.5,0.5,0.5,1))
M=cor(data)
corrplot(M, order = "AOE", type = "upper", tl.pos = "lt")
corrplot(M, add = TRUE, type = "lower", method = "number", order = "AOE",
         col = "black", diag = FALSE, tl.pos = "n", cl.pos = "n")
dev.off()

