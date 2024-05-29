library(reshape2)
library(ggpubr)
expFile="sum.txt"           
data1=read.table(expFile, header=T, sep="\t", check.names=F,row.names = 1)
bioCol=c("#0066FF","#FF0000")
bioCol=bioCol[1:length(levels(factor(data1[,"cluster"])))]
p=ggboxplot(data1, x="type", y="NLRP1", fill="cluster", 
            ylab="Expression of NLRP1",
            xlab="",
            legend.title="Group",
            palette=bioCol)
p+stat_compare_means(aes(group=cluster),
                     method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
