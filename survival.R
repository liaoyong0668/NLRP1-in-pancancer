


library(survival)
library(survminer)

inputFile="expTime.txt"      
pFilter=0.05                 
col=c("red","blue")          
setwd("D:\\survival")                
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime/30    


outTab=data.frame()
for(gene in colnames(rt)[3:(ncol(rt)-2)]){
	
	for(CancerType in levels(factor(rt[,"CancerType"]))){
		rt1=rt[(rt[,"CancerType"]==CancerType),]
		group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
		diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
		pValue=1-pchisq(diff$chisq,df=1)
		if(pValue<pFilter){
			outVector=cbind(gene,CancerType,pValue)
			outTab=rbind(outTab,outVector)
			if(pValue<0.001){
				pValue="p<0.001"
			}else{
				pValue=paste0("p=",sprintf("%.03f",pValue))
			}
			fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
			
			surPlot=ggsurvplot(fit, 
					    data=rt1,
					    title=paste0("Cancer: ",CancerType),
					    pval=pValue,
					    pval.size=6,
					    legend.labs=c("high","low"),
					    legend.title=paste0(gene," levels"),
					    font.legend=12,
					    xlab="Time(years)",
					    palette=col,
					    break.time.by = 2,
					    conf.int=T,
					    fontsize=4,
					    risk.table=TRUE,
					    ylab="Overall survival",
					    risk.table.title="",
					    risk.table.height=.25)
			pdf(file=paste0(gene,"_",CancerType,".pdf"),onefile = FALSE,
					    width = 6,         #ͼƬ?Ŀ???
					    height =5)         #ͼƬ?ĸ߶?
			print(surPlot)
			dev.off()
		}
	}
}

write.table(outTab,file="survival.result.xls",sep="\t",row.names=F,quote=F)
