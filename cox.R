


library(survival)

setwd("D:\\cox")                       
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime/365

allGene=colnames(rt)[3:(ncol(rt)-2)]
allType=levels(factor(rt[,"CancerType"]))
pdf(file="forest.pdf", width = 15,height =10)
geneNum=length(allGene)
n=length(allType)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(1:(geneNum+1),nc=geneNum+1))
forestCol=rainbow(geneNum)

xlim = c(0,1)
par(mar=c(4.5,1,2,2))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text(0.5,n:1,allType,adj=0,cex=1.8)


i=0
for(gene in colnames(rt)[3:(ncol(rt)-2)]){
	i=i+1
	
	outTab=data.frame()
	for(CancerType in levels(factor(rt[,"CancerType"]))){
		rt1=rt[(rt[,"CancerType"]==CancerType),]
		cox=coxph(Surv(futime, fustat) ~ rt1[,gene], data = rt1)
		coxSummary = summary(cox)
		coxCoef=coxSummary$coefficients
	    coxP=coxCoef[,"Pr(>|z|)"]
		outTab=rbind(outTab,
		             cbind(cancer=CancerType,
		                   HR=coxSummary$conf.int[,"exp(coef)"],
		                   HR.95L=coxSummary$conf.int[,"lower .95"],
		                   HR.95H=coxSummary$conf.int[,"upper .95"],
				           pvalue=coxP) )
	}
	outFile=paste0(gene,".cox.txt")
	write.table(outTab,file=outFile,sep="\t",row.names=F,quote=F)   
		

	coxRT <- read.table(outFile,header=T,sep="\t",row.names=1,check.names=F)
	hr <- sprintf("%.3f",coxRT$"HR")
	hr[as.numeric(hr)<0.001]=0.001
	hrLow  <- sprintf("%.3f",coxRT$"HR.95L")
	hrLow[as.numeric(hrLow)<0.0001]=0.0001
	hrHigh <- sprintf("%.3f",coxRT$"HR.95H")
	hrHigh[as.numeric(hrHigh)>1000]=1000
	
	
	par(mar=c(4.5,1,2,2),mgp=c(2,0.5,0))
	LOGindex = 10 
	hrLow = log(as.numeric(hrLow),LOGindex)
	hrHigh = log(as.numeric(hrHigh),LOGindex)
	hr = log(as.numeric(hr),LOGindex)
	xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.01,col=forestCol[i],lwd=3)
	segments(log(1,LOGindex), 0, x1 = log(1,LOGindex), y1 = n)
	points(as.numeric(hr), n:1, pch = 15, col = forestCol[i], cex=2)
	a1 = axis(1,labels=F,tick=F)
	axis(1,a1,10^a1)
	text(log(1,LOGindex),n+1,gene,cex=2,xpd=T)
	if(i==round(geneNum/2,0)){
		text(log(1,LOGindex),-1.6,"Hazard Ratio",cex=2.5,xpd=T)
	}
}
dev.off()
