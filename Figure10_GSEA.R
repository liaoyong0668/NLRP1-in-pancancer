library(clusterProfiler)
library(limma)
library(ggplot2)
library(data.table)
library(ggpubr)
library(GSVA)

samAnno <- read.table("easy_input_sample_annotation.txt", sep = "\t",header = T, check.names = F)
expr <- fread("geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
expr$gene <- gene
expr <- expr[!duplicated(expr$gene),]
rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]

twoclasslimma <- function(subtype   = NULL,
                          featmat   = NULL,
                          treatVar  = NULL,
                          ctrlVar   = NULL,
                          prefix    = NULL,
                          overwt    = FALSE,
                          sort.p    = TRUE,
                          verbose   = TRUE,
                          res.path  = getwd()) {
  
  library(limma)
  if(!is.element("condition", colnames(subtype))) {
    stop("argument of subtype must contain a column named with 'condition'!")
  }
  # Create comparison list for differential analysis between two classes.
  createList  <- function(subtype = NULL) {
    
    tumorsam <- rownames(subtype)
    sampleList <- list()
    treatsamList <- list()
    treatnameList <- c()
    ctrlnameList <- c()
    
    sampleList[[1]] <- tumorsam
    treatsamList[[1]] <- intersect(tumorsam, rownames(subtype[which(subtype$condition == treatVar),,drop = F]))
    treatnameList[1] <- treatVar
    ctrlnameList[1] <- ctrlVar
    
    return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
  }
  
  complist <- createList(subtype = subtype)
  
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(featmat)
  
  # log transformation
  if(max(featmat) < 25 | (max(featmat) >= 25 & min(featmat) < 0)) {
    message("--feature matrix seems to have been standardised (z-score or log transformation), no more action will be performed.")
    gset <- featmat
  }
  if(max(featmat) >= 25 & min(featmat) >= 0){
    message("--log2 transformation done for feature matrix.")
    gset <- log2(featmat + 1)
  }
  
  options(warn = 1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    
    compname <- paste(treatname, "_vs_", ctrlname, sep="")
    tmp <- rep("others", times = length(allsamples))
    names(tmp) <- allsamples
    tmp[samples] <- "control"
    tmp[treatsam] <- "treatment"
    
    if(!is.null(prefix)) {
      outfile <- file.path(res.path, paste(prefix, "_limma_test_result.", compname, ".txt", sep = ""))
    } else {
      outfile <- file.path(res.path, paste("limma_test_result.", compname, ".txt", sep = ""))
    }
    
    if (file.exists(outfile) & (overwt == FALSE)) {
      cat(paste0("limma of ",compname, " exists and skipped...\n"))
      next
    }
    
    pd <- data.frame(Samples = names(tmp),
                     Group = as.character(tmp),
                     stringsAsFactors = FALSE)
    
    design <-model.matrix(~ -1 + factor(pd$Group, levels = c("treatment","control")))
    colnames(design) <- c("treatment","control")
    
    fit <- limma::lmFit(gset, design = design);
    contrastsMatrix <- limma::makeContrasts(treatment - control, levels = c("treatment", "control"))
    fit2 <- limma::contrasts.fit(fit, contrasts = contrastsMatrix)
    fit2 <- limma::eBayes(fit2, 0.01)
    resData <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", number = 100000)
    resData <- as.data.frame(subset(resData, select=c("logFC","t","B","P.Value","adj.P.Val")))
    resData$id <- rownames(resData)
    colnames(resData) <- c("log2fc","t","B","pvalue","padj","id")
    resData$fc <- 2^resData$log2fc
    
    if(sort.p) {
      resData <- resData[order(resData$padj),]
    } else {
      resData <- as.data.frame(resData)
    }
    if(verbose) {
      resData <- resData[,c("id","fc","log2fc","pvalue","padj")]
    } else {
      resData <- resData[,c("id","fc","log2fc","t","B","pvalue","padj")]
    }
    write.table(resData, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
    cat(paste0("limma of ",compname, " done...\n"))
  }
  options(warn = 0)
}
es <- log2(expr[rownames(expr) == "NLRP1",] + 0.01)
msigdb.hallmark <- read.gmt("h.all.v7.2.symbols.gmt")
pct <- 0.3 
gseaTab <- NULL
tumors <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD", #6
            "DLBC","ESCA","GBM","HNSC","KICH","KIRC",
            "KIRP","LAML","LGG","LIHC","LUAD","LUSC",
            "MESO","OV","PAAD","PCPG","PRAD","READ",
            "SARC","SKCM","STAD","TGCT","THCA","THYM",
            "UCEC","UCS","UVM")
for (i in tumors) {
  message("--",i,"...")
  sam <- samAnno[which(samAnno$`cancer type` == i),"simple_barcode"]
  comsam <- intersect(colnames(es), sam) 
  tumsam <- comsam[substr(comsam,14,14) == "1"] 
  es_subset <- as.data.frame(t(es[,tumsam]))
  es_subset <- es_subset[order(es_subset$NLRP1,decreasing = T),,drop = F] 
  high.es <- rownames(es_subset[1:(nrow(es_subset) * pct),,drop = F]) # 取前30%为高???
  low.es <- rownames(es_subset[nrow(es_subset):(nrow(es_subset) - nrow(es_subset) * pct + 1),,drop = F])
  
  subt <- data.frame(condition = rep(c("high","low"),c(length(high.es),length(low.es))),
                     row.names = c(high.es, low.es),
                     stringsAsFactors = F)
  gset <- log2(na.omit(expr[,rownames(subt)]) + 1)
  twoclasslimma(subtype  = subt, 
                featmat  = gset, 
                treatVar = "high", 
                ctrlVar  = "low", 
                prefix   = paste0("TCGA_enrichment_score_",i), 
                overwt   = T,
                sort.p   = F, 
                verbose  = TRUE, 
                res.path = ".") 
  res <- read.table( paste0("TCGA_enrichment_score_",i,"_limma_test_result.high_vs_low.txt"), sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  
  res <- res[order(res$log2fc, decreasing = T),]
  glist <- res$log2fc; names(glist) <- rownames(res)
  
  set.seed(20211114)
  gsea <- GSEA(geneList = glist,
               pvalueCutoff = 1, 
               seed = TRUE,
               TERM2GENE = msigdb.hallmark)
  gc()
  gsea.df <- as.data.frame(gsea) 
  write.table(gsea.df,file = "output_gsea between high and low group of enrichment score in pancancer.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
  
  gseaTab <- rbind.data.frame(gseaTab,
                              data.frame(term = gsea.df$ID,
                                         NES = gsea.df$NES,
                                         FDR = gsea.df$p.adjust,
                                         tumor = i,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
}

darkblue <- "#303B7F"
darkred <- "#D51113"
yellow <- "#EECA1F"
tmp <- gseaTab
tmp$term <- gsub("HALLMARK_","",tmp$term)
my_palette <- colorRampPalette(c(darkblue,yellow,darkred), alpha=TRUE)(n=128)
ggplot(tmp, aes(x=tumor,y=term)) +
  geom_point(aes(size=-log10(FDR),color=NES)) +
  scale_color_gradientn('NES', 
                        colors=my_palette) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), "lines"))
