library(ggplot2)
library(dplyr)

source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
rownames(sampleMeta) = sampleMeta$sampleID
# sampleMeta = fmtMeta(sampleMeta) #DEseq do not accept ordered factor



getDiffByPeriod = function(data,meta){
  res.S1S2 = NULL
  res.S2S3 = NULL
  for(OTU_id in rownames(data)){
    S1.data = t(data[OTU_id,meta$period=='S1'])
    S2.data = t(data[OTU_id,meta$period=='S2'])
    S3.data = t(data[OTU_id,meta$period=='S3'])
    
    res.S1S2 = rbind(res.S1S2,
                     data.frame(OTU_id=OTU_id,
                                LFC = log2((mean(S2.data)+.1)/(mean(S1.data)+.1)),
                                p_value = wilcox.test(S1.data,S2.data)$p.value,
                                abund = max(mean(S2.data),mean(S1.data))/1e3,
                                row.names = OTU_id))
    res.S2S3 = rbind(res.S2S3,
                     data.frame(OTU_id=OTU_id,
                                LFC = log2((mean(S3.data)+.1)/(mean(S2.data)+.1)),
                                p_value = wilcox.test(S3.data,S2.data)$p.value,
                                abund = max(mean(S3.data),mean(S2.data))/1e3,
                                row.names = OTU_id))
  }
  res.S1S2 = res.S1S2[res.S1S2$abund>0.1,]
  res.S2S3 = res.S2S3[res.S2S3$abund>0.1,]
  list(res12=res.S1S2,res23=res.S2S3)
}

# plotdiffgene <- function(data,res,pvalue_=0.01,fold_=2,cluster_cols = F,cluster_rows = T){
#   diff<-subset(res,pvalue<pvalue_&padj<pvalue_&abs(log2FoldChange)>fold_)
#   diffUp<-subset(res,pvalue<pvalue_&padj<pvalue_&log2FoldChange>fold_)
#   diffDn<-subset(res,pvalue<pvalue_&padj<pvalue_&log2FoldChange<(0-fold_))
#   res$tag='NdR'
#   res$tag[which(row.names(res) %in% row.names(diffUp))]='UpR'
#   res$tag[which(row.names(res) %in% row.names(diffDn))]='DnR'
#   h.data=data[which(row.names(data) %in% row.names(diff)),]
#   resplot=pheatmap(h.data,scale='row',cluster_rows = cluster_rows,cluster_cols = cluster_cols)
#   
#   if(cluster_rows){
#     order_row = resplot$tree_row$order
#     h.data = data.frame(h.data[order_row,])
#   }
#   if(cluster_cols){
#     order_col = resplot$tree_col$order
#     h.data = data.frame(h.data[,order_col])
#   }
#   m=apply(h.data,1,mean,na.rm=T)
#   s=apply(h.data,1,sd,na.rm=T)
#   h.data=(h.data-m)/s
#   result<-list(v=res,h=h.data,plot=resplot)
#   return(result)
# }

matchTaxon = function(OTUs,classify,levl = 'class'){
  Taxon = classify[OTUs,levl]
  Taxon = ifelse(is.na(Taxon),'Unclassified',Taxon)
  Taxon
}

getSignifTaxon = function(res,classify){
  signif.OTU = row.names(res)[res$p_value<=0.05 & abs(res$LFC)>=1]
  matchTaxon(signif.OTU,classify)
}
Diffplot = function(res,classify,colorFrame){
  drawData = data.frame(OTU_id = res$OTU_id,
                        LFC = res$LFC,
                        p = -log10(res$p_value),
                        Taxon = matchTaxon(res$OTU_id,classify),
                        signif = ifelse(res$p_value<=0.05 & abs(res$LFC)>=1,'T','F'))
  drawData$TaxonShow = colorFrame[drawData$Taxon,'show']
  drawData$TaxonShow = ifelse(is.na(drawData$TaxonShow),'Others',drawData$TaxonShow)
  n.total = nrow(res)
  n.enriched = sum(!(is.na(res$p_value)) & res$p_value<=0.05 & res$LFC>=1)
  n.depleted = sum(!(is.na(res$p_value)) & res$p_value<=0.05 & res$LFC<=-1)
  plot = ggplot(data = drawData)+
    geom_point(mapping = aes(x=p,y=LFC,size=signif,color=TaxonShow))+
    scale_x_continuous(expand = c(0,0,10/100,0),
                       limits = function(o){c(0,max(o[2],8))},
                       name = expression(paste("-",log[10],"(P value)")))+
    scale_y_continuous(name=expression(paste(log[2],"(Fold Change)")),
                       expand = c(0.1,0,0.1,0),
                       limits = function(o){print(o);c(min(o[1],-10),max(o[2],10))})+
    scale_color_manual(values = colorFrame$color,breaks = colorFrame$show)+
    theme_bw()+genelateBWtheme()+
    scale_size_manual(values = c(.symbolSize/3,.symbolSize),breaks = c('F','T'),guide="none")+
    geom_hline(yintercept = c(1,-1),color="#bbb",lty=2)+
    geom_vline(xintercept = -log10(0.05),color="#bbb",lty=2)+
    annotate(geom="text",x=Inf,y=-Inf,hjust=1,vjust=0,size=realLineWidth(.fontsize),
             label=paste0("Increased: ",sprintf("%.1f",n.enriched/n.total*100),"%\n",
                          "Decreased: ",sprintf("%.1f",n.depleted/n.total*100),"%"))
  plot
}
