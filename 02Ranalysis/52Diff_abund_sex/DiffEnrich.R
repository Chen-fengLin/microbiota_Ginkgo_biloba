library(ggplot2)
library(dplyr)

source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
rownames(sampleMeta) = sampleMeta$sampleID
# sampleMeta = fmtMeta(sampleMeta) #DEseq do not accept ordered factor

SexDiff = function(data,meta){
  data.clean<-data[which(rowSums(data)>0),]
  OTUList = row.names(data.clean)
  res=NULL
  for (OTU_id in OTUList) {
    Female.data = t(data.clean[OTU_id,which(meta$sex=='F')])
    Male.data = t(data.clean[OTU_id,which(meta$sex=='M')])
    p.value = wilcox.test(Female.data,Male.data,exact=TRUE)$p.value
    res=rbind(res,data.frame(
      row.names = OTU_id,
      OTU_id=OTU_id,
      F.Mean = mean(Female.data)/1e3,
      M.Mean = mean(Male.data)/1e3,
      Great = c('F','M')[which.max(c(mean(Female.data),mean(Male.data)))],
      LFC = log2((mean(Male.data)+.1)/(mean(Female.data)+.1)),
      p.value = p.value
    ))
  }
  res[order(res$p.value),]
}

matchTaxon = function(OTUs,classify,levl){
  Taxon = classify[OTUs,levl]
  Taxon = ifelse(is.na(Taxon),'Unclassified',Taxon)
  Taxon
}

getSignifTaxon = function(res,classify,levl = 'class'){
  signif.OTU = row.names(res)[res$p.value<=0.05 & apply(res[,c("F.Mean","M.Mean")],1,max)>=0.01]
  matchTaxon(signif.OTU,classify,levl=levl)
}

sexDiffPlot = function(res,type,stage,classify,colorFrame){
  res = res[[type]][[stage]]
  drawData.all = data.frame(OTU_id = res$OTU_id,
                        p = -log10(res$p.value),
                        group=res$Great,
                        LFC=res$LFC,
                        Taxon = matchTaxon(row.names(res),classify,levl = 'class'),
                        abund = apply(res[,c("F.Mean","M.Mean")],1,max))
  # drawData = drawData.all[drawData.all$abund>=0.01 & drawData.all$p>=-log10(0.05),]
  drawData = drawData.all[drawData.all$p>=-log10(0.05),]
  # drawData$TaxonShow = colorFrame[drawData$Taxon,'show']
  n.total = sum(drawData.all$abund>0)
  n.male = sum(drawData$group=="M")
  n.female = sum(drawData$group=="F")
  plot = ggplot(data = drawData)+
    geom_point(mapping = aes(x=abund,y=p,shape=group),size=.symbolSize,
               color=colorList[which(c('L','R','S')==type)])+
    scale_x_log10(name = "Relative abundance (%)",limits = function(old){c(old[1],max(old[2],10))},labels=function(b){sprintf('%s',b)})+
    scale_y_continuous(name=expression(paste("-",log[10],"(P value)")),limits = function(old){c(-log10(0.05),max(old[2],4))})+
    scale_shape_manual(values = shapeList2[which(c('S1','S2','S3','S1','S2','S3')==stage)],
                       breaks = c('F','M'),
                       labels = c('Female biased','Male biased'))+
    theme_bw()+genelateBWtheme()+
    annotate(geom="text",x=Inf,y=Inf,hjust=1,vjust=1,size=realLineWidth(.fontsize),
             label=paste0("Male: ",sprintf("%.1f",n.male/n.total*100),"%\n",
                          "Female: ",sprintf("%.1f",n.female/n.total*100),"%"))
  plot
}
