source("DiffEnrich.R")

fgi.alldata = read.table(file="../../00data/fungi/feature_table_20norm.tsv")
fgi.sampleMeta = sampleMeta[colnames(fgi.alldata),]

resData = list()

for (type in c('L','R','S')) {
  for(stage in c('S1','S2','S3')){
    data.res = SexDiff(data = fgi.alldata[,fgi.sampleMeta$type==type & fgi.sampleMeta$period==stage],
                                meta = fgi.sampleMeta[fgi.sampleMeta$type==type & fgi.sampleMeta$period==stage,])
    n.signif = sum(data.res$p.value<=0.05 & apply(data.res[,c("F.Mean","M.Mean")],1,max)>=0.01)
    print(paste(type,stage,n.signif))
    write.table(data.res,file=paste0("Diff_res_fgi/",type,".",stage,".sexDiff.tsv"),sep="\t",quote = FALSE)
    resData[[type]][[stage]]=data.res
  }
}

save(resData,file = "fgi.diffSex.sex.Rdata")
####plot####
source("DiffEnrich.R")
library(rlang)
loadS = load("fgi.diffSex.sex.Rdata")
diff.sex = eval(sym(loadS))

fgi.classify = read.table("../../00data/fungi/classified_table_10clean.tsv",sep="\t")
row.names(fgi.classify) = fgi.classify$OTU_id

# 
# ##generate suitable color panel
# signifOTU = c(getSignifTaxon(diff.sex[['L']][['S1']],bac.classify),
#               getSignifTaxon(diff.sex[['L']][['S2']],bac.classify),
#               getSignifTaxon(diff.sex[['L']][['S3']],bac.classify),
#               getSignifTaxon(diff.sex[['R']][['S1']],bac.classify),
#               getSignifTaxon(diff.sex[['R']][['S2']],bac.classify),
#               getSignifTaxon(diff.sex[['R']][['S3']],bac.classify),
#               getSignifTaxon(diff.sex[['S']][['S1']],bac.classify),
#               getSignifTaxon(diff.sex[['S']][['S2']],bac.classify),
#               getSignifTaxon(diff.sex[['S']][['S3']],bac.classify))
# signifOTU.table = as.data.frame(table(signifOTU))
# signifOTU.table = signifOTU.table[order(signifOTU.table[,2],decreasing = TRUE),]
# colnames(signifOTU.table)=c('taxon','Freq')
# signifOTU.table$taxon = as.character(signifOTU.table$taxon);row.names(signifOTU.table)=NULL
# loadS = load("../41DESeq2/bac.assigned_OTU_color.Rdata")
# bac.before.colorFrame = eval(sym(loadS))
# row.names(bac.before.colorFrame)=bac.before.colorFrame$taxon
# signifOTU.table$show = bac.before.colorFrame[signifOTU.table$taxon,'show']
# signifOTU.table$color = bac.before.colorFrame[signifOTU.table$taxon,'color']
# # 
# # usedColor = c("#72d7c7","#b8b2de","#ff6c57","#65add8","#ffad38",
# #               "#ace039","#fec6e3","#c4edba","#ffeb33")#寻找其他合适的颜色
# moreColors = c('orchid1','seagreen4','burlywood2')
# # 
# highlightIndex = which(signifOTU.table$Freq>=4 & signifOTU.table$show=='Others')
# signifOTU.table$show[highlightIndex] = signifOTU.table$taxon[highlightIndex]
# signifOTU.table$color[highlightIndex] = moreColors[1:length(highlightIndex)]
# 
# saveData = signifOTU.table[c(1,3,4)]
# save(saveData,file = "bac.assigned_OTU_color.Rdata")
# 
# weakIndex = which(signifOTU.table$Freq<4 & signifOTU.table$show!='Others')
# signifOTU.table$show[weakIndex] = "Others"
# signifOTU.table$color[weakIndex] = "#666666"
# row.names(signifOTU.table)=signifOTU.table$taxon

##开始作图
fgi.resplot = list(
  L.S1 = sexDiffPlot(res=diff.sex,type='L',stage='S1',classify=fgi.classify,colorFrame=NA),
  L.S2 = sexDiffPlot(res=diff.sex,type='L',stage='S2',classify=fgi.classify,colorFrame=NA),
  L.S3 = sexDiffPlot(res=diff.sex,type='L',stage='S3',classify=fgi.classify,colorFrame=NA),
  R.S1 = sexDiffPlot(res=diff.sex,type='R',stage='S1',classify=fgi.classify,colorFrame=NA),
  R.S2 = sexDiffPlot(res=diff.sex,type='R',stage='S2',classify=fgi.classify,colorFrame=NA),
  R.S3 = sexDiffPlot(res=diff.sex,type='R',stage='S3',classify=fgi.classify,colorFrame=NA),
  S.S1 = sexDiffPlot(res=diff.sex,type='S',stage='S1',classify=fgi.classify,colorFrame=NA),
  S.S2 = sexDiffPlot(res=diff.sex,type='S',stage='S2',classify=fgi.classify,colorFrame=NA),
  S.S3 = sexDiffPlot(res=diff.sex,type='S',stage='S3',classify=fgi.classify,colorFrame=NA)
)

# 
# ##图片拼合
fmtPlotAxis = function(plot){
  plot+theme(legend.position = "none")
}
library(ggpubr)
ggpubr::ggarrange(fmtPlotAxis(fgi.resplot$L.S1)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(fgi.resplot$L.S2)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.resplot$L.S3)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.resplot$R.S1)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(fgi.resplot$R.S2)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.resplot$R.S3)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.resplot$S.S1),
                  fmtPlotAxis(fgi.resplot$S.S2)+theme(axis.title.y = element_blank()),
                  fmtPlotAxis(fgi.resplot$S.S3)+theme(axis.title.y = element_blank()),
                  nrow=3,ncol=3,align = "hv")
# 
ggsave(filename = "fgi.L123R123S123.pdf",
       width = 19,height = 19,units = "cm")

# ##生成图例
# legend.table = unique(bac.signifOTU.table[,c('show','color')])
# row.names(legend.table)<-NULL
# legend.table = legend.table[order(legend.table$show),]
# legend.table = legend.table[c(which(legend.table$show != 'Others' & legend.table$show != 'Unclassified'),
#                               which(legend.table$show == 'Others') , which(legend.table$show == 'Unclassified')),]
# ggplot(data = data.frame(x=1,y=1,taxon = legend.table$show))+
#   geom_point(mapping = aes(x=x,y=y,colour = taxon))+theme_bw()+
#   genelateBWtheme()+
#   scale_color_manual(values = legend.table$color,breaks = legend.table$show,
#                      guide=guide_legend(override.aes = list(shape=16,size=.symbolSize)))
# ggsave('bac.legend.pdf',width=12,height=12,units = "cm")




ggpubr::ggarrange(fmtPlotAxis(bac.resplot$L.S1)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(bac.resplot$L.S2)+theme(axis.title = element_blank()),
                  fmtPlotAxis(bac.resplot$L.S3)+theme(axis.title = element_blank()),
                  fmtPlotAxis(bac.resplot$R.S1)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(bac.resplot$R.S2)+theme(axis.title = element_blank()),
                  fmtPlotAxis(bac.resplot$R.S3)+theme(axis.title = element_blank()),
                  fmtPlotAxis(bac.resplot$S.S1),
                  fmtPlotAxis(bac.resplot$S.S2)+theme(axis.title.y = element_blank()),
                  fmtPlotAxis(bac.resplot$S.S3)+theme(axis.title.y = element_blank()),
                  fmtPlotAxis(fgi.resplot$L.S1)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(fgi.resplot$L.S2)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.resplot$L.S3)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.resplot$R.S1)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(fgi.resplot$R.S2)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.resplot$R.S3)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.resplot$S.S1),
                  fmtPlotAxis(fgi.resplot$S.S2)+theme(axis.title.y = element_blank()),
                  fmtPlotAxis(fgi.resplot$S.S3)+theme(axis.title.y = element_blank()),
                  nrow=6,ncol=3,align = "hv")
# 
ggsave(filename = "bac+fgi.L123R123S123.pdf",
       width = 19,height = 38,units = "cm")




##计算统计表
source("DiffEnrich.R")
library(rlang)
loadS = load("fgi.diffSex.sex.Rdata")
diff.sex = eval(sym(loadS))
signifList = list()
for(type in names(diff.sex)){
  signifList[[type]]=list()
  typeData = diff.sex[[type]]
  for(period in names(typeData)){
    periodData.all = typeData[[period]]
    periodData = periodData.all[periodData.all$p <= 0.05,]
    signifList[[type]][[period]]=list(
      Male   = periodData$OTU_id[periodData$Great == 'M'],
      Female = periodData$OTU_id[periodData$Great == 'F'],
      Total  = periodData.all$OTU_id[periodData.all$F.Mean+periodData.all$M.Mean>0]
    )
    
  }
}

resTable = NULL
for(type in names(signifList)){
  typeData = signifList[[type]]
  for(period in names(typeData)){
    periodData = typeData[[period]]
    resTable = rbind(resTable,data.frame(
      type = type,
      stage = period,
      nMale = length(periodData$Male),
      nFemale = length(periodData$Female),
      nTotal = length(unique(c(periodData$Male,periodData$Female))),
      denominator = length(periodData$Total)
    ))
  }
  resTable = rbind(resTable,data.frame(
    type = type,
    stage = 'Total',
    nMale = length(unique(c(typeData$S1$Male,typeData$S2$Male,typeData$S3$Male))),
    nFemale = length(unique(c(typeData$S1$Female,typeData$S2$Female,typeData$S3$Female))),
    nTotal = length(unique(c(typeData$S1$Male,typeData$S2$Male,typeData$S3$Male,typeData$S1$Female,typeData$S2$Female,typeData$S3$Female))),
    denominator = length(unique(c(typeData$S1$Total,typeData$S2$Total,typeData$S3$Total)))
  ))
}

resTable$pMale = resTable$nMale/resTable$denominator*100
resTable$pFemale = resTable$nFemale/resTable$denominator*100
resTable$pTotal = resTable$nTotal/resTable$denominator*100


resTable$sMale = sprintf("%d (%.1f%%)",resTable$nMale,resTable$pMale)
resTable$sFemale = sprintf("%d (%.1f%%)",resTable$nFemale,resTable$pFemale)
resTable$sTotal = sprintf("%d (%.1f%%)",resTable$nTotal,resTable$pTotal)

write.table(resTable,file = "fgi.TableS.tsv",sep="\t",quote = F,row.names = F)
