source("Diffabund.R")

fgi.alldata = read.table(file="../../00data/fungi/feature_table_20norm.tsv")
fgi.sampleMeta = sampleMeta[colnames(fgi.alldata),]

###Soil
fgi.soil.res = getDiffByPeriod(data = fgi.alldata[,fgi.sampleMeta$type=='S'],
                            meta = fgi.sampleMeta[fgi.sampleMeta$type=='S',])
write.table(fgi.soil.res$res12,file = "DiffRes/fgi.soil.S1S2diff.tsv",sep = "\t",quote = FALSE)
write.table(fgi.soil.res$res23,file = "DiffRes/fgi.soil.S2S3diff.tsv",sep = "\t",quote = FALSE)
###Root
fgi.root.res = getDiffByPeriod(data = fgi.alldata[,fgi.sampleMeta$type=='R'],
                                   meta = fgi.sampleMeta[fgi.sampleMeta$type=='R',])
write.table(fgi.root.res$res12,file = "DiffRes/fgi.root.S1S2diff.tsv",sep = "\t",quote = FALSE)
write.table(fgi.root.res$res23,file = "DiffRes/fgi.root.S2S3diff.tsv",sep = "\t",quote = FALSE)
###Leaf
fgi.leaf.res = getDiffByPeriod(data = fgi.alldata[,fgi.sampleMeta$type=='L'],
                                   meta = fgi.sampleMeta[fgi.sampleMeta$type=='L',])
write.table(fgi.leaf.res$res12,file = "DiffRes/fgi.leaf.S1S2diff.tsv",sep = "\t",quote = FALSE)
write.table(fgi.leaf.res$res23,file = "DiffRes/fgi.leaf.S2S3diff.tsv",sep = "\t",quote = FALSE)


####plot####
source("Diffabund.R")
fgi.classify = read.table("../../00data/fungi/classified_table_10clean.tsv",sep="\t")
row.names(fgi.classify) = fgi.classify$OTU_id
fgi.soil.res12 = read.table("DiffRes/fgi.soil.S1S2diff.tsv")
fgi.soil.res23 = read.table("DiffRes/fgi.soil.S2S3diff.tsv")
fgi.root.res12 = read.table("DiffRes/fgi.root.S1S2diff.tsv")
fgi.root.res23 = read.table("DiffRes/fgi.root.S2S3diff.tsv")
fgi.leaf.res12 = read.table("DiffRes/fgi.leaf.S1S2diff.tsv")
fgi.leaf.res23 = read.table("DiffRes/fgi.leaf.S2S3diff.tsv")
##generate suitable color panel
fgi.signifOTU = c(getSignifTaxon(fgi.soil.res12,fgi.classify),
                  getSignifTaxon(fgi.soil.res23,fgi.classify),
                  getSignifTaxon(fgi.root.res12,fgi.classify),
                  getSignifTaxon(fgi.root.res23,fgi.classify),
                  getSignifTaxon(fgi.leaf.res12,fgi.classify),
                  getSignifTaxon(fgi.leaf.res23,fgi.classify))
fgi.signifOTU.table = as.data.frame(table(fgi.signifOTU))
fgi.signifOTU.table = fgi.signifOTU.table[order(fgi.signifOTU.table[,2],decreasing = TRUE),]
colnames(fgi.signifOTU.table)=c('taxon','Freq')
fgi.signifOTU.table$taxon = as.character(fgi.signifOTU.table$taxon);row.names(fgi.signifOTU.table)=NULL
loadS = load("../30occuNet/fungi.colorPanel.Rdata")
fgi.before.colorFrame = eval(as.symbol(loadS))
fgi.signifOTU.table$show = fgi.before.colorFrame[fgi.signifOTU.table$taxon,'show']
fgi.signifOTU.table$color = fgi.before.colorFrame[fgi.signifOTU.table$taxon,'color']
fgi.signifOTU.table[is.na(fgi.signifOTU.table$show),c('show')]=c('Others')
fgi.signifOTU.table[is.na(fgi.signifOTU.table$color),c('color')]=c('#666666')

usedColor = c("#72d7c7","#b8b2de","#ff6c57","#65add8","#ffad38",
              "#ace039","#fec6e3","#c4edba","#ffeb33")#寻找其他合适的颜色
moreColors = c('tan3','purple2','orchid1','seagreen4','burlywood2')

highlightIndex = which(fgi.signifOTU.table$Freq>=4 & fgi.signifOTU.table$show=='Others')
fgi.signifOTU.table$show[highlightIndex] = fgi.signifOTU.table$taxon[highlightIndex]
fgi.signifOTU.table$color[highlightIndex] = moreColors[1:length(highlightIndex)]

saveData = fgi.signifOTU.table[c(1,3,4)]
save(saveData,file = "fgi.assigned_OTU_color.Rdata")

weakIndex = which(fgi.signifOTU.table$Freq<4 & fgi.signifOTU.table$show!='Others')
fgi.signifOTU.table$show[weakIndex] = "Others"
fgi.signifOTU.table$color[weakIndex] = "#666666"
row.names(fgi.signifOTU.table)=fgi.signifOTU.table$taxon

##开始作图
fgi.leafS12.plot = Diffplot(fgi.leaf.res12,fgi.classify,colorFrame=fgi.signifOTU.table)
fgi.leafS23.plot = Diffplot(fgi.leaf.res23,fgi.classify,colorFrame=fgi.signifOTU.table)
fgi.rootS12.plot = Diffplot(fgi.root.res12,fgi.classify,colorFrame=fgi.signifOTU.table)
fgi.rootS23.plot = Diffplot(fgi.root.res23,fgi.classify,colorFrame=fgi.signifOTU.table)
fgi.soilS12.plot = Diffplot(fgi.soil.res12,fgi.classify,colorFrame=fgi.signifOTU.table)
fgi.soilS23.plot = Diffplot(fgi.soil.res23,fgi.classify,colorFrame=fgi.signifOTU.table)

##图片拼合
fmtPlotAxis = function(plot){
  plot+theme(legend.position = "none")
}
library(ggpubr)
ggpubr::ggarrange(fmtPlotAxis(fgi.leafS12.plot)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(fgi.rootS12.plot)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.soilS12.plot)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.leafS23.plot),
                  fmtPlotAxis(fgi.rootS23.plot)+theme(axis.title.y = element_blank()),
                  fmtPlotAxis(fgi.soilS23.plot)+theme(axis.title.y = element_blank()),
                  nrow=2,ncol=3,align = "hv")

ggsave(filename = "fgi.L12-R12-S12+L23-R23-S23.pdf",
       width = 18,height = 12,units = "cm")

##生成图例
legend.table = unique(fgi.signifOTU.table[,c('show','color')])
row.names(legend.table)<-NULL
legend.table = legend.table[order(legend.table$show),]
legend.table = legend.table[c(which(legend.table$show != 'Others' & legend.table$show != 'Unclassified'),
                              which(legend.table$show == 'Others') , which(legend.table$show == 'Unclassified')),]
ggplot(data = data.frame(x=1,y=1,taxon = legend.table$show))+
  geom_point(mapping = aes(x=x,y=y,colour = taxon))+theme_bw()+
  genelateBWtheme()+
  scale_color_manual(values = legend.table$color,breaks = legend.table$show,
                     guide=guide_legend(override.aes = list(shape=16,size=.symbolSize)))
ggsave('fgi.legend.pdf',width=12,height=12,units = "cm")




##合并细菌真菌图
ggpubr::ggarrange(fmtPlotAxis(bac.leafS12.plot)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(bac.leafS23.plot)+theme(axis.title = element_blank()),
                  fmtPlotAxis(bac.rootS12.plot)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(bac.rootS23.plot)+theme(axis.title = element_blank()),
                  fmtPlotAxis(bac.soilS12.plot),
                  fmtPlotAxis(bac.soilS23.plot)+theme(axis.title.y = element_blank()),
                  
                  fmtPlotAxis(fgi.leafS12.plot)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(fgi.leafS23.plot)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.rootS12.plot)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(fgi.rootS23.plot)+theme(axis.title = element_blank()),
                  fmtPlotAxis(fgi.soilS12.plot),
                  fmtPlotAxis(fgi.soilS23.plot)+theme(axis.title.y = element_blank()),
                  nrow=6,ncol=2,align = "hv")

ggsave(filename = "bac+fgi.L12-R12-S12+L23-R23-S23.pdf",
       width = 10,height = 27,units = "cm")
