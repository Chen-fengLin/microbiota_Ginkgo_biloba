source("Diffabund.R")

bac.alldata = read.table(file="../../00data/bacterial/feature_table_20norm.tsv")
bac.sampleMeta = sampleMeta[colnames(bac.alldata),]

###Soil
bac.soil.res = getDiffByPeriod(data = bac.alldata[,bac.sampleMeta$type=='S'],
                            meta = bac.sampleMeta[bac.sampleMeta$type=='S',])
write.table(bac.soil.res$res12,file = "DiffRes/bac.soil.S1S2diff.tsv",sep = "\t",quote = FALSE)
write.table(bac.soil.res$res23,file = "DiffRes/bac.soil.S2S3diff.tsv",sep = "\t",quote = FALSE)
###Root
bac.root.res = getDiffByPeriod(data = bac.alldata[,bac.sampleMeta$type=='R'],
                                   meta = bac.sampleMeta[bac.sampleMeta$type=='R',])
write.table(bac.root.res$res12,file = "DiffRes/bac.root.S1S2diff.tsv",sep = "\t",quote = FALSE)
write.table(bac.root.res$res23,file = "DiffRes/bac.root.S2S3diff.tsv",sep = "\t",quote = FALSE)
###Leaf
bac.leaf.res = getDiffByPeriod(data = bac.alldata[,bac.sampleMeta$type=='L'],
                                   meta = bac.sampleMeta[bac.sampleMeta$type=='L',])
write.table(bac.leaf.res$res12,file = "DiffRes/bac.leaf.S1S2diff.tsv",sep = "\t",quote = FALSE)
write.table(bac.leaf.res$res23,file = "DiffRes/bac.leaf.S2S3diff.tsv",sep = "\t",quote = FALSE)


####plot####
source("Diffabund.R")
bac.classify = read.table("../../00data/bacterial/classified_table_10clean.tsv",sep="\t")
row.names(bac.classify) = bac.classify$OTU_id
bac.soil.res12 = read.table("DiffRes/bac.soil.S1S2diff.tsv")
bac.soil.res23 = read.table("DiffRes/bac.soil.S2S3diff.tsv")
bac.root.res12 = read.table("DiffRes/bac.root.S1S2diff.tsv")
bac.root.res23 = read.table("DiffRes/bac.root.S2S3diff.tsv")
bac.leaf.res12 = read.table("DiffRes/bac.leaf.S1S2diff.tsv")
bac.leaf.res23 = read.table("DiffRes/bac.leaf.S2S3diff.tsv")
##generate suitable color panel
bac.signifOTU = c(getSignifTaxon(bac.soil.res12,bac.classify),
                  getSignifTaxon(bac.soil.res23,bac.classify),
                  getSignifTaxon(bac.root.res12,bac.classify),
                  getSignifTaxon(bac.root.res23,bac.classify),
                  getSignifTaxon(bac.leaf.res12,bac.classify),
                  getSignifTaxon(bac.leaf.res23,bac.classify))
bac.signifOTU.table = as.data.frame(table(bac.signifOTU))
bac.signifOTU.table = bac.signifOTU.table[order(bac.signifOTU.table[,2],decreasing = TRUE),]
colnames(bac.signifOTU.table)=c('taxon','Freq')
bac.signifOTU.table$taxon = as.character(bac.signifOTU.table$taxon);row.names(bac.signifOTU.table)=NULL
loadS = load("../30occuNet/bacterial.colorPanel.Rdata")
bac.before.colorFrame = eval(as.symbol(loadS))
bac.signifOTU.table$show = bac.before.colorFrame[bac.signifOTU.table$taxon,'show']
bac.signifOTU.table$color = bac.before.colorFrame[bac.signifOTU.table$taxon,'color']
bac.signifOTU.table[is.na(bac.signifOTU.table$show),c('show')]=c('Others')
bac.signifOTU.table[is.na(bac.signifOTU.table$color),c('color')]=c('#666666')

usedColor = c("#72d7c7","#b8b2de","#ff6c57","#65add8","#ffad38",
              "#ace039","#fec6e3","#c4edba","#ffeb33")#寻找其他合适的颜色
moreColors = c('tan3','purple2','orchid1','seagreen4','burlywood2')

highlightIndex = which(bac.signifOTU.table$Freq>=4 & bac.signifOTU.table$show=='Others')
bac.signifOTU.table$show[highlightIndex] = bac.signifOTU.table$taxon[highlightIndex]
bac.signifOTU.table$color[highlightIndex] = moreColors[1:length(highlightIndex)]

saveData = bac.signifOTU.table[c(1,3,4)]
save(saveData,file = "bac.assigned_OTU_color.Rdata")

weakIndex = which(bac.signifOTU.table$Freq<4 & bac.signifOTU.table$show!='Others')
bac.signifOTU.table$show[weakIndex] = "Others"
bac.signifOTU.table$color[weakIndex] = "#666666"
row.names(bac.signifOTU.table)=bac.signifOTU.table$taxon

##开始作图
bac.leafS12.plot = Diffplot(bac.leaf.res12,bac.classify,colorFrame=bac.signifOTU.table)
bac.leafS23.plot = Diffplot(bac.leaf.res23,bac.classify,colorFrame=bac.signifOTU.table)
bac.rootS12.plot = Diffplot(bac.root.res12,bac.classify,colorFrame=bac.signifOTU.table)
bac.rootS23.plot = Diffplot(bac.root.res23,bac.classify,colorFrame=bac.signifOTU.table)
bac.soilS12.plot = Diffplot(bac.soil.res12,bac.classify,colorFrame=bac.signifOTU.table)
bac.soilS23.plot = Diffplot(bac.soil.res23,bac.classify,colorFrame=bac.signifOTU.table)

##图片拼合
fmtPlotAxis = function(plot){
  plot+theme(legend.position = "none")
}
library(ggpubr)
ggpubr::ggarrange(fmtPlotAxis(bac.leafS12.plot)+theme(axis.title.x = element_blank()),
                  fmtPlotAxis(bac.rootS12.plot)+theme(axis.title = element_blank()),
                  fmtPlotAxis(bac.soilS12.plot)+theme(axis.title = element_blank()),
                  fmtPlotAxis(bac.leafS23.plot),
                  fmtPlotAxis(bac.rootS23.plot)+theme(axis.title.y = element_blank()),
                  fmtPlotAxis(bac.soilS23.plot)+theme(axis.title.y = element_blank()),
                  nrow=2,ncol=3,align = "hv")

ggsave(filename = "bac.L12-R12-S12+L23-R23-S23.pdf",
       width = 18,height = 12,units = "cm")

##生成图例
legend.table = unique(bac.signifOTU.table[,c('show','color')])
row.names(legend.table)<-NULL
legend.table = legend.table[order(legend.table$show),]
legend.table = legend.table[c(which(legend.table$show != 'Others' & legend.table$show != 'Unclassified'),
                              which(legend.table$show == 'Others') , which(legend.table$show == 'Unclassified')),]
ggplot(data = data.frame(x=1,y=1,taxon = legend.table$show))+
  geom_point(mapping = aes(x=x,y=y,colour = taxon))+theme_bw()+
  genelateBWtheme()+
  scale_color_manual(values = legend.table$color,breaks = legend.table$show,
                     guide=guide_legend(override.aes = list(shape=16,size=.symbolSize)))
ggsave('bac.legend.pdf',width=12,height=12,units = "cm")

