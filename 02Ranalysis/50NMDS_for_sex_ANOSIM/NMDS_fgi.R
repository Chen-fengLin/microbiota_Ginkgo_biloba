source("NMDS.R")
alldata = read.table(file="../../00data/fungi/feature_table_20norm.tsv")
sampleMeta = sampleMeta[colnames(alldata),]

NMDS.sex <- function(type,period,try=10,trymax=500,maxit=200,position=c(-1,1)){
  data.clean = getCleanData(alldata[,which(sampleMeta$type==type & sampleMeta$period==period)])
  data.NMDS = calculateNMDS(type,period,data.clean,distance="bray",dim=2,try=try,trymax=trymax,maxit=maxit,
                            model="global",autotransform=TRUE,stress=1,noshare=FALSE)
  write.table(data.NMDS$points,quote = FALSE, col.names = T,row.names = T,sep = "\t",
              file = paste0("NMDS_fgi/",type,"_",period,"_NMDSForSex_stress",round(data.NMDS$stress,5),".tsv"))
  data.NMDSplot = plotNMDS(type,period,data.NMDS$msd_data,ellipse.params=list(linewidth=.linewidth,color=color.dir[[type]]))
  data.NMDS.plot = data.NMDSplot$p+genelateBWtheme()
  data.anosim = calculateanosim(data.clean,'sex')

  data.NMDS.plot3 = addtext(data.NMDS.plot,position = position,stress = data.NMDS$stress,ANOSIM = data.anosim)
  data.NMDS.plot3
}

Leaf.S1.plot = NMDS.sex('L','S1',position=c(-1,1))
Leaf.S2.plot = NMDS.sex('L','S2',position=c(-1,1))
Leaf.S3.plot = NMDS.sex('L','S3',position=c(-1,1))
Root.S1.plot = NMDS.sex('R','S1',position=c(-1,1))
Root.S2.plot = NMDS.sex('R','S2',position=c(-1,1))
Root.S3.plot = NMDS.sex('R','S3',position=c(-1,1))
Soil.S1.plot = NMDS.sex('S','S1',position=c(-1,1))
Soil.S2.plot = NMDS.sex('S','S2',position=c(-1,1))
Soil.S3.plot = NMDS.sex('S','S3',position=c(-1,1))



fgi.plotList = list(Leaf.S1.plot+theme(axis.title.x = element_blank()),
                  Leaf.S2.plot+theme(axis.title = element_blank()),
                  Leaf.S3.plot+theme(axis.title = element_blank()),
                  Root.S1.plot+theme(axis.title.x = element_blank()),
                  Root.S2.plot+theme(axis.title = element_blank()),
                  Root.S3.plot+theme(axis.title = element_blank()),
                  Soil.S1.plot,
                  Soil.S2.plot+theme(axis.title.y = element_blank())+scale_y_continuous(expand = c(0.1,0,0.1,0)),
                  Soil.S3.plot+theme(axis.title.y = element_blank()))

save(fgi.plotList,file="fgi.plotList.R")


library(rlang)
loadS = load("bac.plotList.Rdata")
bac.plot.List = eval(sym(loadS))
loadS = load("fgi.plotList.Rdata")
fgi.plot.List = eval(sym(loadS))

library(ggpubr)
ggpubr::ggarrange(plotlist = c(bac.plot.List,fgi.plot.List),
                  nrow = 6,ncol = 3,align = "hv")
ggsave(filename = "bac.fgi.sex.pdf",width = 19,height = 36,units = "cm")
