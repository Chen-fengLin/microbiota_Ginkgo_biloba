source("NMDS.R")
alldata = read.table(file="../../00data/bacterial/feature_table_20norm.tsv")
sampleMeta = sampleMeta[colnames(alldata),]

NMDS.sex <- function(type,period,try=10,trymax=500,maxit=200,position=c(-1,1)){
  data.clean = getCleanData(alldata[,which(sampleMeta$type==type & sampleMeta$period==period)])
  data.NMDS = calculateNMDS(type,period,data.clean,distance="bray",dim=2,try=try,trymax=trymax,maxit=maxit,
                            model="global",autotransform=TRUE,stress=1,noshare=FALSE)
  write.table(data.NMDS$points,quote = FALSE, col.names = T,row.names = T,sep = "\t",
              file = paste0("NMDS_bac/",type,"_",period,"_NMDSForSex_stress",round(data.NMDS$stress,5),".tsv"))
  data.NMDSplot = plotNMDS(type,period,data.NMDS$msd_data,ellipse.params=list(linewidth=.linewidth,color=color.dir[[type]]))
  data.NMDS.plot = data.NMDSplot$p+genelateBWtheme()
  data.anosim = calculateanosim(data.clean,'sex')
  # plot(data.anosim)
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

bac.plotList = list(Leaf.S1.plot+theme(axis.title.x = element_blank()),
                  Leaf.S2.plot+theme(axis.title = element_blank())+scale_y_continuous(expand = c(0.1,0,0.4,0)),
                  Leaf.S3.plot+theme(axis.title = element_blank())+scale_y_continuous(expand = c(0.1,0,0.4,0)),
                  Root.S1.plot+theme(axis.title.x = element_blank())+scale_x_continuous(expand = c(0.2,0,0.1,0)),
                  Root.S2.plot+theme(axis.title = element_blank()),
                  Root.S3.plot+theme(axis.title = element_blank())+scale_y_continuous(expand = c(0.1,0,0.3,0)),
                  Soil.S1.plot+scale_x_continuous(expand = c(0.2,0,0.1,0)),
                  Soil.S2.plot+theme(axis.title.y = element_blank()),
                  Soil.S3.plot+theme(axis.title.y = element_blank()))

save(bac.plotList,file="bac.plotList.R")
