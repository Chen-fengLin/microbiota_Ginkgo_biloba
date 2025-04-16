source("NMDS.R")
alldata = read.table(file="../../00data/fungi/feature_table_20norm.tsv")
sampleMeta = sampleMeta[colnames(alldata),]

##all sample
alldata.clean = getCleanData(alldata)
alldata.NMDS = calculateNMDS(alldata.clean,distance="bray",dim=2,try=10,trymax=500,maxit=200,
                             model="global",autotransform=TRUE,stress=1,noshare=FALSE)
write.table(alldata.NMDS$points,quote = FALSE, col.names = T,row.names = T,sep = "\t",
            file = paste0("NMDS_fgi/00all_NMDS_stress",round(alldata.NMDS$stress,5),".tsv"))
alldata.NMDSplot = plotNMDS(alldata.NMDS$msd_data,
                            aesList = aes(color=typeByName,shape=periodSexByName),size=.symbolSize,
                            ellipseAES=aes(color=typeByName),ellipse.params=list(linewidth=.linewidth))
alldata.NMDS.plot = alldata.NMDSplot$p+
  genelateBWtheme()+
  scale_color_manual(name="",values=colorList,
                     guide=guide_legend(override.aes = list(alpha = 1,shape=21)))+
  scale_shape_manual(name="",values =shapeList2,
                     guide=guide_legend(override.aes = list(alpha = 1)))

alldata.permanova = calculatePERMANOVA(alldata.clean,alldata.clean~type*period*sex)
alldata.permanova

fgi.alldata.NMDS.plot3 = addtext(alldata.NMDS.plot,fontsize=.fontsize,position = c(-1,1),
                             stress = alldata.NMDS$stress,PERMANOVA = alldata.permanova)

save(fgi.alldata.NMDS.plot3,file = "./NMDS_fgi/00alldata.plot.Rdata")


##Soil sample
soildata.clean = getCleanData(alldata[,which(sampleMeta$type=='S')])
soildata.NMDS = calculateNMDS(soildata.clean,distance="bray",dim=2,try=10,trymax=500,maxit=200,
                             model="global",autotransform=TRUE,stress=1,noshare=FALSE)
write.table(soildata.NMDS$points,quote = FALSE, col.names = T,row.names = T,sep = "\t",
            file = paste0("NMDS_fgi/01soil_NMDS_stress",round(soildata.NMDS$stress,5),".tsv"))
soildata.NMDSplot = plotNMDS(soildata.NMDS$msd_data,
                            aesList = aes(shape=period,fill=sexByName),size=.symbolSize,color=soilColor,show.legend=FALSE,
                            ellipseAES=aes(linetype=period),ellipse.params=list(linewidth=.linewidth,color=soilColor))
soildata.NMDS.plot = soildata.NMDSplot$p+
  genelateBWtheme()+
  scale_shape_manual(name="",values =shapeList)+
  scale_fill_manual(name="",values=c('transparent',soilColor))

# soildata.NMDS.plotWithEllipse = soildata.NMDSplot$ellipse+
#   genelateBWtheme(linewidth=0.5,linewidthBold = 0.75,fontsizeTitle = 10,fontsize = 8)+
#   scale_shape_manual(name="",values =c(21,24,22),
#                      guide=guide_legend(override.aes = list(alpha = 1,color="brown")))+
#   scale_fill_manual(name="",values=c('#ffffff00',alpha('brown',0.8)),
#                     guide=guide_legend(override.aes = list(alpha = 1,shape=21,color="brown")))

soildata.permanova = calculatePERMANOVA(soildata.clean,soildata.clean~period*sex)
soildata.permanova


fgi.soildata.NMDS.plot3 = addtext(soildata.NMDS.plot,fontsize=.fontsize,
                             stress = soildata.NMDS$stress,PERMANOVA = soildata.permanova)



##root sample
rootdata.clean = getCleanData(alldata[,which(sampleMeta$type=='R')])
rootdata.NMDS = calculateNMDS(rootdata.clean,distance="bray",dim=2,try=10,trymax=500,maxit=200,
                             model="global",autotransform=TRUE,stress=1,noshare=FALSE)
write.table(rootdata.NMDS$points,quote = FALSE, col.names = T,row.names = T,sep = "\t",
            file = paste0("NMDS_fgi/02root_NMDS_stress",round(rootdata.NMDS$stress,5),".tsv"))
rootdata.NMDSplot = plotNMDS(rootdata.NMDS$msd_data,
                            aesList = aes(shape=period,fill=sexByName),size=.symbolSize,color=rootColor,show.legend=FALSE,
                            ellipseAES=aes(linetype=period),ellipse.params=list(linewidth=.linewidth,color=rootColor))
rootdata.NMDS.plot = rootdata.NMDSplot$p+
  genelateBWtheme()+
  scale_shape_manual(name="",values =shapeList)+
  scale_fill_manual(name="",values=c('transparent',rootColor))

# rootdata.NMDS.plotWithEllipse = rootdata.NMDSplot$ellipse+
#   genelateBWtheme(linewidth=0.5,linewidthBold = 0.75,fontsizeTitle = 10,fontsize = 8)+
#   scale_shape_manual(name="",values =c(21,24,22),
#                      guide=guide_legend(override.aes = list(alpha = 1,color="#772ba6")))+
#   scale_fill_manual(name="",values=c('#ffffff00',alpha('#772ba6',0.8)),
#                     guide=guide_legend(override.aes = list(alpha = 1,shape=21,color="#772ba6")))

rootdata.permanova = calculatePERMANOVA(rootdata.clean,rootdata.clean~period*sex)
rootdata.permanova


fgi.rootdata.NMDS.plot3 = addtext(rootdata.NMDS.plot,fontsize=.fontsize,position = c(-1,1),
                             stress = rootdata.NMDS$stress,PERMANOVA = rootdata.permanova)


##leaf sample
leafdata.clean = getCleanData(alldata[,which(sampleMeta$type=='L')])
leafdata.NMDS = calculateNMDS(leafdata.clean,distance="bray",dim=2,try=10,trymax=500,maxit=200,
                             model="global",autotransform=TRUE,stress=1,noshare=FALSE)
write.table(leafdata.NMDS$points,quote = FALSE, col.names = T,row.names = T,sep = "\t",
            file = paste0("NMDS_fgi/03leaf_NMDS_stress",round(leafdata.NMDS$stress,5),".tsv"))
leafdata.NMDSplot = plotNMDS(leafdata.NMDS$msd_data,
                            aesList = aes(shape=period,fill=sexByName),size=.symbolSize,color=leafColor,show.legend=FALSE,
                            ellipseAES=aes(linetype=period),ellipse.params=list(linewidth=.linewidth,color=leafColor))
leafdata.NMDS.plot = leafdata.NMDSplot$p+
  genelateBWtheme()+
  scale_shape_manual(name="",values =shapeList)+
  scale_fill_manual(name="",values=c('transparent',leafColor))

leafdata.permanova = calculatePERMANOVA(leafdata.clean,leafdata.clean~period*sex)
leafdata.permanova

fgi.leafdata.NMDS.plot3 = addtext(leafdata.NMDS.plot,fontsize=.fontsize,position = c(-1,1),
                             stress = leafdata.NMDS$stress,PERMANOVA = leafdata.permanova)



# library(ggpubr)
# arrange_plot = ggpubr::ggarrange(leafdata.NMDS.plot3,
#                   rootdata.NMDS.plot3+theme(axis.title.y = element_blank()),
#                   soildata.NMDS.plot3+theme(axis.title.y = element_blank()),
#                   nrow = 1,labels = c('(c)'),font.label=labelFont,align='hv')
# ggsave(file="./NMDS_fgi/10_type_arranged_plot.pdf",width=14,height=4.4,units="cm")
# save(arrange_plot,file = "./NMDS_fgi/10_type_arranged_plot.Rdata")

fgiplot = list(leaf=fgi.leafdata.NMDS.plot3,
               root=fgi.rootdata.NMDS.plot3,
               soil=fgi.soildata.NMDS.plot3)
save(fgiplot,file = "./NMDS_fgi/smallfgiplot.Rdata")