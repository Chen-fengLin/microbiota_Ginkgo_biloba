source("../../public/public.R")
sampleMeta = read.table(file = "../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)
sampleMeta = sampleMeta[order(sampleMeta$sampleID),]
library(ggplot2)
library(ggsignif)

bac.icamp = read.csv("../02.3iCAMP_confidence_updated/bac.iCAMP.process.CbMPDiCBraya.csv",header = T)[-1]

bac.btnNiche.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(period in levels(sampleMeta$period)){
    sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period,]$sampleID
    bac.btnNiche.icamp = rbind(bac.btnNiche.icamp,bac.icamp[bac.icamp$sample1 %in% sampleList & bac.icamp$sample2 %in% sampleList,])
  }
}

bac.btnNiche.icamp.bind = 
  data.frame(bac.btnNiche.icamp[1:2],contri=(bac.btnNiche.icamp$Homogeneous.Selection)*100,
             sampleMeta[bac.btnNiche.icamp$sample1,c("type","period","typeByName")])

# bac.plotdata = tidyr::gather(bac.btnNiche.icamp.bind,key = "process",value="contri",Heterogeneous.Selection:Drift.and.Others)
# process = c('Heterogeneous.Selection','Homogeneous.Selection','Dispersal.Limitation','Homogenizing.Dispersal','Drift.and.Others')
# bac.plotdata$process = factor(bac.plotdata$process,levels = process,ordered = T)

# bac.plotdata.mean = aggregate.data.frame(bac.plotdata$contri,by = list(type=bac.plotdata$typeByName,period=bac.plotdata$period,process=bac.plotdata$process),mean)
bac.plot=ggplot(data=bac.btnNiche.icamp.bind)+
  geom_violin(aes(x=period,y=contri,fill=typeByName,color=typeByName),
              linewidth=.linewidth)+
  geom_jitter(aes(x=period,y=contri,shape = period,color=typeByName),width=0.1,
             fill="transparent",size = .symbolSize)+
  geom_boxplot(aes(x=period,y=contri),width=0.1,
               outlier.size = .symbolSize,outlier.shape=NA)+
  geom_signif(aes(x=period,y=contri),step_increase=0.1,
              map_signif_level = function(p){sprintf("p = %.2g", p)},
              textsize = realLineWidth(.fontsize),size=.linewidth,
              tip_length = 0.03,
              comparisons = list(c(1,2),c(2,3),c(1,3)))+
  facet_wrap(~typeByName)+
  theme_bw()+genelateBWtheme()+theme_bw()+genelateBWtheme()+
  scale_y_continuous(name = "Relative importance of homogeneous selection",
                     labels = function(p){paste0(p,'%')},
                     expand = c(0,0,0.1,0),
                     breaks = c(0,25,50,75,100))+
  scale_fill_manual(values = alpha(colorList,0.2),aesthetics = c('color','fill'),guide="none")+
  scale_shape_manual(values = shapeList,guide="none")+
  
  theme(axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0,'pt'),
        # legend.title = element_text(face = "bold",size = .fontsizeTitlePT),
        # axis.title.y = element_text(colour = "#000000",size = .fontsizePT),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",size = .fontsizeTitlePT))




fgi.icamp = read.csv("../02.3iCAMP_confidence_updated/fgi.iCAMP.process.CbMPDiCBraya.csv",header = T)[-1]

fgi.btnNiche.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(period in levels(sampleMeta$period)){
    sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period,]$sampleID
    fgi.btnNiche.icamp = rbind(fgi.btnNiche.icamp,fgi.icamp[fgi.icamp$sample1 %in% sampleList & fgi.icamp$sample2 %in% sampleList,])
  }
}

fgi.btnNiche.icamp.bind = 
  data.frame(fgi.btnNiche.icamp[1:2],contri=(fgi.btnNiche.icamp$Homogeneous.Selection)*100,
             sampleMeta[fgi.btnNiche.icamp$sample1,c("type","period","typeByName")])

# fgi.plotdata = tidyr::gather(fgi.btnNiche.icamp.bind,key = "process",value="contri",Heterogeneous.Selection:Drift.and.Others)
# process = c('Heterogeneous.Selection','Homogeneous.Selection','Dispersal.Limitation','Homogenizing.Dispersal','Drift.and.Others')
# fgi.plotdata$process = factor(fgi.plotdata$process,levels = process,ordered = T)

# fgi.plotdata.mean = aggregate.data.frame(fgi.plotdata$contri,by = list(type=fgi.plotdata$typeByName,period=fgi.plotdata$period,process=fgi.plotdata$process),mean)
fgi.plot=ggplot(data=fgi.btnNiche.icamp.bind)+
  geom_violin(aes(x=period,y=contri,fill=typeByName,color=typeByName),
              linewidth=.linewidth)+
  geom_jitter(aes(x=period,y=contri,shape = period,color=typeByName),width=0.1,
              fill="transparent",size = .symbolSize)+
  geom_boxplot(aes(x=period,y=contri),width=0.1,
               outlier.size = .symbolSize,outlier.shape=NA)+
  geom_signif(aes(x=period,y=contri),step_increase=0.1,
              map_signif_level = function(p){sprintf("p = %.2g", p)},
              textsize = realLineWidth(.fontsize),size=.linewidth,
              tip_length = 0.03,
              comparisons = list(c(1,2),c(2,3),c(1,3)))+
  facet_wrap(~typeByName)+
  theme_bw()+genelateBWtheme()+theme_bw()+genelateBWtheme()+
  scale_y_continuous(name = "Relative importance of homogeneous selection",
                     labels = function(p){paste0(p,'%')},
                     expand = c(0,0,0.1,0),
                     breaks = c(0,25,50,75,100))+
  scale_fill_manual(values = alpha(colorList,0.2),aesthetics = c('color','fill'),guide="none")+
  scale_shape_manual(values = shapeList,guide="none")+
  
  theme(axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0,'pt'),
        # legend.title = element_text(face = "bold",size = .fontsizeTitlePT),
        # axis.title.y = element_text(colour = "#000000",size = .fontsizePT),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",size = .fontsizeTitlePT))

library(ggpubr)
ggarrange(bac.plot,fgi.plot,nrow = 2,ncol = 1,align = "hv",labels = c('(a)','(b)'),
          common.legend = T,legend = "bottom",font.label = labelFont)
ggsave("bac.fgi.selection.pdf",width = 18,height = 14,units = 'cm')
