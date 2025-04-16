source("../../public/public.R")
sampleMeta = read.table(file = "../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)
sampleMeta = sampleMeta[order(sampleMeta$sampleID),]
library(ggplot2)

bac.icamp = read.csv("../03iCAMP_SES.RC_add/bac.iCAMP.process.bNRIiRCa.csv")[-1]

bac.btnNiche.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(period in levels(sampleMeta$period)){
    sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period,]$sampleID
    bac.btnNiche.icamp = rbind(bac.btnNiche.icamp,bac.icamp[bac.icamp$sample1 %in% sampleList & bac.icamp$sample2 %in% sampleList,])
  }
}

bac.btnNiche.icamp.bind = 
  data.frame(bac.btnNiche.icamp,
             sampleMeta[bac.btnNiche.icamp$sample1,c("type","period","typeByName")])

bac.plotdata = tidyr::gather(bac.btnNiche.icamp.bind,key = "process",value="contri",Heterogeneous.Selection:Drift.and.Others)
process = c('Heterogeneous.Selection','Homogeneous.Selection','Dispersal.Limitation','Homogenizing.Dispersal','Drift.and.Others')
bac.plotdata$process = factor(bac.plotdata$process,levels = process,ordered = T)

bac.plotdata.mean = aggregate.data.frame(bac.plotdata$contri,by = list(type=bac.plotdata$typeByName,period=bac.plotdata$period,process=bac.plotdata$process),mean)

bar.bac.plot = ggplot(bac.plotdata.mean)+
  geom_bar(aes(x=period,y=x*100,fill = process),stat = "identity")+
  facet_wrap(~type)+theme_bw()+genelateBWtheme()+
  scale_y_continuous(name = "Relative importance",expand = c(0,0),
                     labels = function(p){paste0(p,'%')})+
  scale_fill_manual(name="Ecological process",values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#666"),
                    breaks = process,
                    labels=c('Heterogeneous selection','Homogeneous selection','Dispersal limitation','Homogenizing dispersal','Drift and others'))+
  theme(axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0,'pt'),
        # legend.title = element_text(face = "bold",size = .fontsizeTitlePT),
        axis.title.y = element_text(colour = "#000000",size = .fontsizePT),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",size = .fontsizePT),
        legend.key.size = unit(.symbolSizePT*3,units = 'pt'))
# 
# bac.nich.select = NULL
# for(type in levels(sampleMeta$type)){
#   sampleList = sampleMeta[sampleMeta$type==type,]$sampleID
#   bac.nich.select = rbind(bac.nich.select,bac.icamp[bac.icamp$sample1 %in% sampleList & bac.icamp$sample2 %in% sampleList,])
# }
# bac.nich.select.bind = 
#   data.frame(bac.nich.select,
#              sampleMeta[bac.nich.select$sample1,c("type","typeByName")])
# 
# bac.nich.select.bind.mean = tapply(bac.nich.select.bind$Homogeneous.Selection,bac.nich.select.bind$typeByName,mean)
# bac.nich.select.bind.sd = tapply(bac.nich.select.bind$Homogeneous.Selection,bac.nich.select.bind$typeByName,sd)
# bac.nich.select.bind.n = tapply(bac.nich.select.bind$Homogeneous.Selection,bac.nich.select.bind$typeByName,length)
# bac.nich.select.bind.se = bac.nich.select.bind.sd / sqrt(bac.nich.select.bind.n)
# bac.nich.select.bind.barplot = data.frame(type=factor(names(bac.nich.select.bind.mean),levels = c('Leaf','Root','Rhizosphere soil'),ordered = T),
#                                           mean=bac.nich.select.bind.mean,
#                                           sd=bac.nich.select.bind.sd,
#                                           n=bac.nich.select.bind.n,
#                                           se=bac.nich.select.bind.se)
#   
# ggplot()+
#   geom_bar(data=bac.nich.select.bind.barplot,aes(x=type,y=mean*100),fill="#1F78B4",stat="identity")+
#   geom_errorbar(data = bac.nich.select.bind.barplot,aes(x=type,ymax=mean*100+se*100,ymin=mean*100-se*100),width=0.5)+
#   theme_bw()+genelateBWtheme()+
#   scale_y_continuous(name = "Relative importance (%)",expand = c(0,0,0.05,0))+
#   theme(axis.title.x = element_blank(),
#         legend.title = element_text(face = "bold",size = .fontsizeTitlePT),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold",size = .fontsizeTitlePT),
#         legend.key.size = unit(.symbolSizePT*4,units = 'pt'))


fgi.icamp = read.csv("../03iCAMP_SES.RC_add/fgi.iCAMP.process.bNRIiRCa.csv")[-1]

fgi.btnNiche.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(period in levels(sampleMeta$period)){
    sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period,]$sampleID
    fgi.btnNiche.icamp = rbind(fgi.btnNiche.icamp,fgi.icamp[fgi.icamp$sample1 %in% sampleList & fgi.icamp$sample2 %in% sampleList,])
  }
}

fgi.btnNiche.icamp.bind = 
  data.frame(fgi.btnNiche.icamp,
             sampleMeta[fgi.btnNiche.icamp$sample1,c("type","period","typeByName")])

fgi.plotdata = tidyr::gather(fgi.btnNiche.icamp.bind,key = "process",value="contri",Heterogeneous.Selection:Drift.and.Others)
process = c('Heterogeneous.Selection','Homogeneous.Selection','Dispersal.Limitation','Homogenizing.Dispersal','Drift.and.Others')
fgi.plotdata$process = factor(fgi.plotdata$process,levels = process,ordered = T)

fgi.plotdata.mean = aggregate.data.frame(fgi.plotdata$contri,by = list(type=fgi.plotdata$typeByName,period=fgi.plotdata$period,process=fgi.plotdata$process),mean)

bar.fgi.plot = ggplot(fgi.plotdata.mean)+
  geom_bar(aes(x=period,y=x*100,fill = process),stat = "identity")+
  facet_wrap(~type)+theme_bw()+genelateBWtheme()+
  scale_y_continuous(name = "Relative importance",expand = c(0,0),
                     labels = function(p){paste0(p,'%')})+
  scale_fill_manual(name="Ecological process",values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#666"),
                    breaks = process,
                    labels=c('Heterogeneous selection','Homogeneous selection','Dispersal limitation','Homogenizing dispersal','Drift and others'))+
  theme(axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0,'pt'),
        # legend.title = element_text(face = "bold",size = .fontsizeTitlePT),
        axis.title.y = element_text(colour = "#000000",size = .fontsizePT),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",size = .fontsizePT),
        legend.key.size = unit(.symbolSizePT*3,units = 'pt'))

library(ggpubr)
ggarrange(bar.bac.plot,bar.fgi.plot,nrow = 2,ncol = 1,align = "hv",labels = c('(a)','(b)'),
          common.legend = T,legend = "bottom",font.label = labelFont)
ggsave("bac.fgi.arrangeProcess_add_SES.RC.pdf",width = 8,height = 7,units = 'cm')
