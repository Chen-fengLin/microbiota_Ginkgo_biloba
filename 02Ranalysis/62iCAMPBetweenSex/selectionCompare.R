source("../../public/public.R")
sampleMeta = read.table(file = "../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)
sampleMeta = sampleMeta[order(sampleMeta$sampleID),]
library(ggplot2)
compareSignif = function(x,data,value,group,method){
  this.value=data[x,value]
  this.meta=factor(data[x,group])
  if(length(levels(this.meta))-2){return(NA)}
  this.value1 = this.value[this.meta==levels(this.meta)[1]]
  this.value2 = this.value[this.meta==levels(this.meta)[2]]
  if(missing(method)){
    return(max(mean(this.value1)+sd(this.value1)/sqrt(length(this.value1)),
               mean(this.value2)+sd(this.value2)/sqrt(length(this.value2))))
  }
  method(this.value1,this.value2)$p.value
}
aggregate.signif <-function(data,by,value,group,method=t.test){
  res = aggregate.data.frame(1:nrow(data),by=by,FUN = compareSignif,
                       data=data,value=value,group=group,method=method)
  res.ypos = aggregate.data.frame(1:nrow(data),by=by,FUN = compareSignif,
                             data=data,value=value,group=group)
  data.frame(
    res[-ncol(res)],
    pvalue=res$x,
    ypos=res.ypos$x,
    star=star(res$x)
  )
  
}


bac.icamp = read.csv("../02iCAMP_confidence/bac.iCAMP.process.CbMPDiCBraya.csv",header = T)[-1]

bac.btnNicheStage.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(period in levels(sampleMeta$period)){
    for(sex in levels(sampleMeta$sex)){
      sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period & sampleMeta$sex==sex,]$sampleID
      bac.btnNicheStage.icamp = rbind(bac.btnNicheStage.icamp,bac.icamp[bac.icamp$sample1 %in% sampleList & bac.icamp$sample2 %in% sampleList,])
    }
  }
}

bac.btnNicheStage.icamp.bind = 
  data.frame(bac.btnNicheStage.icamp[1:2],contri=(bac.btnNicheStage.icamp$Homogeneous.Selection)*100,
             sampleMeta[bac.btnNicheStage.icamp$sample1,c("type","period","typeByName","sex","sexByName")])

##添加一个total stage
bac.btnNicheTotal.icamp=NULL
for(type in levels(sampleMeta$type)){
    for(sex in levels(sampleMeta$sex)){
      sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$sex==sex,]$sampleID
      bac.btnNicheTotal.icamp = rbind(bac.btnNicheTotal.icamp,bac.icamp[bac.icamp$sample1 %in% sampleList & bac.icamp$sample2 %in% sampleList,])
    }
}

bac.btnNicheTotal.icamp.bind = 
  data.frame(bac.btnNicheTotal.icamp[1:2],contri=(bac.btnNicheTotal.icamp$Homogeneous.Selection)*100,
             type=sampleMeta[bac.btnNicheTotal.icamp$sample1,c("type")],
             period="Total",
             sampleMeta[bac.btnNicheTotal.icamp$sample1,c("typeByName","sex","sexByName")])
##

bac.btnNiche.icamp.bind = rbind(bac.btnNicheStage.icamp.bind,bac.btnNicheTotal.icamp.bind)

bac.btnNiche.icamp.bind$typesex = paste(bac.btnNiche.icamp.bind$typeByName,bac.btnNiche.icamp.bind$sexByName)
bac.btnNiche.icamp.bind$typesex = factor(bac.btnNiche.icamp.bind$typesex,
                                         levels = c('Leaf Female','Leaf Male','Root Female','Root Male','Rhizosphere soil Female','Rhizosphere soil Male'),
                                         ordered = T)
bac.btnNiche.icamp.bind.mean = aggregate.data.frame(bac.btnNiche.icamp.bind$contri,by = list(typeByName=bac.btnNiche.icamp.bind$typeByName,period=bac.btnNiche.icamp.bind$period,typesex=bac.btnNiche.icamp.bind$typesex),mean)
bac.btnNiche.icamp.bind.sd = aggregate.data.frame(bac.btnNiche.icamp.bind$contri,by = list(typeByName=bac.btnNiche.icamp.bind$typeByName,period=bac.btnNiche.icamp.bind$period,typesex=bac.btnNiche.icamp.bind$typesex),sd)
bac.btnNiche.icamp.bind.n = aggregate.data.frame(bac.btnNiche.icamp.bind$contri,by = list(typeByName=bac.btnNiche.icamp.bind$typeByName,period=bac.btnNiche.icamp.bind$period,typesex=bac.btnNiche.icamp.bind$typesex),length)
bac.btnNiche.icamp.bind.summary = data.frame(bac.btnNiche.icamp.bind.mean[1:3],
                                             mean=bac.btnNiche.icamp.bind.mean$x,
                                             sd=bac.btnNiche.icamp.bind.sd$x,
                                             n=bac.btnNiche.icamp.bind.n$x,
                                             se=bac.btnNiche.icamp.bind.sd$x/sqrt(bac.btnNiche.icamp.bind.n$x))

bac.btnNiche.icamp.bind.summary$period = 
  factor(bac.btnNiche.icamp.bind.summary$period,
         levels = c('Total','S1','S2','S3'),ordered = T)

bac.signifdata=aggregate.signif(bac.btnNiche.icamp.bind,
                            by = list(typeByName=bac.btnNiche.icamp.bind$typeByName,
                                      period=bac.btnNiche.icamp.bind$period),
                            value='contri',group='sex')

bac.signifdata$period = factor(bac.signifdata$period,
                               levels = c('Total','S1','S2','S3'),ordered = T)




bac.plot=ggplot()+
  geom_bar(data=bac.btnNiche.icamp.bind.summary,
           aes(x=period,y=mean,fill=typesex,color=typesex),
           linewidth=.linewidth,stat = 'identity',
           position = position_dodge2(width = 0.8))+
  geom_errorbar(data=bac.btnNiche.icamp.bind.summary,
                aes(x=period,ymax=mean+se,ymin = mean-se,color=typesex),
                position = position_dodge2(padding= 0.5),show.legend = F)+
  geom_text(data=bac.signifdata,aes(x=period,y=ypos,label=star),vjust=-1)+
  facet_wrap(~typeByName)+
  theme_bw()+genelateBWtheme()+
  scale_y_continuous(name = "Relative importance",
                     labels = function(p){paste0(p,'%')},
                     expand = c(0,0,0.1,0))+
  scale_fill_manual(values = alpha(c('transparent',leafColor,'transparent',rootColor,'transparent',soilColor),0.8))+
  scale_color_manual(values = alpha(rep(colorList,each=2),1))+
  theme(axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0,'pt'),
        # legend.title = element_text(face = "bold",size = .fontsizeTitlePT),
        axis.title.y = element_text(colour = "#000000",size = .fontsizePT),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",size = .fontsizePT))





fgi.icamp = read.csv("../02iCAMP_confidence/fgi.iCAMP.process.CbMPDiCBraya.csv",header = T)[-1]

fgi.btnNicheStage.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(period in levels(sampleMeta$period)){
    for(sex in levels(sampleMeta$sex)){
      sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period & sampleMeta$sex==sex,]$sampleID
      fgi.btnNicheStage.icamp = rbind(fgi.btnNicheStage.icamp,fgi.icamp[fgi.icamp$sample1 %in% sampleList & fgi.icamp$sample2 %in% sampleList,])
    }
  }
}

fgi.btnNicheStage.icamp.bind = 
  data.frame(fgi.btnNicheStage.icamp[1:2],contri=(fgi.btnNicheStage.icamp$Homogeneous.Selection)*100,
             sampleMeta[fgi.btnNicheStage.icamp$sample1,c("type","period","typeByName","sex","sexByName")])

##添加一个total stage
fgi.btnNicheTotal.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(sex in levels(sampleMeta$sex)){
    sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$sex==sex,]$sampleID
    fgi.btnNicheTotal.icamp = rbind(fgi.btnNicheTotal.icamp,fgi.icamp[fgi.icamp$sample1 %in% sampleList & fgi.icamp$sample2 %in% sampleList,])
  }
}

fgi.btnNicheTotal.icamp.bind = 
  data.frame(fgi.btnNicheTotal.icamp[1:2],contri=(fgi.btnNicheTotal.icamp$Homogeneous.Selection)*100,
             type=sampleMeta[fgi.btnNicheTotal.icamp$sample1,c("type")],
             period="Total",
             sampleMeta[fgi.btnNicheTotal.icamp$sample1,c("typeByName","sex","sexByName")])
##

fgi.btnNiche.icamp.bind = rbind(fgi.btnNicheStage.icamp.bind,fgi.btnNicheTotal.icamp.bind)

fgi.btnNiche.icamp.bind$typesex = paste(fgi.btnNiche.icamp.bind$typeByName,fgi.btnNiche.icamp.bind$sexByName)
fgi.btnNiche.icamp.bind$typesex = factor(fgi.btnNiche.icamp.bind$typesex,
                                         levels = c('Leaf Female','Leaf Male','Root Female','Root Male','Rhizosphere soil Female','Rhizosphere soil Male'),
                                         ordered = T)
fgi.btnNiche.icamp.bind.mean = aggregate.data.frame(fgi.btnNiche.icamp.bind$contri,by = list(typeByName=fgi.btnNiche.icamp.bind$typeByName,period=fgi.btnNiche.icamp.bind$period,typesex=fgi.btnNiche.icamp.bind$typesex),mean)
fgi.btnNiche.icamp.bind.sd = aggregate.data.frame(fgi.btnNiche.icamp.bind$contri,by = list(typeByName=fgi.btnNiche.icamp.bind$typeByName,period=fgi.btnNiche.icamp.bind$period,typesex=fgi.btnNiche.icamp.bind$typesex),sd)
fgi.btnNiche.icamp.bind.n = aggregate.data.frame(fgi.btnNiche.icamp.bind$contri,by = list(typeByName=fgi.btnNiche.icamp.bind$typeByName,period=fgi.btnNiche.icamp.bind$period,typesex=fgi.btnNiche.icamp.bind$typesex),length)
fgi.btnNiche.icamp.bind.summary = data.frame(fgi.btnNiche.icamp.bind.mean[1:3],
                                             mean=fgi.btnNiche.icamp.bind.mean$x,
                                             sd=fgi.btnNiche.icamp.bind.sd$x,
                                             n=fgi.btnNiche.icamp.bind.n$x,
                                             se=fgi.btnNiche.icamp.bind.sd$x/sqrt(fgi.btnNiche.icamp.bind.n$x))

fgi.btnNiche.icamp.bind.summary$period = 
  factor(fgi.btnNiche.icamp.bind.summary$period,
         levels = c('Total','S1','S2','S3'),ordered = T)

fgi.signifdata=aggregate.signif(fgi.btnNiche.icamp.bind,
                                by = list(typeByName=fgi.btnNiche.icamp.bind$typeByName,
                                          period=fgi.btnNiche.icamp.bind$period),
                                value='contri',group='sex')

fgi.signifdata$period = factor(fgi.signifdata$period,
                               levels = c('Total','S1','S2','S3'),ordered = T)




fgi.plot=ggplot()+
  geom_bar(data=fgi.btnNiche.icamp.bind.summary,
           aes(x=period,y=mean,fill=typesex,color=typesex),
           linewidth=.linewidth,stat = 'identity',
           position = position_dodge2(width = 0.8))+
  geom_errorbar(data=fgi.btnNiche.icamp.bind.summary,
                aes(x=period,ymax=mean+se,ymin = mean-se,color=typesex),
                position = position_dodge2(padding= 0.5),show.legend = F)+
  geom_text(data=fgi.signifdata,aes(x=period,y=ypos,label=star),vjust=-1)+
  facet_wrap(~typeByName)+
  theme_bw()+genelateBWtheme()+
  scale_y_continuous(name = "Relative importance",
                     labels = function(p){paste0(p,'%')},
                     expand = c(0,0,0.1,0))+
  scale_fill_manual(values = alpha(c('transparent',leafColor,'transparent',rootColor,'transparent',soilColor),0.8))+
  scale_color_manual(values = alpha(rep(colorList,each=2),1))+
  theme(axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0,'pt'),
        # legend.title = element_text(face = "bold",size = .fontsizeTitlePT),
        axis.title.y = element_text(colour = "#000000",size = .fontsizePT),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",size = .fontsizePT))








library(ggpubr)
ggarrange(bac.plot+theme(legend.key.size = unit(.symbolSizePT*3,'pt')),
          fgi.plot+theme(legend.key.size = unit(.symbolSizePT*3,'pt')),
          nrow = 2,ncol = 1,align = "hv",labels = c('(a)','(b)'),
          common.legend = T,legend = "bottom",font.label = labelFont)
ggsave("bac.fgi.selection.pdf",width = 9,height = 8,units = 'cm')
