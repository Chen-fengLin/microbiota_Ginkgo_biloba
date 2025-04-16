source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",sep="\t",header = T)
sampleMeta = fmtMeta(sampleMeta)
library(ggplot2)
library(ggsignif)

calculateSignif <- function(combindData,method=wilcox.test){
  yposres = NULL
  for (type in levels(combindData$typeByName)) {
    for(period in levels(combindData$period)){
      maleData = combindData[combindData$typeByName==type & combindData$period==period & combindData$sex == 'M',"abund"]
      femaleData = combindData[combindData$typeByName==type & combindData$period==period & combindData$sex == 'F',"abund"]
      t.pvalue = t.test(maleData,femaleData)$p.value
      w.pvalue = wilcox.test(maleData,femaleData)$p.value
      yposres = rbind(yposres,data.frame(
        type=type,
        period=period,
        sexByName = c('Female','Male'),
        ypos=c(mean(femaleData)+sd(femaleData)/sqrt(length(femaleData)),mean(maleData)+sd(maleData)/sqrt(length(maleData)))
      ))
    }
  }
  signifres=list()
  for (type in levels(combindData$typeByName)) {
    typeData = data.frame(
      S1F = combindData[combindData$typeByName==type & combindData$period=='S1'&combindData$sex=='F',"abund"],
      S1M = combindData[combindData$typeByName==type & combindData$period=='S1'&combindData$sex=='M',"abund"],
      S2F = combindData[combindData$typeByName==type & combindData$period=='S2'&combindData$sex=='F',"abund"],
      S2M = combindData[combindData$typeByName==type & combindData$period=='S2'&combindData$sex=='M',"abund"],
      S3F = combindData[combindData$typeByName==type & combindData$period=='S3'&combindData$sex=='F',"abund"],
      S3M = combindData[combindData$typeByName==type & combindData$period=='S3'&combindData$sex=='M',"abund"]
    )
    signifmatrix = matrix(nrow = 6,ncol = 6,dimnames = list(colnames(typeData),colnames(typeData)))
    for (x in colnames(typeData)) {
      for(y in colnames(typeData)){
        signifmatrix[x,y]=star(method(typeData[[x]],typeData[[y]])$p.value)
      }
    }
    signifres[[type]]=signifmatrix
  }
  list(ypos=yposres,signif=signifres)
}


OTUCompare.Global <- function(feature,OTUList,sampleMeta,diff.signif){
  OTUList = OTUList[OTUList %in% diff.signif$OTU_id]
  feature.func = as.data.frame(t(feature[OTUList,]))
  feature.sum = as.data.frame(apply(feature.func,1,sum)/1e3)
  colnames(feature.sum)='abund'
  OTUCompare(feature.sum,sampleMeta)
}
OTUCompare<- function(feature.sum,sampleMeta){
  combindData = data.frame(sampleMeta[rownames(feature.sum),],feature.sum)
  plotdata.mean = aggregate(combindData$abund,by=list(combindData$typeByName,combindData$period,combindData$sexByName),FUN=mean)
  plotdata.sd = aggregate(combindData$abund,by=list(combindData$typeByName,combindData$period,combindData$sexByName),FUN=sd)
  plotdata.n = aggregate(combindData$abund,by=list(combindData$typeByName,combindData$period,combindData$sexByName),FUN=length)
  plotdata = data.frame(type=plotdata.mean[[1]],
                        period = plotdata.mean[[2]],
                        sex=plotdata.mean[[3]],
                        mean = plotdata.mean$x,
                        sd=plotdata.sd$x,
                        n=plotdata.n$x,
                        se = plotdata.sd$x / sqrt(plotdata.n$x))
  
  plotdata.use = plotdata;plotdata.use$typesex = paste(plotdata.use$type,plotdata.use$sex)
  plotdata.use$type = factor(plotdata.use$type,levels = c("Leaf","Root","Rhizosphere soil"),ordered = TRUE)
  plotdata.use$typesex = factor(plotdata.use$typesex,levels = c("Leaf Female","Leaf Male","Root Female","Root Male","Rhizosphere soil Female","Rhizosphere soil Male"),ordered = TRUE)
  signifData = calculateSignif(combindData)
  
  plot = ggplot(data=plotdata.use)+
    geom_bar(aes(x=period,y=mean,fill=typesex,color=typesex),
             stat = "identity",position = position_dodge2(),linewidth=.linewidth)+
    facet_wrap(~type)+
    theme_bw()+genelateBWtheme()+
    geom_errorbar(mapping = aes(x=period,ymax=mean+se,ymin=mean-se,color=typesex),position = position_dodge2(padding=0.5),show.legend = F,linewidth=.linewidth)+
    # geom_text(data=signifData,mapping = aes(x=x,y=ypos,label = w.star),
    #           hjust=0.5,vjust=0.5)+
    scale_color_manual(values = rep(c(leafColor,rootColor,soilColor),each=2))+
    scale_fill_manual(values = alpha(c('transparent',leafColor,'transparent',rootColor,'transparent',soilColor),0.8))+
    scale_x_discrete()+
    scale_y_continuous(expand = c(0,0,0.1,0),name = 'Relative abundance',
                       labels = function(b){paste0(b,'%')})+
    theme(axis.ticks.x = element_blank(),
          axis.ticks.length.x = unit(0,'pt'),
          axis.title.x = element_blank(),
          legend.key.height = unit(.symbolSizePT*4,'pt'),legend.key.width = unit(.symbolSizePT*6,'pt'),
          strip.background = element_blank())
  list(alldata = combindData,cleandata=plotdata,plot=plot,sig=signifData)
}
# OTUCompare<- function(feature.sum,sampleMeta){
#   combindData = data.frame(sampleMeta[rownames(feature.sum),],feature.sum)
#   plotdata.mean = aggregate(combindData$abund,by=list(combindData$typeByName,combindData$period,combindData$sexByName),FUN=mean)
#   plotdata.sd = aggregate(combindData$abund,by=list(combindData$typeByName,combindData$period,combindData$sexByName),FUN=sd)
#   plotdata.n = aggregate(combindData$abund,by=list(combindData$typeByName,combindData$period,combindData$sexByName),FUN=length)
#   plotdata = data.frame(type=plotdata.mean[[1]],
#                         period = plotdata.mean[[2]],
#                         sex=plotdata.mean[[3]],
#                         mean = plotdata.mean$x,
#                         sd=plotdata.sd$x,
#                         n=plotdata.n$x,
#                         se = plotdata.sd$x / sqrt(plotdata.n$x))
#   
#   plotdata.use = plotdata;plotdata.use$typesex = paste(plotdata.use$type,plotdata.use$sex)
#   plotdata.use$type = factor(plotdata.use$type,levels = c("Leaf","Root","Rhizosphere soil"),ordered = TRUE)
#   plotdata.use$typesex = factor(plotdata.use$typesex,levels = c("Leaf Female","Leaf Male","Root Female","Root Male","Rhizosphere soil Female","Rhizosphere soil Male"),ordered = TRUE)
#   signifData = calculateSignif(combindData)
#   
#   combindData.use = combindData
#   combindData.use$type = combindData.use$typeByName
#   combindData.use$typesex  = paste(combindData.use$typeByName,combindData.use$sexByName)
#   plot = ggplot()+
#     geom_bar(data=plotdata.use,aes(x=period,y=mean,fill=typesex,color=typesex),
#              stat = "identity",position = position_dodge2(),linewidth=.linewidth)+
#     facet_wrap(~type)+
#     theme_bw()+genelateBWtheme()+
#     geom_errorbar(data=plotdata.use,mapping = aes(x=period,ymax=mean+se,ymin=mean-se,color=typesex),position = position_dodge2(padding=0.5),show.legend = F,linewidth=.linewidth)+
#     # geom_text(data=signifData,mapping = aes(x=x,y=ypos,label = w.star),
#     #           hjust=0.5,vjust=0.5)+
#     scale_color_manual(values = rep(c(leafColor,rootColor,soilColor),each=2))+
#     scale_fill_manual(values = alpha(c('transparent',leafColor,'transparent',rootColor,'transparent',soilColor),0.4))+
#     scale_x_discrete()+
#     scale_y_continuous(expand = c(0,0,0.1,0),name = 'Relative abundance',
#                        labels = function(b){paste0(b,'%')})+
#     geom_jitter(data=combindData.use,aes(x=period,y = abund,color=typesex),show.legend = F,shape=16,
#                 position = position_jitterdodge(jitter.width = 0.1),size=.symbolSizePT/4,alpha=0.8)+
#     scale_shape_manual(values = c(21,21))+
#     theme(axis.ticks.x = element_blank(),
#           axis.ticks.length.x = unit(0,'pt'),
#           axis.title.x = element_blank(),
#           legend.key.height = unit(.symbolSizePT*4,'pt'),legend.key.width = unit(.symbolSizePT*6,'pt'),
#           strip.background = element_blank())
#   list(alldata = combindData,cleandata=plotdata,plot=plot,sig=signifData)
# }
addsignif <- function(res,signif){
  res$sig$ypos$signif = signif
  res$sig$ypos$type = factor(res$sig$ypos$type,levels = c("Leaf","Root","Rhizosphere soil"),ordered = TRUE)
  res$sig$ypos$xpos = match(res$sig$ypos$period,c('S1','S2','S3')) + (match(res$sig$ypos$sexByName,c('Female','Male'))-1.5)/2
  res$plot = res$plot+
    geom_text(data=res$sig$ypos,
              aes(x=xpos,y=ypos,label = signif),hjust=0.5,vjust=-0.5,
              size.unit = 'pt',size=.fontsizePT)
  res
}

#bacterial
bac.classify = read.table("../../00data/bacterial/classified_table_10clean.tsv",sep="\t")
bac.feature = read.table("../../00data/bacterial/feature_table_20norm.tsv",sep="\t")
bac.function = read.table("../../00data/bacterial/function_table_10clean.tsv")
bac.sampleMeta = sampleMeta[colnames(bac.feature),]
bac.diff = read.table("../53Diff_abund_sex_enrich/bac.diffData.LRS.tsv")
bac.diff.signif = bac.diff[bac.diff$p.value<=0.05,]


#aerobic_chemoheterotrophy" "chemoheterotrophy"
OTU.chemoheterotrophy = rownames(bac.function)[bac.function$chemoheterotrophy>0]

chemoheterotrophy.res = OTUCompare.Global(bac.feature,OTU.chemoheterotrophy,bac.sampleMeta,bac.diff.signif)
chemoheterotrophy.global = addsignif(chemoheterotrophy.res,c('a','d','cd','bcd','b','bc',
                                                             'ab','ab','ab','ab','a','b',
                                                             'a','a','a','a','a','a'))
# OTUCompare.Global(bac.feature,OTU.aerobic_chemoheterotrophy,bac.sampleMeta,bac.diff.signif)
OTU.photoautotrophy = rownames(bac.function)[bac.function$photoautotrophy>0]
photoautotrophy.res = OTUCompare.Global(bac.feature,OTU.photoautotrophy,bac.sampleMeta,bac.diff.signif)
photoautotrophy.global = addsignif(photoautotrophy.res,c('a','a','a','a','a','a',
                                                         'a','a','a','a','b','b',
                                                         'a','ab','ab','b','ab','b'))

OTU.nitrogen_fixation = rownames(bac.function)[bac.function$nitrogen_fixation>0]
nitrogen_fixation.res = OTUCompare.Global(bac.feature,OTU.nitrogen_fixation,bac.sampleMeta,bac.diff.signif)
nitrogen_fixation.global = addsignif(photoautotrophy.res,c('a','c','c','b','c','c',
                                                           'a','a','a','a','a','a',
                                                           'ab','a','c','bc','c','abc'))

OTU.nitrogen_respiration = 
# rownames(bac.function)[bac.function$nitrate_respiration>0]
# rownames(bac.function)[bac.function$nitrite_respiration>0]
  rownames(bac.function)[bac.function$nitrogen_respiration>0]
nitrogen_respiration.res = OTUCompare.Global(bac.feature,OTU.nitrogen_respiration,bac.sampleMeta,bac.diff.signif)
nitrogen_respiration.global = addsignif(nitrogen_respiration.res,c('a','b','b','b','b','b',
                                                                   'a','a','a','a','b','c',
                                                                   'a','a','a','a','a','a'))


ggsave('chemoheterotrophy.pdf',
       chemoheterotrophy.global$plot+
         scale_y_continuous(expand = c(0,0,0.2,0),
                            name = 'Relative abundance',
                            labels = function(b){paste0(b,'%')})+
         theme(legend.key.height = unit(.symbolSizePT*2,'pt'),
               legend.key.width = unit(.symbolSizePT*3,'pt'),
               strip.text = element_blank()),
       width = 110,height = 35,units = 'mm')

ggsave('photoautotrophy.pdf',
       photoautotrophy.global$plot,
       width = 180,height = 70,units = 'mm')

ggsave('nitrogen_fixation.pdf',
       nitrogen_fixation.global$plot,
       width = 180,height = 70,units = 'mm')

ggsave('nitrogen_respiration.pdf',
       nitrogen_respiration.global$plot,
       width = 180,height = 70,units = 'mm')


leaf.signifOTU =unique(bac.diff.signif[bac.diff.signif$type=='L',]$OTU_id)
root.signifOTU =unique(bac.diff.signif[bac.diff.signif$type=='R',]$OTU_id)
soil.signifOTU =unique(bac.diff.signif[bac.diff.signif$type=='S',]$OTU_id)

sum(leaf.signifOTU  %in% root.signifOTU)
sum(soil.signifOTU  %in% root.signifOTU)
sum(soil.signifOTU  %in% leaf.signifOTU)
sum(soil.signifOTU  %in% leaf.signifOTU & soil.signifOTU  %in% root.signifOTU)

length(leaf.signifOTU)
length(root.signifOTU)
length(soil.signifOTU)


OTU.chemoheterotrophy.L = OTU.chemoheterotrophy[OTU.chemoheterotrophy %in% leaf.signifOTU]
OTU.chemoheterotrophy.R = OTU.chemoheterotrophy[OTU.chemoheterotrophy %in% root.signifOTU]
OTU.chemoheterotrophy.S = OTU.chemoheterotrophy[OTU.chemoheterotrophy %in% soil.signifOTU]

sum(OTU.chemoheterotrophy.L  %in% OTU.chemoheterotrophy.R)
sum(OTU.chemoheterotrophy.L  %in% OTU.chemoheterotrophy.S)
sum(OTU.chemoheterotrophy.R  %in% OTU.chemoheterotrophy.S)
sum(OTU.chemoheterotrophy.R  %in% OTU.chemoheterotrophy.S & OTU.chemoheterotrophy.R  %in% OTU.chemoheterotrophy.L)

length(OTU.chemoheterotrophy.L)
length(OTU.chemoheterotrophy.R)
length(OTU.chemoheterotrophy.S)


write.table(bac.classify[bac.classify$OTU_id %in% OTU.chemoheterotrophy.L,],
            file = "LeafSexDiff.chemoheterotrophy.tsv",row.names = F,quote = F,
            sep="\t",na = '-')
write.table(bac.classify[bac.classify$OTU_id %in% OTU.chemoheterotrophy.R,],
            file = "RootSexDiff.chemoheterotrophy.tsv",row.names = F,quote = F,
            sep="\t",na = '-')





#fungi
fgi.classify = read.table("../../00data/fungi/classified_table_10clean.tsv",sep="\t")
fgi.feature = read.table("../../00data/fungi/feature_table_20norm.tsv",sep="\t")
fgi.function = read.table("../../00data/fungi/function_table_10clean.tsv")
fgi.sampleMeta = sampleMeta[colnames(fgi.feature),]
fgi.diff = read.table("../53Diff_abund_sex_enrich/fgi.diffData.LRS.tsv")
fgi.diff.signif = fgi.diff[fgi.diff$p.value<=0.05,]


#aerobic_chemoheterotrophy" "chemoheterotrophy"
OTU.allsaprotroph = rownames(fgi.function)[rowSums(fgi.function[,colnames(fgi.function)[grep('saprotroph',colnames(fgi.function))]])>0]
OTU.unspecified_saprotroph = rownames(fgi.function)[fgi.function$unspecified_saprotroph>0]
OTU.litter_saprotroph = rownames(fgi.function)[fgi.function$litter_saprotroph>0]

#三种分析策略
# 1、将符合该功能的OTU全部选出来--不佳
# 2、将符合该功能的差异OTU全部选出来-挺好
# OTU.aerobic_chemoheterotrophy = rownames(bac.function)[bac.function$aerobic_chemoheterotrophy>0]
allsaprotroph.res = OTUCompare.GroupByTypePeriod(fgi.feature,OTU.allsaprotroph,fgi.sampleMeta,fgi.diff.signif)
allsaprotroph.global = addsignif(allsaprotroph.res,c('ab','ab','a','b','c','c'))

unspecified_saprotroph.res = OTUCompare.GroupByTypePeriod(fgi.feature,OTU.unspecified_saprotroph,fgi.sampleMeta,fgi.diff.signif)
unspecified_saprotroph.global = addsignif(unspecified_saprotroph.res,c('a','bc','c','c','bc','b',
                                                                       'a','a','a','a','b','b',
                                                                       'b','a','b','ab','b','ab'))

litter_saprotroph.res = OTUCompare.GroupByTypePeriod(fgi.feature,OTU.litter_saprotroph,fgi.sampleMeta,fgi.diff.signif)
litter_saprotroph.global = addsignif(litter_saprotroph.res,c('ab','ab','a','b','c','c',
                                                             'a','a','a','a','a','a',
                                                             'a','b','a','a','a','ab'))

OTU.plant_pathogen = rownames(fgi.function)[fgi.function$plant_pathogen>0]
plant_pathogen.res = OTUCompare.GroupByTypePeriod(fgi.feature,OTU.plant_pathogen,fgi.sampleMeta,fgi.diff.signif)
plant_pathogen.global = addsignif(plant_pathogen.res,c('a','ab','a','b','c','c',
                                                       'ac','a','bc','ab','b','b',
                                                       'ab','a','b','ab','ab','ab'))

ggpubr::ggarrange(allsaprotroph.global$plot,litter_saprotroph.global$plot,plant_pathogen.global$plot,
                  nrow=3,ncol = 1,align = "hv",labels = c('(a)','(b)','(c)'),
                  font.label = labelFont,common.legend = T,legend = 'bottom')
ggsave("allsaprotroph+litter_saprotroph+plant_pathogen.pdf",width=18,height=18,units = 'cm')





leaf.signifOTU =unique(fgi.diff.signif[fgi.diff.signif$type=='L',]$OTU_id)
root.signifOTU =unique(fgi.diff.signif[fgi.diff.signif$type=='R',]$OTU_id)
soil.signifOTU =unique(fgi.diff.signif[fgi.diff.signif$type=='S',]$OTU_id)

sum(leaf.signifOTU  %in% root.signifOTU)
sum(soil.signifOTU  %in% root.signifOTU)
sum(soil.signifOTU  %in% leaf.signifOTU)
sum(soil.signifOTU  %in% leaf.signifOTU & soil.signifOTU  %in% root.signifOTU)

length(leaf.signifOTU)
length(root.signifOTU)
length(soil.signifOTU)



OTU.allsaprotroph.L = OTU.allsaprotroph[OTU.allsaprotroph %in% leaf.signifOTU]
OTU.allsaprotroph.R = OTU.allsaprotroph[OTU.allsaprotroph %in% root.signifOTU]
OTU.allsaprotroph.S = OTU.allsaprotroph[OTU.allsaprotroph %in% soil.signifOTU]

OTU.plant_pathogen.L = OTU.plant_pathogen[OTU.plant_pathogen %in% leaf.signifOTU]
OTU.plant_pathogen.R = OTU.plant_pathogen[OTU.plant_pathogen %in% root.signifOTU]
OTU.plant_pathogen.S = OTU.plant_pathogen[OTU.plant_pathogen %in% soil.signifOTU]



write.table(fgi.classify[fgi.classify$OTU_id %in% OTU.allsaprotroph.L,],
            file = "LeafSexDiff.allsaprotroph.tsv",row.names = F,quote = F,
            sep="\t",na = '-')
write.table(fgi.classify[fgi.classify$OTU_id %in% OTU.allsaprotroph.R,],
            file = "RootSexDiff.allsaprotroph.tsv",row.names = F,quote = F,
            sep="\t",na = '-')

write.table(fgi.classify[fgi.classify$OTU_id %in% OTU.plant_pathogen.L,],
            file = "LeafSexDiff.plant_pathogen.tsv",row.names = F,quote = F,
            sep="\t",na = '-')
write.table(fgi.classify[fgi.classify$OTU_id %in% OTU.plant_pathogen.R,],
            file = "RootSexDiff.plant_pathogen.tsv",row.names = F,quote = F,
            sep="\t",na = '-')
