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
      w.pvalue = method(maleData,femaleData)$p.value
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


OTUCompare.Global <- function(feature,OTUList,sampleMeta){
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
  
  plot = ggplot(data=plotdata.use[plotdata.use$type=='Rhizosphere soil',])+
    geom_bar(aes(x=period,y=mean,fill=typesex,color=typesex),
             stat = "identity",position = position_dodge2(),linewidth=.linewidth)+
    facet_wrap(~type)+
    theme_bw()+genelateBWtheme()+
    geom_errorbar(mapping = aes(x=period,ymax=mean+se,ymin=mean-se,color=typesex),position = position_dodge2(padding=0.5),show.legend = F,linewidth=.linewidth)+
    # geom_text(data=signifData,mapping = aes(x=x,y=ypos,label = w.star),
    #           hjust=0.5,vjust=0.5)+
    scale_color_manual(values = rep(c(soilColor),each=2))+
    scale_fill_manual(values = alpha(c('transparent',soilColor),0.8))+
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

addsignif <- function(res,signif){
  res$sig$ypos$signif = signif
  res$sig$ypos$type = factor(res$sig$ypos$type,levels = c("Leaf","Root","Rhizosphere soil"),ordered = TRUE)
  res$sig$ypos$xpos = match(res$sig$ypos$period,c('S1','S2','S3')) + (match(res$sig$ypos$sexByName,c('Female','Male'))-1.5)/2*0.9
  res$plot = res$plot+
    geom_text(data=res$sig$ypos[res$sig$ypos$type=='Rhizosphere soil',],
              aes(x=xpos,y=ypos,label = signif),hjust=0.5,vjust=-0.5,
              size.unit = 'pt',size=.fontsizePT)
  res
}


#bacterial
bac.classify = read.table("../../00data/bacterial/classified_table_10clean.tsv",sep="\t")
bac.feature = read.table("../../00data/bacterial/feature_table_20norm.tsv",sep="\t")
bac.function = read.table("../../00data/bacterial/function_table_10clean.tsv")
bac.sampleMeta = sampleMeta[colnames(bac.feature),]


OTU.Actinobacteria = bac.classify[!is.na(bac.classify$class) & bac.classify$class=='Actinobacteria',"OTU_id"]

Actinobacteria.res = OTUCompare.Global(bac.feature,OTU.Actinobacteria,bac.sampleMeta)
Actinobacteria.global = addsignif(Actinobacteria.res,c('a','b','b','b','b','b',
                                                             'a','a','a','a','a','a',
                                                             'ab','a','b','ab','b','b'))



OTU.Solirubrobacterales = bac.classify[!is.na(bac.classify$order) & bac.classify$order=='Solirubrobacterales',"OTU_id"]
Solirubrobacterales.res = OTUCompare.Global(bac.feature,OTU.Solirubrobacterales,bac.sampleMeta)
Solirubrobacterales.global = addsignif(Solirubrobacterales.res,c('a','bd','bc','c','bd','d',
                                                       'a','abc','bc','ab','c','c',
                                                       'a','b','ab','b','a','a'))


OTU.Acidimicrobiales = bac.classify[!is.na(bac.classify$order) & bac.classify$order=='Acidimicrobiales',"OTU_id"]
Acidimicrobiales.res = OTUCompare.Global(bac.feature,OTU.Acidimicrobiales,bac.sampleMeta)
Acidimicrobiales.global = addsignif(Acidimicrobiales.res,c('a','b','bc','c','b','bc',
                                                           'a','a','a','a','a','a',
                                                           'a','a','ab','a','b','ab'))


OTU.Rhizomicrobium = bac.classify[!is.na(bac.classify$family) & bac.classify$family=='Rhizomicrobium',"OTU_id"]
Rhizomicrobium.res = OTUCompare.Global(bac.feature,OTU.Rhizomicrobium,bac.sampleMeta)
Rhizomicrobium.global = addsignif(Rhizomicrobium.res,c('a','a','a','a','a','a',
                                                       'ab','ab','ab','a','bc','c',
                                                       'a','b','b','b','b','b'))

ggpubr::ggarrange(Actinobacteria.global$plot+theme(strip.background = element_blank(),strip.text = element_blank()),
                  Solirubrobacterales.global$plot+theme(strip.background = element_blank(),strip.text = element_blank()),
                  Acidimicrobiales.global$plot+theme(strip.background = element_blank(),strip.text = element_blank()),
                  Rhizomicrobium.global$plot+theme(strip.background = element_blank(),strip.text = element_blank()),
                  align = "hv",nrow=2,ncol=2,
                  common.legend = T,legend = "top")

ggsave("Actinobacteria+Solirubrobacterales+Acidimicrobiales+Rhizomicrobium.pdf",
       width = 12,height =8 ,units = "cm")
