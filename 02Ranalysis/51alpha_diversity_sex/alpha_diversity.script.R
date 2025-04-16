library(ggplot2)
library(dplyr)
library(ape)
library(vegan)
library(picante)
source("../../public/public.R")

sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)

alpha.bind = read.table("../32alpha_diversity/alpha_diversity.tsv",sep="\t")
alpha.bind$period = sampleMeta[alpha.bind$sampleID,'period']
alpha.bind$sexByName = sampleMeta[alpha.bind$sampleID,'sexByName']

alpha.bind.mean = aggregate(alpha.bind$value, 
                            list(OTUtype = alpha.bind$OTUtype,
                                 params = alpha.bind$params,
                                 stage = alpha.bind$period,
                                 sexByName = alpha.bind$sexByName,
                                 typeByName = alpha.bind$typeByName),
                            mean)
alpha.bind.sd = aggregate(alpha.bind$value, 
                          list(OTUtype = alpha.bind$OTUtype,
                               params = alpha.bind$params,
                               stage = alpha.bind$period,
                               sexByName = alpha.bind$sexByName,
                               typeByName = alpha.bind$typeByName),
                          sd)
alpha.bind.n = aggregate(alpha.bind$value, 
                         list(OTUtype = alpha.bind$OTUtype,
                              params = alpha.bind$params,
                              stage = alpha.bind$period,
                              sexByName = alpha.bind$sexByName,
                              typeByName = alpha.bind$typeByName),
                         length)

alpha.drawdata = data.frame(OTUtype = alpha.bind.mean$OTUtype,
                            params = alpha.bind.mean$params,
                            stage = alpha.bind.mean$stage,
                            sexByName = alpha.bind.mean$sexByName,
                            typeByName = alpha.bind.mean$typeByName,
                            mean = alpha.bind.mean$x,
                            sd = alpha.bind.sd$x,
                            n=alpha.bind.n$x)

alpha.drawdata$se = alpha.drawdata$sd / sqrt(alpha.drawdata$n)

alpha.drawdata$OTUtype = factor(alpha.drawdata$OTUtype,levels = c('bacterial','fungi'),ordered = T)
alpha.drawdata$params = factor(alpha.drawdata$params,levels=c('Shannon','Chao1','PD'),ordered = T)
alpha.drawdata$sexByName = factor(alpha.drawdata$sexByName,levels=c('Female','Male'),ordered = T)
alpha.drawdata$typeByName = factor(alpha.drawdata$typeByName,levels = c('Leaf','Root','Rhizosphere soil'),ordered = TRUE)


star = function(pvalue){
  if(pvalue<0.001)return('***');
  if(pvalue<0.01)return('**');
  if(pvalue<0.05)return('*');
  return(NA)
}


signifData = NULL
for (OTUtype in unique(alpha.bind$OTUtype)) {
  for(params in unique(alpha.bind$params)){
    for(typeByName in unique(alpha.bind$typeByName)){
      for(stage in unique(alpha.bind$period)){
        values=alpha.bind[alpha.bind$OTUtype==OTUtype &
                            alpha.bind$params==params &
                            alpha.bind$typeByName==typeByName&
                            alpha.bind$period==stage,]
        value1 = values[values$sexByName=='Male','value']
        value2 = values[values$sexByName=='Female','value']
        t.pvalue = t.test(value1,value2)$p.value
        w.pvalue = wilcox.test(value1,value2)$p.value
        signifData = rbind(signifData,data.frame(OTUtype=OTUtype,params=params,typeByName=typeByName,stage=stage,
                                         yPos = mean(c(value1,value2)),
                                         t.pvalue = t.pvalue,t.star = star(t.pvalue),
                                         w.pvalue = w.pvalue,w.star = star(w.pvalue)))
      }
    }
  }
}


signifData$OTUtype = factor(signifData$OTUtype,levels = c('bacterial','fungi'),ordered = T)
signifData$params = factor(signifData$params,levels=c('Shannon','Chao1','PD'),ordered = T)
signifData$stage = factor(signifData$stage,levels = c('S1','S2','S3'),ordered = T)
signifData$typeByName = factor(signifData$typeByName,levels = c('Leaf','Root','Rhizosphere soil'),ordered = TRUE)
alpha.drawdata$typeSex = paste0(alpha.drawdata$type,alpha.drawdata$sex)
plot = ggplot(data = alpha.drawdata)+
  geom_errorbar(mapping = aes(x=stage,ymax=mean+se,ymin=mean-se,color=typeByName),
                width=0.2,linewidth=.linewidth,alpha=0.6)+
  geom_point(mapping = aes(x=stage,y=mean,color=typeByName,fill=typeSex,shape=sexByName),
             size=.symbolSize)+
  geom_line(mapping = aes(x=stage,y=mean,color=typeByName,group=sexByName,linetype=sexByName),
            linewidth=.linewidth)+
  geom_text(data=signifData,mapping = aes(x=stage,y=yPos,label=t.star,color=typeByName),
            size = .fontsizeTitle+2,size.unit = 'pt',fontface="bold",show.legend = F)+
  facet_wrap(~OTUtype+typeByName+params,scales="free",ncol = 3,labeller = function(x){return(x[3])})+
  scale_color_manual(values = colorList)+
  scale_linetype_manual(values = c(2,1))+
  scale_shape_manual(values = c(23,23),guide=guide_legend(override.aes = list(fill=c("transparent","#000"))))+
  scale_fill_manual(values = c("transparent",leafColor,"transparent",soilColor,"transparent",rootColor),
                    guide="none")+
  theme_bw()+genelateBWtheme()+
  theme(axis.title= element_blank(),legend.position = 'bottom')+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = '#000',size = .fontsizeTitle,face = "bold"))

plot
ggsave(filename = "alpha_diversity_byTypeStare.pdf",width = 18,height = 36,units = 'cm')

plot+theme(plot.margin = unit(c(0,10,0,10),units = "mm"))
ggsave(filename = "FigS.pdf",width = 18,height = 26,units = 'cm')
