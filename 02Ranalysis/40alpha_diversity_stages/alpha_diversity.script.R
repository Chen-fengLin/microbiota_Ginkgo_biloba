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

alpha.bind.mean = aggregate(alpha.bind$value, 
                            list(OTUtype = alpha.bind$OTUtype,
                                 params = alpha.bind$params,
                                 stage = alpha.bind$period,
                                 typeByName = alpha.bind$typeByName),
                            mean)
alpha.bind.sd = aggregate(alpha.bind$value, 
                          list(OTUtype = alpha.bind$OTUtype,
                               params = alpha.bind$params,
                               stage = alpha.bind$period,
                               typeByName = alpha.bind$typeByName),
                          sd)
alpha.bind.n = aggregate(alpha.bind$value, 
                         list(OTUtype = alpha.bind$OTUtype,
                              params = alpha.bind$params,
                              stage = alpha.bind$period,
                              typeByName = alpha.bind$typeByName),
                         length)

alpha.drawdata = data.frame(OTUtype = alpha.bind.mean$OTUtype,
                            params = alpha.bind.mean$params,
                            stage = alpha.bind.mean$stage,
                            typeByName = alpha.bind.mean$typeByName,
                            mean = alpha.bind.mean$x,
                            sd = alpha.bind.sd$x,
                            n=alpha.bind.n$x)

alpha.drawdata$se = alpha.drawdata$sd / sqrt(alpha.drawdata$n)

alpha.drawdata$OTUtype = factor(alpha.drawdata$OTUtype,levels = c('bacterial','fungi'),ordered = T)
alpha.drawdata$params = factor(alpha.drawdata$params,levels=c('Shannon','Chao1','PD'),ordered = T)
alpha.drawdata$typeByName = factor(alpha.drawdata$typeByName,levels = c('Leaf','Root','Rhizosphere soil'),ordered = TRUE)


star = function(pvalue){
  if(pvalue<0.001)return('***');
  if(pvalue<0.01)return('**');
  if(pvalue<0.05)return('*');
  return(NA)
}
calculateSignif = function(alpha1,alpha2){
  result = NULL
  alpha1 = alpha1[order(alpha1$sampleID),]
  alpha2 = alpha2[order(alpha2$sampleID),]
  for (OTUtype in unique(alpha1$OTUtype)) {
    for(params in unique(alpha1$params)){
      for(typeByName in unique(alpha1$typeByName)){
        value1 = alpha1[alpha1$OTUtype==OTUtype&alpha1$params==params&alpha1$typeByName==typeByName,'value']
        value2 = alpha2[alpha2$OTUtype==OTUtype&alpha2$params==params&alpha2$typeByName==typeByName,'value']
        t.pvalue = t.test(value1,value2,paired = TRUE)$p.value
        w.pvalue = wilcox.test(value1,value2,paired = TRUE)$p.value
        result = rbind(result,data.frame(OTUtype=OTUtype,params=params,typeByName=typeByName,
                                         yPos = mean(c(value1,value2)),
                                         t.pvalue = t.pvalue,t.star = star(t.pvalue),
                                         w.pvalue = w.pvalue,w.star = star(w.pvalue)))
      }
    }
  }
  result
}

P1P2Signif = calculateSignif(alpha.bind[alpha.bind$period=='S1',],alpha.bind[alpha.bind$period=='S2',])
P2P3Signif = calculateSignif(alpha.bind[alpha.bind$period=='S2',],alpha.bind[alpha.bind$period=='S3',])



signifData = rbind(data.frame(xValue=1.5,P1P2Signif),data.frame(xValue=2.5,P2P3Signif))

signifData$OTUtype = factor(signifData$OTUtype,levels = c('bacterial','fungi'),ordered = T)
signifData$params = factor(signifData$params,levels=c('Shannon','Chao1','PD'),ordered = T)
signifData$typeByName = factor(signifData$typeByName,levels = c('Leaf','Root','Rhizosphere soil'),ordered = TRUE)

ggplot(data = alpha.drawdata)+
  geom_errorbar(mapping = aes(x=stage,ymax=mean+se,ymin=mean-se,color=typeByName),
                width=0.2,linewidth=.linewidth,alpha=0.6)+
  geom_point(mapping = aes(x=stage,y=mean,color=typeByName),
             shape=18,size=.symbolSize)+
  geom_line(mapping = aes(x=stage,y=mean,color=typeByName,group=typeByName),
            linewidth=.linewidth)+
  
  geom_text(data=signifData,mapping = aes(x=xValue,y=yPos,label=t.star,color=typeByName),
            size = .fontsizeTitle,size.unit = 'pt',fontface="bold",show.legend = F)+
  facet_wrap(~OTUtype+params,scales="free")+
  scale_y_continuous(expand = c(0,0,0.2,0),limits = function(raw){return(c(0,raw[2]));})+
  scale_color_manual(values = colorList)+
  
  theme_bw()+genelateBWtheme()+
  theme(axis.title= element_blank(),legend.position = 'bottom')+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = '#000',size = .fontsizeTitle,face = "bold"))

ggsave(filename = "alpha_diversity_byTypeStare.pdf",width = 12,height = 11,units = 'cm')




shannon_plot = ggplot(data = alpha.drawdata[alpha.drawdata$params=='Shannon',])+
  geom_errorbar(mapping = aes(x=stage,ymax=mean+se,ymin=mean-se,color=typeByName),
                width=0.2,linewidth=.linewidth,alpha=0.6)+
  geom_point(mapping = aes(x=stage,y=mean,color=typeByName),
             shape=18,size=.symbolSize)+
  geom_line(mapping = aes(x=stage,y=mean,color=typeByName,group=typeByName),
            linewidth=.linewidth)+
  
  geom_text(data=signifData[signifData$params=='Shannon',],mapping = aes(x=xValue,y=yPos,label=t.star,color=typeByName),
            size = .fontsizeTitle,size.unit = 'pt',fontface="bold",show.legend = F)+
  facet_wrap(~OTUtype,scales="free",nrow=2)+
  scale_y_continuous(expand = c(0,0,0.2,0),name = "Shannon diversity index",
                     limits = function(raw){return(c(0,raw[2]));})+
  scale_color_manual(values = colorList)+
  theme_bw()+genelateBWtheme()+
  theme(axis.title.x= element_blank(),legend.position = 'right')+
  theme(strip.background = element_blank(),strip.text = element_blank())

chao1_plot = ggplot(data = alpha.drawdata[alpha.drawdata$params=='Chao1',])+
  geom_errorbar(mapping = aes(x=stage,ymax=mean+se,ymin=mean-se,color=typeByName),
                width=0.2,linewidth=.linewidth,alpha=0.6)+
  geom_point(mapping = aes(x=stage,y=mean,color=typeByName),
             shape=18,size=.symbolSize)+
  geom_line(mapping = aes(x=stage,y=mean,color=typeByName,group=typeByName),
            linewidth=.linewidth)+
  
  geom_text(data=signifData[signifData$params=='Chao1',],mapping = aes(x=xValue,y=yPos,label=t.star,color=typeByName),
            size = .fontsizeTitle,size.unit = 'pt',fontface="bold",show.legend = F)+
  facet_wrap(~OTUtype,scales="free",nrow=2)+
  scale_y_continuous(expand = c(0,0,0.2,0),name = "Chao1 richness",
                     limits = function(raw){return(c(0,raw[2]));})+
  scale_color_manual(values = colorList)+
  theme_bw()+genelateBWtheme()+
  theme(axis.title.x= element_blank(),legend.position = 'right')+
  theme(strip.background = element_blank(),strip.text = element_blank())

PD_plot = ggplot(data = alpha.drawdata[alpha.drawdata$params=='PD',])+
  geom_errorbar(mapping = aes(x=stage,ymax=mean+se,ymin=mean-se,color=typeByName),
                width=0.2,linewidth=.linewidth,alpha=0.6)+
  geom_point(mapping = aes(x=stage,y=mean,color=typeByName),
             shape=18,size=.symbolSize)+
  geom_line(mapping = aes(x=stage,y=mean,color=typeByName,group=typeByName),
            linewidth=.linewidth)+
  
  geom_text(data=signifData[signifData$params=='PD',],mapping = aes(x=xValue,y=yPos,label=t.star,color=typeByName),
            size = .fontsizeTitle,size.unit = 'pt',fontface="bold",show.legend = F)+
  facet_wrap(~OTUtype,scales="free",nrow=2)+
  scale_y_continuous(expand = c(0,0,0.2,0),name = "Faith's phylogenetic diversity",
                     limits = function(raw){return(c(0,raw[2]));})+
  scale_color_manual(values = colorList)+
  theme_bw()+genelateBWtheme()+
  theme(axis.title.x= element_blank(),legend.position = 'right')+
  theme(strip.background = element_blank(),strip.text = element_blank())


ggsave(shannon_plot,filename = "shannon_byTypeStar.pdf",width = 7.5,height =8,units = 'cm')
