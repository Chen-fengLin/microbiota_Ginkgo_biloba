library(rlang)
source("../../public/public.R")
sample_meta = read.table("../../00data/sample_meta.tsv",sep="\t",header = T)
sample_meta = fmtMeta(sample_meta)
library(ggplot2)
colorList = c('#F8766D','#00BA38','#619CFF')

filt_distance = function(matrix,from_ids,to_ids){
  from_ids=from_ids[order(from_ids)]
  to_ids=to_ids[order(to_ids)]
  used_matrix = matrix[from_ids,to_ids]
  number = diag(used_matrix)
  number
}

loadS = load("../20Distance_matrix/bacterial_BC_dismatrix.Rdata")
bacterial_distance = as.matrix(eval(sym(loadS)))
bac.sample_meta = sample_meta[rownames(bacterial_distance),]
bac.stagePairData = NULL
for(type in list('L','R','S')){
  for(stagePair in list(c('S1','S2'),c('S2','S3'),c('S1','S3'))){
    numberData = filt_distance(bacterial_distance,
                               rownames(bac.sample_meta)[bac.sample_meta$type==type & bac.sample_meta$period==stagePair[1]],
                               rownames(bac.sample_meta)[bac.sample_meta$type==type & bac.sample_meta$period==stagePair[2]])
    bac.stagePairData = rbind(bac.stagePairData,data.frame(type=type,stagePair = paste(stagePair,collapse = '-'),distance=numberData))
  }
}
bac.stagePairData$type = factor(bac.stagePairData$type,levels = c('L','R','S'),ordered = TRUE)
bac.stagePairData$stagePair = factor(bac.stagePairData$stagePair,levels = c('S1-S2','S2-S3','S1-S3'),ordered = TRUE)

###计算显著差异###
bac.signif = NULL
for(type1 in c('L','R','S')){
  for(stagePair1 in c('S1-S2','S2-S3','S1-S3')){
    for(type2 in c('L','R','S')){
      for(stagePair2 in c('S1-S2','S2-S3','S1-S3')){
        group1Name = paste(type1,stagePair1);group2Name = paste(type2,stagePair2);
        group1 = bac.stagePairData[bac.stagePairData$type==type1&bac.stagePairData$stagePair==stagePair1,'distance']
        group2 = bac.stagePairData[bac.stagePairData$type==type2&bac.stagePairData$stagePair==stagePair2,'distance']
        p_value = wilcox.test(group1,group2)$p.value
        if(p_value<=0.05 & group1Name<group2Name){
          bac.signif = rbind(bac.signif,data.frame(group1=group1Name,group2=group2Name,p_value=p_value))
        }
      }
    }
  }
}

bac.signif.frame=data.frame(type=rep(c(1,2,3),each=3),
                            stagePair=rep(c(1,2,3),3),
                            label = c('a','a','a','b','bc','c','d','d','d'),
                            y=aggregate(bac.stagePairData$distance,list(bac.stagePairData$stagePair,bac.stagePairData$type),max)$x)

bac.plot = ggplot(data=bac.stagePairData)+
  geom_boxplot(aes(x=type,y=distance,color=stagePair),outlier.shape = 21)+
  theme_bw()+genelateBWtheme()+
  scale_x_discrete(labels=c('Leaf','Root','Rhizosphere soil'))+
  scale_y_continuous(limits = c(0,1),name = 'Bray–Curtis dissimilarity',expand = c(0.05,0,0.1,0))+
  scale_color_manual(values = colorList)+
  theme(axis.title.x = element_blank(),
        legend.position = "inside",legend.position.inside = c(1,1),legend.justification = c(1,1))+
  geom_text(data = bac.signif.frame,
            aes(x=type+(stagePair-2)/4,y=y,label = label),vjust = -0.5,
            size.unit = "pt",size = .fontsizeTitlePT)




loadS = load("../20Distance_matrix/fungi_BCdismatrix.Rdata")
fungi_distance = as.matrix(eval(sym(loadS)))
fgi.sample_meta = sample_meta[rownames(fungi_distance),]
fgi.stagePairData = NULL
for(type in list('L','R','S')){
  for(stagePair in list(c('S1','S2'),c('S2','S3'),c('S1','S3'))){
    numberData = filt_distance(fungi_distance,
                               rownames(fgi.sample_meta)[fgi.sample_meta$type==type & fgi.sample_meta$period==stagePair[1]],
                               rownames(fgi.sample_meta)[fgi.sample_meta$type==type & fgi.sample_meta$period==stagePair[2]])
    fgi.stagePairData = rbind(fgi.stagePairData,data.frame(type=type,stagePair = paste(stagePair,collapse = '-'),distance=numberData))
  }
}
fgi.stagePairData$type = factor(fgi.stagePairData$type,levels = c('L','R','S'),ordered = TRUE)
fgi.stagePairData$stagePair = factor(fgi.stagePairData$stagePair,levels = c('S1-S2','S2-S3','S1-S3'),ordered = TRUE)

###计算显著差异###
fgi.signif = NULL
for(type1 in c('L','R','S')){
  for(stagePair1 in c('S1-S2','S2-S3','S1-S3')){
    for(type2 in c('L','R','S')){
      for(stagePair2 in c('S1-S2','S2-S3','S1-S3')){
        group1Name = paste(type1,stagePair1);group2Name = paste(type2,stagePair2);
        group1 = fgi.stagePairData[fgi.stagePairData$type==type1&fgi.stagePairData$stagePair==stagePair1,'distance']
        group2 = fgi.stagePairData[fgi.stagePairData$type==type2&fgi.stagePairData$stagePair==stagePair2,'distance']
        p_value = wilcox.test(group1,group2)$p.value
        if(p_value<=0.05 & group1Name<group2Name){
          fgi.signif = rbind(fgi.signif,data.frame(group1=group1Name,group2=group2Name,p_value=p_value))
        }
      }
    }
  }
}

fgi.signif.frame=data.frame(type=rep(c(1,2,3),each=3),
                            stagePair=rep(c(1,2,3),3),
                            label = c('a','b','b','c','c','c','d','d','d'),
                            y=aggregate(fgi.stagePairData$distance,list(fgi.stagePairData$stagePair,fgi.stagePairData$type),max)$x)

fgi.plot=ggplot(data=fgi.stagePairData)+
  geom_boxplot(aes(x=type,y=distance,color=stagePair),outlier.shape = 21)+
  theme_bw()+genelateBWtheme()+
  scale_x_discrete(labels=c('Leaf','Root','Rhizosphere soil'))+
  scale_y_continuous(limits = c(0,1),name = 'Bray–Curtis dissimilarity',expand = c(0.05,0,0.1,0))+
  scale_color_manual(values = colorList)+
  theme(axis.title.x = element_blank(),
        legend.position = "inside",legend.position.inside = c(1,1),legend.justification = c(1,1))+
  geom_text(data = fgi.signif.frame,
            aes(x=type+(stagePair-2)/4,y=y,label = label),vjust = -0.5,
            size.unit = "pt",size = .fontsizeTitlePT)


library(ggpubr)
ggpubr::ggarrange(bac.plot,fgi.plot,nrow=1,ncol = 2,labels = c('(a)','(b)'),
                  font.label = labelFont,align = "hv")

ggsave(filename = "distanceBetWeenStage_SameTreeCompare.pdf",width = 18,height = 8,units = 'cm')
