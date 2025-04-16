source("../../public/public.R")
sampleMeta = read.table(file = "../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)
sampleMeta = sampleMeta[order(sampleMeta$sampleID),]
library(ggplot2)
alpha = read.table("../32alpha_diversity/alpha_diversity.tsv",sep = "\t")
colorLists = c('#F8766D','#00BA38','#619CFF')

bac.icamp = read.csv("../02iCAMP_confidence/bac.iCAMP.process.CbMPDiCBraya.csv",header = T)[-1]
loadS = load("../20Distance_matrix/bacterial_BC_dismatrix.Rdata")
bac.bc = as.matrix(eval(as.symbol(loadS)))
bac.shannon = alpha[alpha$OTUtype=='bacterial' & alpha$params=='Shannon',]
rownames(bac.shannon)=bac.shannon$sampleID

bac.btnNiche.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(period in levels(sampleMeta$period)){
    # for(sex in levels(sampleMeta$sex)){
      # sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period & sampleMeta$sex==sex,]$sampleID
    sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period,]$sampleID
    
          bac.btnNiche.icamp = rbind(bac.btnNiche.icamp,bac.icamp[bac.icamp$sample1 %in% sampleList & bac.icamp$sample2 %in% sampleList,])
    # }
  }
}

bac.btnNiche.icamp.bind = 
  data.frame(bac.btnNiche.icamp[1:2],contri=(bac.btnNiche.icamp$Homogeneous.Selection)*100,
             sampleMeta[bac.btnNiche.icamp$sample1,c("type","typeByName","period")],
             sexCom = paste(sampleMeta[bac.btnNiche.icamp$sample1,"sexByName"],sampleMeta[bac.btnNiche.icamp$sample2,"sexByName"],sep = "-"))

bac.btnNiche.icamp.bind$bc=diag(bac.bc[bac.btnNiche.icamp.bind$sample1,bac.btnNiche.icamp.bind$sample2])
bac.btnNiche.icamp.bind$shannon = abs(bac.shannon[bac.btnNiche.icamp.bind$sample1,"value"] - bac.shannon[bac.btnNiche.icamp.bind$sample2,"value"])

bac.btnNiche.icamp.lm = NULL
for(type in levels(bac.btnNiche.icamp.bind$typeByName)){
  typedata = bac.btnNiche.icamp.bind[bac.btnNiche.icamp.bind$typeByName==type,]
  linear = summary(lm(bc~contri,data = typedata))
  bac.btnNiche.icamp.lm = 
    rbind(bac.btnNiche.icamp.lm,
           data.frame(typeByName = type,
                      k = linear$coefficients["contri","Estimate"],
                      b = linear$coefficients["(Intercept)","Estimate"],
                      P = linear$coefficients["contri","Pr(>|t|)"],
                      R2= linear$r.squared))
}
bac.btnNiche.icamp.lm$typeByName = factor(bac.btnNiche.icamp.lm$typeByName,levels = c("Leaf","Root","Rhizosphere soil"),ordered = T)

line.bac.plot = ggplot(data = bac.btnNiche.icamp.bind)+
  geom_point(aes(x=contri,y=bc,color=type),shape=21,size = .symbolSize/2,alpha=0.3)+
  geom_smooth(aes(x=contri,y=bc,color=type,fill = type),
              method = 'lm',linewidth=.linewidthBold,alpha = 0.2)+
  geom_text(data=bac.btnNiche.icamp.lm,aes(x=Inf,y=Inf,label = sprintf("R2 = %.3f\np = %.2g",R2,P)),
            size = .fontsizePT,size.unit = 'pt',hjust=1,vjust=1.5)+
  facet_wrap(~typeByName,scales = "free")+
  theme_bw()+genelateBWtheme()+
  scale_color_manual(values = alpha(colorList,0.8),guide="none")+
  scale_fill_manual(values = alpha(colorList,0.1),guide="none")+
  scale_x_continuous(name = "Relative contribution of homogeneous selection (%)",expand = c(0.02,0))+
  scale_y_continuous(name="Bray–Curtis dissimilarity",limits = function(o){c(0,o[2])},expand = c(0.02,0))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = .fontsizePT,colour = "#000000",face="bold"),
        axis.title = element_text(size = .fontsizePT,colour = "#000000"))







fgi.icamp = read.csv("../02iCAMP_confidence/fgi.iCAMP.process.CbMPDiCBraya.csv",header = T)[-1]
loadS = load("../20Distance_matrix/fungi_BCdismatrix.Rdata")
fgi.bc = as.matrix(eval(as.symbol(loadS)))
fgi.shannon = alpha[alpha$OTUtype=='fungi' & alpha$params=='Shannon',]
rownames(fgi.shannon)=fgi.shannon$sampleID

fgi.btnNiche.icamp=NULL
for(type in levels(sampleMeta$type)){
  for(period in levels(sampleMeta$period)){
    # for(sex in levels(sampleMeta$sex)){
    # sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period & sampleMeta$sex==sex,]$sampleID
    sampleList = sampleMeta[sampleMeta$type==type & sampleMeta$period==period,]$sampleID
    
    fgi.btnNiche.icamp = rbind(fgi.btnNiche.icamp,fgi.icamp[fgi.icamp$sample1 %in% sampleList & fgi.icamp$sample2 %in% sampleList,])
    # }
  }
}

fgi.btnNiche.icamp.bind = 
  data.frame(fgi.btnNiche.icamp[1:2],contri=(fgi.btnNiche.icamp$Homogeneous.Selection)*100,
             sampleMeta[fgi.btnNiche.icamp$sample1,c("type","typeByName","period")],
             sexCom = paste(sampleMeta[fgi.btnNiche.icamp$sample1,"sexByName"],sampleMeta[fgi.btnNiche.icamp$sample2,"sexByName"],sep = "-"))

fgi.btnNiche.icamp.bind$bc=diag(fgi.bc[fgi.btnNiche.icamp.bind$sample1,fgi.btnNiche.icamp.bind$sample2])
fgi.btnNiche.icamp.bind$shannon = abs(fgi.shannon[fgi.btnNiche.icamp.bind$sample1,"value"] - fgi.shannon[fgi.btnNiche.icamp.bind$sample2,"value"])

fgi.btnNiche.icamp.lm = NULL
for(type in levels(fgi.btnNiche.icamp.bind$typeByName)){
  typedata = fgi.btnNiche.icamp.bind[fgi.btnNiche.icamp.bind$typeByName==type,]
  linear = summary(lm(bc~contri,data = typedata))
  fgi.btnNiche.icamp.lm = 
    rbind(fgi.btnNiche.icamp.lm,
          data.frame(typeByName = type,
                     k = linear$coefficients["contri","Estimate"],
                     b = linear$coefficients["(Intercept)","Estimate"],
                     P = linear$coefficients["contri","Pr(>|t|)"],
                     R2= linear$r.squared))
}
fgi.btnNiche.icamp.lm$typeByName = factor(fgi.btnNiche.icamp.lm$typeByName,levels = c("Leaf","Root","Rhizosphere soil"),ordered = T)

line.fgi.plot = ggplot(data = fgi.btnNiche.icamp.bind)+
  geom_point(aes(x=contri,y=bc,color=type),shape=21,size = .symbolSize/2,alpha=0.3)+
  geom_smooth(aes(x=contri,y=bc,color=type,fill = type),
              method = 'lm',linewidth=.linewidthBold,alpha = 0.2)+
  geom_text(data=fgi.btnNiche.icamp.lm,aes(x=Inf,y=Inf,label = sprintf("R2 = %.3f\np = %.2g",R2,P)),
            size = .fontsizePT,size.unit = 'pt',hjust=1,vjust=1.5)+
  facet_wrap(~typeByName,scales = "free")+
  theme_bw()+genelateBWtheme()+
  scale_color_manual(values = alpha(colorList,0.8),guide="none")+
  scale_fill_manual(values = alpha(colorList,0.1),guide="none")+
  scale_x_continuous(name = "Relative contribution of homogeneous selection (%)",expand = c(0.02,0))+
  scale_y_continuous(name="Bray–Curtis dissimilarity",limits = function(o){c(0,o[2])},expand = c(0.02,0))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = .fontsizePT,colour = "#000000",face="bold"),
        axis.title = element_text(size = .fontsizePT,colour = "#000000"))

ggpubr::ggarrange(line.bac.plot,line.fgi.plot,nrow=2,ncol=1,align = "hv",
                  labels = c('(a)','(b)'),font.label = labelFont)

ggsave(filename = "HoS2BC.line.pdf",width = 9,height = 8,units = "cm")
