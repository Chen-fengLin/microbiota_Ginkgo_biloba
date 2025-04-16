library(ggplot2)
library(dplyr)
library(ape)
library(vegan)
library(picante)
source("../../public/public.R")

getCleanData <- function(rawdata){
  cleandata = rawdata[which(rowSums(rawdata)>0),]
  cleandata
}

sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)



##bacterial
bac.alldata = read.table(file="../../00data/bacterial/feature_table_10clean.tsv")
bac.alldata = as.data.frame(t(rrarefy(t(bac.alldata),min(colSums(bac.alldata)))))

bac.Shannon <- diversity(getCleanData(bac.alldata),index = "shannon",MARGIN = 2)
bac.estimate <- estimateR(t(getCleanData(bac.alldata)))

bac.tree = read.tree("../../00data/bacterial/bac.tree")
bac.PD = pd(t(getCleanData(bac.alldata)),bac.tree,include.root = FALSE)


nSample = ncol(bac.alldata)
bac.alpha = data.frame(OTUtype="bacterial",
                       value = c(bac.Shannon,bac.estimate[2,],bac.PD$PD),
                       sampleID = c(names(bac.Shannon),colnames(bac.estimate),row.names(bac.PD)),
                       params = rep(c('Shannon','Chao1','PD'),each=nSample))

##fungi
fgi.alldata = read.table(file="../../00data/fungi/feature_table_10clean.tsv")
fgi.alldata = as.data.frame(t(rrarefy(t(fgi.alldata),min(colSums(fgi.alldata)))))

fgi.Shannon <- diversity(getCleanData(fgi.alldata),index = "shannon",MARGIN = 2)
fgi.estimate <- estimateR(t(getCleanData(fgi.alldata)))

fgi.tree = read.tree("../../00data/fungi/fgi.tree")
fgi.PD = pd(t(getCleanData(fgi.alldata)),fgi.tree,include.root = FALSE)


nSample = ncol(fgi.alldata)
fgi.alpha = data.frame(OTUtype="fungi",
                       value = c(fgi.Shannon,fgi.estimate[2,],fgi.PD$PD),
                       sampleID = c(names(fgi.Shannon),colnames(fgi.estimate),row.names(fgi.PD)),
                       params = rep(c('Shannon','Chao1','PD'),each=nSample))


alpha.bind = rbind(bac.alpha,fgi.alpha)

alpha.bind$typeByName = sampleMeta[alpha.bind$sampleID,'typeByName']

alpha.bind.mean = aggregate(alpha.bind$value, 
       list(OTUtype = alpha.bind$OTUtype,
            params = alpha.bind$params,
            typeByName = alpha.bind$typeByName),
       mean)
alpha.bind.sd = aggregate(alpha.bind$value, 
                          list(OTUtype = alpha.bind$OTUtype,
                               params = alpha.bind$params,
                               typeByName = alpha.bind$typeByName),
                          sd)
alpha.bind.n = aggregate(alpha.bind$value, 
                         list(OTUtype = alpha.bind$OTUtype,
                              params = alpha.bind$params,
                              typeByName = alpha.bind$typeByName),
                         length)
alpha.drawdata = data.frame(OTUtype = alpha.bind.mean$OTUtype,
                            params = alpha.bind.mean$params,
                            typeByName = alpha.bind.mean$typeByName,
                            mean = alpha.bind.mean$x,
                            sd = alpha.bind.sd$x,
                            n=alpha.bind.n$x)
alpha.drawdata$se = alpha.drawdata$sd / sqrt(alpha.drawdata$n)

alpha.drawdata$OTUtype = factor(alpha.drawdata$OTUtype,levels = c('bacterial','fungi'),ordered = T)
alpha.drawdata$params = factor(alpha.drawdata$params,levels=c('Shannon','Chao1','PD'),ordered = T)

write.table(alpha.bind,"alpha_diversity.tsv",sep="\t",quote = FALSE)
write.table(alpha.drawdata,"alpha_diversity_meanByType.tsv",sep="\t",quote=FALSE)

star = function(pvalue){
  if(pvalue<0.001)return('***');
  if(pvalue<0.01)return('**');
  if(pvalue<0.05)return('*');
  return("ns")
}
calculateSignif = function(alpha1,alpha2){
  result = NULL
  for (OTUtype in unique(alpha1$OTUtype)) {
    for(params in unique(alpha1$params)){
      value1 = alpha1[alpha1$OTUtype==OTUtype&alpha1$params==params,'value']
      value2 = alpha2[alpha2$OTUtype==OTUtype&alpha2$params==params,'value']
      t.pvalue = t.test(value1,value2)$p.value
      w.pvalue = wilcox.test(value1,value2)$p.value
      result = rbind(result,data.frame(OTUtype=OTUtype,params=params,
                                       t.pvalue = t.pvalue,t.star = star(t.pvalue),
                                       w.pvalue = w.pvalue,w.star = star(w.pvalue)))
    }
  }
  result
}

leafRootSignif = calculateSignif(alpha.bind[alpha.bind$typeByName=='Leaf',],alpha.bind[alpha.bind$typeByName=='Root',])
rootSoilSignif = calculateSignif(alpha.bind[alpha.bind$typeByName=='Root',],alpha.bind[alpha.bind$typeByName=='Rhizosphere soil',])
signifData = rbind(data.frame(xValue=1.5,leafRootSignif),data.frame(xValue=2.5,rootSoilSignif))

signifData$OTUtype = factor(signifData$OTUtype,levels = c('bacterial','fungi'),ordered = T)
signifData$params = factor(signifData$params,levels=c('Shannon','Chao1','PD'),ordered = T)

ggplot(data = alpha.drawdata)+
  geom_bar(mapping = aes(x=typeByName,y=mean,fill=typeByName),
           stat = 'identity',show.legend = F,width=0.6)+
  geom_errorbar(mapping = aes(x=typeByName,ymax=mean+se,ymin=mean-se),
                stat = 'identity',width=0.3,linewidth=.linewidth)+
  geom_text(data=signifData,mapping = aes(x=xValue,y=-Inf,label=t.star),
            size = .fontsize,size.unit = 'pt',fontface="bold")+
  facet_wrap(~OTUtype+params,scales="free",labeller = function(r){r[2]})+
  scale_x_discrete(labels=c('L','R','S'))+
  scale_y_continuous(expand = c(0,0,0.2,0))+
  theme_bw()+genelateBWtheme()+
  theme(axis.title= element_blank())+
  scale_fill_manual(values = colorList)+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = '#000',size = .fontsizeTitle,face = "bold"))

ggsave(filename = "alpha_diversity_byType.pdf",width = 8.4,height = 6.8,units = 'cm')
