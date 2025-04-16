# library(pheatmap)
library(ggplot2)
library(dplyr)
library(vegan)
source("../../public/public.R")

NSIMUL=1000
NPARALLEL=64

sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)


##bacterial
bac.alldata = read.table(file="../../00data/bacterial/feature_table_10clean.tsv")
bac.sampleMeta = sampleMeta[colnames(bac.alldata),]


bac.tempRes = nestedtemp(t(bac.alldata))
bac.10data = bac.tempRes$comm
bac.temp = bac.tempRes$statistic
bac.temp.smooth = as.data.frame(bac.tempRes$smooth)
bac.temp.check = oecosimu(t(bac.alldata), nestedtemp,method = "c0",
                          nsimul = NSIMUL,parallel = NPARALLEL)
bac.temp.p = bac.temp.check$oecosimu$pval
print("####bac temp check finished")
print(bac.temp.check)

bac.nodfRes = nestednodf(t(bac.alldata))
bac.nodf = bac.nodfRes$statistic['NODF']
bac.nodf.check = oecosimu(t(bac.alldata), nestednodf,method = "c0",
                          nsimul = NSIMUL,parallel = NPARALLEL)
bac.nodf.p = bac.nodf.check$oecosimu$pval[3]
print("####bac nodf check finished")
print(bac.nodf.check)

bac.10data.wide = data.frame(sampleID = rownames(bac.10data),bac.10data,row.names = NULL)
bac.10data.long = tidyr::gather(bac.10data.wide,key = "OTU_id",value = "value",-1)
bac.10data.long$type=sampleMeta[bac.10data.long$sampleID,'type']
bac.10data.long$sampleID = factor(
  bac.10data.long$sampleID,
  levels = rev(bac.10data.wide$sampleID)[order(sampleMeta[rev(bac.10data.wide$sampleID),'type'])],
  ordered = T
)

bac.10data.long$OTU_id = factor(
  bac.10data.long$OTU_id,
  levels = colnames(bac.10data),
  ordered = T
)

bac.plot = ggplot(data=bac.10data.long[bac.10data.long$value>0,])+
  geom_tile(aes(x=OTU_id,y=sampleID,fill = type),show.legend = FALSE)+
  scale_y_discrete(expand = expansion(0,0),labels=c('Leaf','Root','Soil'),
                   breaks=levels(bac.10data.long$sampleID)[c(30,90,150)])+
  scale_x_discrete(expand = expansion(0,0))+
  scale_fill_manual(values = colorList)+
  theme_bw()+genelateBWtheme()+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(colour = '#000000',size = .fontsizeTitle,face = "bold"))+
  geom_line(data = bac.temp.smooth,aes(x=x*ncol(bac.10data),y=(1-y)*nrow(bac.10data)),
            color="#000000",linewidth=.linewidthBold)+
  annotate(geom="text",x=Inf,y=-Inf,vjust=-2.5,hjust=1,size=realLineWidth(.fontsize),
           label=paste0('Nested Temp = ',round(bac.temp,2),", p = ",round(bac.temp.p,3)))+
  annotate(geom="text",x=Inf,y=-Inf,vjust=-0.5,hjust=1,size=realLineWidth(.fontsize),
           label=paste0('Nested NODF = ',round(bac.nodf,2),", p = ",round(bac.nodf.p,3)))






##fungi
fgi.alldata = read.table(file="../../00data/fungi/feature_table_10clean.tsv")
fgi.sampleMeta = sampleMeta[colnames(fgi.alldata),]


fgi.tempRes = nestedtemp(t(fgi.alldata))
fgi.10data = fgi.tempRes$comm
fgi.temp = fgi.tempRes$statistic
fgi.temp.smooth = as.data.frame(fgi.tempRes$smooth)
fgi.temp.check = oecosimu(t(fgi.alldata), nestedtemp,method = "c0",
                          nsimul = NSIMUL,parallel = NPARALLEL)
fgi.temp.p = fgi.temp.check$oecosimu$pval
print("####fgi temp check finished")
print(fgi.temp.check)

fgi.nodfRes = nestednodf(t(fgi.alldata))
fgi.nodf = fgi.nodfRes$statistic['NODF']
fgi.nodf.check = oecosimu(t(fgi.alldata), nestednodf,method = "c0",
                          nsimul = NSIMUL,parallel = NPARALLEL)
fgi.nodf.p = fgi.nodf.check$oecosimu$pval[3]
print("####fgi nodf check finished")
print(fgi.nodf.check)

fgi.10data.wide = data.frame(sampleID = rownames(fgi.10data),fgi.10data,row.names = NULL)
fgi.10data.long = tidyr::gather(fgi.10data.wide,key = "OTU_id",value = "value",-1)
fgi.10data.long$type=sampleMeta[fgi.10data.long$sampleID,'type']
fgi.10data.long$sampleID = factor(
  fgi.10data.long$sampleID,
  levels = rev(fgi.10data.wide$sampleID)[order(sampleMeta[rev(fgi.10data.wide$sampleID),'type'])],
  ordered = T
)

fgi.10data.long$OTU_id = factor(
  fgi.10data.long$OTU_id,
  levels = colnames(fgi.10data),
  ordered = T
)

fgi.plot=ggplot(data=fgi.10data.long[fgi.10data.long$value>0,])+
  geom_tile(aes(x=OTU_id,y=sampleID,fill = type),show.legend = FALSE)+
  scale_y_discrete(expand = expansion(0,0),labels=c('Leaf','Root','Soil'),
                   breaks=levels(fgi.10data.long$sampleID)[c(30,90,150)])+
  scale_x_discrete(expand = expansion(0,0))+
  scale_fill_manual(values = colorList)+
  theme_bw()+genelateBWtheme()+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(colour = '#000000',size = .fontsizeTitle,face = "bold"))+
  geom_line(data = fgi.temp.smooth,aes(x=x*ncol(fgi.10data),y=(1-y)*nrow(fgi.10data)),
            color="#000000",linewidth=.linewidthBold)+
  annotate(geom="text",x=Inf,y=-Inf,vjust=-2.5,hjust=1,size=realLineWidth(.fontsize),
           label=paste0('Nested Temp = ',round(fgi.temp,2),", p = ",round(fgi.temp.p,3)))+
  annotate(geom="text",x=Inf,y=-Inf,vjust=-0.5,hjust=1,size=realLineWidth(.fontsize),
           label=paste0('Nested NODF = ',round(fgi.nodf,2),", p = ",round(fgi.nodf.p,3)))







library(ggpubr)
resplot = ggpubr::ggarrange(bac.plot,fgi.plot,ncol=1,nrow=2,
                  align = "hv",labels = c('(a)','(b)'),font.label = labelFont)
ggsave(file="bacfgi.filter.pdf",resplot,width = 10,height = 6.5,units = 'cm')
