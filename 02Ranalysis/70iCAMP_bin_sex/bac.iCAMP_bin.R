PREFIX="bac"

#加载全部数据，提取所需
loadS = load(paste0("../02iCAMP_confidence/",PREFIX,".iCAMP.iCAMP.Confidence.detail.rda"))
iCAMP.detail = eval(as.symbol(loadS))
  
OTU2bin = iCAMP.detail$detail$taxabin$sp.bin
OTU2bin = data.frame(OTU_id=rownames(OTU2bin),OTU2bin,bin_id=paste0('bin',OTU2bin$bin.id.new))
OTU2bin = OTU2bin[c(-2,-3)]

bMPDi.wide = iCAMP.detail$detail$SigbMPDi
colnames(bMPDi.wide) <- sub("CbMPDi.","",colnames(bMPDi.wide))
bMPDi <- tidyr::gather(bMPDi.wide,key="bin_id",value="CbMPD",c(-1,-2))

BCa.wide = iCAMP.detail$detail$SigBCa
colnames(BCa.wide) <- sub("Cbraya.","",colnames(BCa.wide))
BCa <- tidyr::gather(BCa.wide,key="bin_id",value="CBray",c(-1,-2))

bin.weight.wide = iCAMP.detail$detail$bin.weight
bin.weight <- tidyr::gather(bin.weight.wide,key="bin_id",value="bin.weight",c(-1,-2))

identical(bMPDi[c(1,3)],BCa[c(1,3)])
identical(bMPDi[,1],bin.weight[,1])
identical(bMPDi[,2],bin.weight[,2])
identical(bMPDi[,3],bin.weight[,3])

iCAMP.data = data.frame(
  bin.weight,
  CbMPD = bMPDi$CbMPD,
  CBray = BCa$CBray
)
save(OTU2bin,iCAMP.data,file = paste0(PREFIX,".01iCAMPdata.Rdata"))
#关闭

#加载需要数据，进行bin注释，检查系统发育树是否有问题
PREFIX="bac"
source("../../public/public.R")
load(paste0(PREFIX,".01iCAMPdata.Rdata"))

OTU.classify = read.table("../../00data/bacterial/classified_table_10clean.tsv",sep="\t")
rownames(OTU.classify)=OTU.classify$OTU_id
bin.classify = data.frame(OTU2bin,OTU.classify[rownames(OTU2bin),])
write.table(bin.classify,file = paste0(PREFIX,".02bin.classify.raw.tsv"),
            sep = "\t",quote = FALSE)
#输出后手动整理

library(ape)
library(ggtree)
OTU.tree=read.tree("../../00data/bacterial/bac.tree")
OTU.tree.groupByPhylum = groupOTU(OTU.tree,split(OTU.classify$OTU_id,OTU.classify$phylum))
#原始树、分类单元和bin的关系
ggtree(OTU.tree.groupByPhylum,layout = "circular",linewidth=.linewidth,branch.length = "none",
       aes(color=group)) %<+% OTU2bin+
  geom_tippoint(aes(color=as.character(bin.id.new)))

#获得只有bin的树
OTU_bin.tree = read.tree(file = paste0(PREFIX,".04bin.mafft.treefile"))
bin.classify = read.table(paste0(PREFIX,".03bin.classify.modified.tsv"))
rownames(bin.classify)=bin.classify$OTU_id
OTU_bin.tree$tip.label = bin.classify[OTU_bin.tree$tip.label,"bin_id"]
bin.classify = data.frame(bin_id=bin.classify$bin_id,bin.classify)

ggtree(OTU_bin.tree,layout = "fan",open.angle = 180,branch.length = "none") %<+%
  bin.classify +
  geom_tippoint(aes(colour = Taxa))

write.tree(OTU_bin.tree,file = paste0(PREFIX,".05bin.tree"))
#关闭



#整理每个bin的数据
PREFIX="bac"
contri = function(x){
  sum(x <= -0.975)*100/length(x)
}
source("../../public/public.R")
load(file = "bac.01iCAMPdata.Rdata")

# 这里考察的是每个bin的贡献，所以不用考虑bin的权重。
# 后续分析中考察不同类群，需要将bin加和，因此对bin要进行权重加和
sampleMeta = read.table("../../00data/sample_meta.tsv",sep="\t",header = T)
sampleMeta = fmtMeta(sampleMeta)

##mean by sex
sexMean.data = NULL
{
  M.sampleList = sampleMeta[sampleMeta$sex=='M',"sampleID"]
  M.sample.iCAMP.data = iCAMP.data[iCAMP.data$samp1 %in% M.sampleList & iCAMP.data$samp2 %in% M.sampleList,]
  M.contri.iCAMP.data = aggregate(x = M.sample.iCAMP.data$CbMPD,by=list(bin_id=M.sample.iCAMP.data$bin_id),contri)
  F.sampleList = sampleMeta[sampleMeta$sex=='F',"sampleID"]
  F.sample.iCAMP.data = iCAMP.data[iCAMP.data$samp1 %in% F.sampleList & iCAMP.data$samp2 %in% F.sampleList,]
  F.contri.iCAMP.data = aggregate(x = F.sample.iCAMP.data$CbMPD,by=list(bin_id=F.sample.iCAMP.data$bin_id),contri)
  if(!identical(M.contri.iCAMP.data$bin_id,F.contri.iCAMP.data$bin_id)){
    stop("not identical")
  }
  sexMean.data = rbind(sexMean.data,
                       data.frame(bin_id=M.contri.iCAMP.data$bin_id,
                                  M.hoS = M.contri.iCAMP.data$x,
                                  F.hoS = F.contri.iCAMP.data$x))
}
write.table(sexMean.data,file = paste0(PREFIX,".06-1sex.bin.hoS.tsv"),sep="\t",quote = F)

##mean by type sex
typesexMean.data = NULL
for(type in levels(sampleMeta$type)){
  M.sampleList = sampleMeta[sampleMeta$sex=='M' & sampleMeta$type==type,"sampleID"]
  M.sample.iCAMP.data = iCAMP.data[iCAMP.data$samp1 %in% M.sampleList & iCAMP.data$samp2 %in% M.sampleList,]
  M.contri.iCAMP.data = aggregate(x = M.sample.iCAMP.data$CbMPD,by=list(bin_id=M.sample.iCAMP.data$bin_id),contri)
  F.sampleList = sampleMeta[sampleMeta$sex=='F' & sampleMeta$type==type,"sampleID"]
  F.sample.iCAMP.data = iCAMP.data[iCAMP.data$samp1 %in% F.sampleList & iCAMP.data$samp2 %in% F.sampleList,]
  F.contri.iCAMP.data = aggregate(x = F.sample.iCAMP.data$CbMPD,by=list(bin_id=F.sample.iCAMP.data$bin_id),contri)
  if(!identical(M.contri.iCAMP.data$bin_id,F.contri.iCAMP.data$bin_id)){
    stop("not identical")
  }
  typesexMean.data = rbind(typesexMean.data,
                       data.frame(type=type,
                                  bin_id=M.contri.iCAMP.data$bin_id,
                                  M.hoS = M.contri.iCAMP.data$x,
                                  F.hoS = F.contri.iCAMP.data$x))
}
write.table(typesexMean.data,file = paste0(PREFIX,".06-2typesex.bin.hoS.tsv"),sep="\t",quote = F)


##mean by type sex period
typeperiodsexMean.data = NULL
for(type in levels(sampleMeta$type)){
for(period in levels(sampleMeta$period)){
  M.sampleList = sampleMeta[sampleMeta$sex=='M' & sampleMeta$type==type & sampleMeta$period==period,"sampleID"]
  M.sample.iCAMP.data = iCAMP.data[iCAMP.data$samp1 %in% M.sampleList & iCAMP.data$samp2 %in% M.sampleList,]
  M.contri.iCAMP.data = aggregate(x = M.sample.iCAMP.data$CbMPD,by=list(bin_id=M.sample.iCAMP.data$bin_id),contri)
  F.sampleList = sampleMeta[sampleMeta$sex=='F' & sampleMeta$type==type & sampleMeta$period==period,"sampleID"]
  F.sample.iCAMP.data = iCAMP.data[iCAMP.data$samp1 %in% F.sampleList & iCAMP.data$samp2 %in% F.sampleList,]
  F.contri.iCAMP.data = aggregate(x = F.sample.iCAMP.data$CbMPD,by=list(bin_id=F.sample.iCAMP.data$bin_id),contri)
  if(!identical(M.contri.iCAMP.data$bin_id,F.contri.iCAMP.data$bin_id)){
    stop("not identical")
  }
  typeperiodsexMean.data = rbind(typeperiodsexMean.data,
                           data.frame(type=type,
                                      period=period,
                                      bin_id=M.contri.iCAMP.data$bin_id,
                                      M.hoS = M.contri.iCAMP.data$x,
                                      F.hoS = F.contri.iCAMP.data$x))
}
}
write.table(typeperiodsexMean.data,file = paste0(PREFIX,".06-3typeperiodsex.bin.hoS.tsv"),sep="\t",quote = F)

#关闭



#计算每个bin的丰度
PREFIX = 'bac'
load(paste0(PREFIX,".01iCAMPdata.Rdata"))
OTU.abundance = read.table("../../00data/bacterial/feature_table_20norm.tsv")

bin_index = OTU2bin[rownames(OTU.abundance),"bin_id"]
sum(is.na(bin_index))
bin.abundance = NULL
for(sample in colnames(OTU.abundance)){
  sampleAbund = as.data.frame(tapply(OTU.abundance[,sample],bin_index,sum))
  colnames(sampleAbund)=sample
  if(class(bin.abundance)==class(NULL)){
    bin.abundance = sampleAbund
  }else{
    if(!identical(rownames(bin.abundance),rownames(sampleAbund))){
      stop('Error not match')
    }else{
      bin.abundance = cbind(bin.abundance,sampleAbund)
    }
  }
}
write.table(bin.abundance,file = paste0(PREFIX,".07-0bin.abundPerSample.tsv"),
            sep="\t",quote = F)

bin.mean.abundance = as.data.frame(apply(bin.abundance,1,mean))
colnames(bin.mean.abundance)='abundance.k'
write.table(bin.mean.abundance,file = paste0(PREFIX,".07-1bin.mean.abund.tsv"),
            sep="\t",quote = F)

# 关闭


PREFIX='bac'
library(ggplot2)
library(ggnewscale)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(tidyr)

source("../../public/public.R")
source("./treePlot.R")


bin_abund = read.table(paste0(PREFIX,".07-1bin.mean.abund.tsv"))

bin.type.hoS = read.table(paste0(PREFIX,".06-2typesex.bin.hoS.tsv"))

bin.type.hoS$max.hoS = apply(bin.type.hoS[c('M.hoS','F.hoS')],1,max)
bin.type.hoS$ratio = log2((bin.type.hoS$M.hoS+1)/(bin.type.hoS$F.hoS+1))

bin.type.hoS.wide = NULL
for (type in unique(bin.type.hoS$type)) {
  type.hoS = bin.type.hoS[bin.type.hoS$type==type,]
  colnames(type.hoS)[3:ncol(type.hoS)] = paste(type,colnames(type.hoS)[3:ncol(type.hoS)],sep = ".")
  if(class(bin.type.hoS.wide)==class(NULL)){
    bin.type.hoS.wide = type.hoS[-1]
  }else{
    if(!identical(bin.type.hoS.wide$bin_id,type.hoS$bin_id)){
      stop("Err not match")
    }
    bin.type.hoS.wide = cbind(bin.type.hoS.wide,type.hoS[c(-1,-2)])
  }
}
identical(bin.type.hoS.wide$bin_id,rownames(bin_abund))
bin.type.hoSAbund = cbind(bin.type.hoS.wide,bin_abund)


bin.tree=read.tree(paste0(PREFIX,".05bin.tree"))
bin.classify = read.table(paste0(PREFIX,".03bin.classify.modified.tsv"))
bin.type.hoSAbund$Taxa = bin.classify[bin.type.hoSAbund$bin_id,"Taxa"]

bin.tree.groupByTaxa = groupOTU(bin.tree,split(bin.classify$bin_id,bin.classify$Taxa))
taxaColor = c(Actinobacteria="#66C2A5",
              Alphaproteobacteria="#A6D854",
              Gammaproteobacteria="#E78AC3",
              Deltaproteobacteria="#FC8D62",
              Firmicutes="#8DA0CB",
              Planctomycetia="#FAC92F",
              Others="#B3B3B3")

ggtree(bin.tree.groupByTaxa,aes(colour = group),layout = "circular",open.angle = 180,
       ladderize = TRUE,branch.length = "none",linewidth=.linewidth)+
  geom_tiplab(size=.fontsizePT/ggplot2::.pt,show.legend=F)+
  scale_color_manual(values = taxaColor,name="Phylum/Class")+genelateTreetheme()
  
ggsave(paste0(PREFIX,'.09labelTree.pdf'),width=10,height=20,units = "cm")


MIN_SCALE = 8.2
# min(bin.type.hoSAbund$abundance.k) of fungi

bin.type.hoSAbund$abund = log10(bin.type.hoSAbund$abundance.k/MIN_SCALE)


# 为保证图形标尺统一，进行重新整理
bin.type.hoSAbund.Leaf.M = bin.type.hoSAbund[,c("bin_id","L.M.hoS","L.ratio")]
bin.type.hoSAbund.Leaf.F = bin.type.hoSAbund[,c("bin_id","L.F.hoS","L.ratio")]
colnames(bin.type.hoSAbund.Leaf.M)[2]=colnames(bin.type.hoSAbund.Leaf.F)[2]="L.hoS"
bin.type.hoSAbund.Leaf.F$L.hoS=-bin.type.hoSAbund.Leaf.F$L.hoS
bin.type.hoSAbund.Leaf = rbind(bin.type.hoSAbund.Leaf.M,bin.type.hoSAbund.Leaf.F)
bin.type.hoSAbund.Leaf$L.ratio = ifelse(bin.type.hoSAbund.Leaf$L.ratio >  2,  2.5, bin.type.hoSAbund.Leaf$L.ratio)
bin.type.hoSAbund.Leaf$L.ratio = ifelse(bin.type.hoSAbund.Leaf$L.ratio < -2, -2.5, bin.type.hoSAbund.Leaf$L.ratio)


bin.type.hoSAbund.Root.M = bin.type.hoSAbund[,c("bin_id","R.M.hoS","R.ratio")]
bin.type.hoSAbund.Root.F = bin.type.hoSAbund[,c("bin_id","R.F.hoS","R.ratio")]
colnames(bin.type.hoSAbund.Root.M)[2]=colnames(bin.type.hoSAbund.Root.F)[2]="R.hoS"
bin.type.hoSAbund.Root.F$R.hoS=-bin.type.hoSAbund.Root.F$R.hoS
bin.type.hoSAbund.Root = rbind(bin.type.hoSAbund.Root.M,bin.type.hoSAbund.Root.F)
bin.type.hoSAbund.Root$R.ratio = ifelse(bin.type.hoSAbund.Root$R.ratio >  2,  2.5, bin.type.hoSAbund.Root$R.ratio)
bin.type.hoSAbund.Root$R.ratio = ifelse(bin.type.hoSAbund.Root$R.ratio < -2, -2.5, bin.type.hoSAbund.Root$R.ratio)

bin.type.hoSAbund.Soil.M = bin.type.hoSAbund[,c("bin_id","S.M.hoS","S.ratio")]
bin.type.hoSAbund.Soil.F = bin.type.hoSAbund[,c("bin_id","S.F.hoS","S.ratio")]
colnames(bin.type.hoSAbund.Soil.M)[2]=colnames(bin.type.hoSAbund.Soil.F)[2]="S.hoS"
bin.type.hoSAbund.Soil.F$S.hoS=-bin.type.hoSAbund.Soil.F$S.hoS
bin.type.hoSAbund.Soil = rbind(bin.type.hoSAbund.Soil.M,bin.type.hoSAbund.Soil.F)
bin.type.hoSAbund.Soil$S.ratio = ifelse(bin.type.hoSAbund.Soil$S.ratio >  2,  2.5, bin.type.hoSAbund.Soil$S.ratio)
bin.type.hoSAbund.Soil$S.ratio = ifelse(bin.type.hoSAbund.Soil$S.ratio < -2, -2.5, bin.type.hoSAbund.Soil$S.ratio)

# 
# ggtree(bin.tree.groupByTaxa,aes(colour = group),layout = "fan",open.angle = 180,
#        ladderize = TRUE,branch.length = "none",linewidth=.linewidth)+
#   scale_color_manual(values = taxaColor,name="Phylum/Class")+
#   genelateTreetheme() + 
#   geom_fruit(data = bin.type.hoSAbund,geom=geom_col,
#              aes(y=bin_id,x=abund,fill=Taxa),
#              pwidth = 0.2,axis.params = list(axis="x",
#                                              text.size=.fontsize/3/ggplot2::.pt,
#                                              line.size=.linewidth/5))+
#   scale_fill_manual(values = taxaColor,name="Phylum/Class")+
#   ggnewscale::new_scale_fill()+
#   geom_fruit(data = bin.type.hoSAbund,geom=geom_col,
#              aes(y=bin_id,x=S.ratio,fill=S.max.hoS),
#              pwidth = 0.2,offset=0.1,axis.params = list(axis="xy",
#                                                         text.size=.fontsize/3/ggplot2::.pt,
#                                                         line.size=.linewidth/5))+
#   scale_fill_gradientn(colors = c(soilColor,soilColor,soilColor),limits=c(0,100),
#                        values = c(0,.5,1))+
#   ggnewscale::new_scale_fill()+
#   geom_fruit(data = bin.type.hoSAbund,geom=geom_col,
#              aes(y=bin_id,x=R.ratio,fill=R.max.hoS),
#              pwidth = 0.2,offset=0.1,axis.params = list(axis="xy",
#                                                         text.size=.fontsize/3/ggplot2::.pt,
#                                                         line.size=.linewidth/5))+
#   scale_fill_gradientn(colors = c(rootColor,rootColor,rootColor),limits=c(0,100),
#                        values = c(0,.5,1))+
#   ggnewscale::new_scale_fill()+
#   geom_fruit(data = bin.type.hoSAbund,geom=geom_col,
#              aes(y=bin_id,x=L.ratio,fill=L.max.hoS),
#              pwidth = 0.2,offset=0.1,axis.params = list(axis="xy",
#                                                         text.size=.fontsize/3/ggplot2::.pt,
#                                                         line.size=.linewidth/5))+
#   scale_fill_gradientn(colors = c(leafColor,leafColor,leafColor),limits=c(0,100),
#                        values = c(0,.5,1))+
#   theme(legend.position = "none")
# 


bac.plot = ggtree(bin.tree.groupByTaxa,aes(colour = group),layout = "circular",open.angle = 180,
       ladderize = TRUE,branch.length = "none",linewidth=.linewidth/2)+
  scale_color_manual(values = taxaColor,name="Phylum/Class")+
  genelateTreetheme() + 
  geom_fruit(data = bin.type.hoSAbund,geom=geom_col,
             aes(y=bin_id,x=abund,fill=Taxa),
             pwidth = 0.2,axis.params = list(axis="x",
                                             text.size=.fontsize/3/ggplot2::.pt,
                                             line.size=.linewidth/5,
                                             limits=c(0,3.5)))+
  scale_fill_manual(values = taxaColor,name="Phylum/Class")+
  ggnewscale::new_scale_fill()+
  geom_fruit(data = bin.type.hoSAbund.Soil,geom=geom_col,
             aes(y=bin_id,x=S.hoS,fill=S.ratio),
             pwidth = 0.3782461,offset=0.3782461,
             axis.params = list(axis="xy",
                                text.size=.fontsize/3/ggplot2::.pt,
                                line.size=.linewidth/5,
                                limits=c(-100,100)))+
  scale_fill_gradientn(colors = alpha(soilColor,c(1,1,0.05,1,1)),limits=c(-3,3),
                       values = scales::rescale(c(-3,-2,0,2,3),from=range(-3,3)),
                       breaks=c(-3,-2,0,2,3),labels=c('',-2,0,2,''))+
  ggnewscale::new_scale_fill()+
  geom_fruit(data = bin.type.hoSAbund.Root,geom=geom_col,
             aes(y=bin_id,x=R.hoS,fill=R.ratio),
             pwidth = 0.3063547,offset=0.3063547,
             axis.params = list(axis="xy",
                                text.size=.fontsize/3/ggplot2::.pt,
                                line.size=.linewidth/5,
                                limits=c(-75,85)))+
  scale_fill_gradientn(colors = alpha(rootColor,c(1,1,0.05,1,1)),limits=c(-3,3),
                       values = scales::rescale(c(-3,-2,0,2,3),from=range(-3,3)),
                       breaks=c(-3,-2,0,2,3),labels=c('',-2,0,2,''))+
  ggnewscale::new_scale_fill()+
  geom_fruit(data = bin.type.hoSAbund.Leaf,geom=geom_col,
             aes(y=bin_id,x=L.hoS,fill=L.ratio),
             pwidth = 0.2,offset=0.2,
             axis.params = list(axis="xy",
                                text.size=.fontsize/3/ggplot2::.pt,
                                line.size=.linewidth/5,
                                limits=c(-45,55)))+
  scale_fill_gradientn(colors = alpha(leafColor,c(1,1,0.05,1,1)),limits=c(-3,3),
                       values = scales::rescale(c(-3,-2,0,2,3),from=range(-3,3)),
                       breaks=c(-3,-2,0,2,3),labels=c('',-2,0,2,''))+
  theme(legend.key.size = unit(.symbolSizePT*3,'pt'))
  

bac.plot


save(bac.plot,file="bac.19.plot.Rdata")


#关闭
PREFIX='bac'
library(ggplot2)
source("../../public/public.R")

bin.classify = read.table(paste0(PREFIX,".03bin.classify.modified.tsv"))
load(file = paste0(PREFIX,".01iCAMPdata.Rdata"))
sampleMeta = read.table("../../00data/sample_meta.tsv",sep="\t",header = T)
sampleMeta = fmtMeta(sampleMeta)

iCAMP.data$weighted.CbMPD = ifelse(iCAMP.data$CbMPD <= -0.975,iCAMP.data$bin.weight,0)
iCAMP.data$Taxa = bin.classify[iCAMP.data$bin_id,"Taxa"]

iCAMP.HoS.byTaxa = aggregate.data.frame(iCAMP.data$weighted.CbMPD,
                                         list(samp1=iCAMP.data$samp1,samp2=iCAMP.data$samp2,
                                              Taxa=iCAMP.data$Taxa),
                                         sum)
iCAMP.weight.byTaxa = aggregate.data.frame(iCAMP.data$bin.weight,
                                        list(samp1=iCAMP.data$samp1,samp2=iCAMP.data$samp2,
                                             Taxa=iCAMP.data$Taxa),
                                        sum)
iCAMP.HoS.byTaxa.frame = data.frame(
  iCAMP.HoS.byTaxa[1:3],
  HoS = iCAMP.HoS.byTaxa$x*100, 
  
  Taxa.weight = iCAMP.weight.byTaxa$x,
  HoS.within.Taxa = ifelse(iCAMP.HoS.byTaxa$x==0,0,iCAMP.HoS.byTaxa$x*100 / iCAMP.weight.byTaxa$x)
)


##mean by type

HoS.byTaxa.typeMean = NULL
for(type in levels(sampleMeta$type)){
  sampleList = sampleMeta[sampleMeta$type==type,"sampleID"]
  sample.HoS.byTaxa = iCAMP.HoS.byTaxa.frame[iCAMP.HoS.byTaxa.frame$samp1 %in% sampleList & 
                                               iCAMP.HoS.byTaxa.frame$samp2 %in% sampleList,]
  sample.HoS.byTaxa.typeMean = aggregate.data.frame(x=sample.HoS.byTaxa$HoS,
                                                    by = list(Taxa = sample.HoS.byTaxa$Taxa),
                                                    mean)
  sample.HoS.byTaxa.typeMean.within = aggregate.data.frame(x=sample.HoS.byTaxa$HoS.within.Taxa,
                                                    by = list(Taxa = sample.HoS.byTaxa$Taxa),
                                                    mean)
  HoS.byTaxa.typeMean = rbind(HoS.byTaxa.typeMean,
                           data.frame(type=type,
                                      Taxa=sample.HoS.byTaxa.typeMean$Taxa,
                                      HoS = sample.HoS.byTaxa.typeMean$x,
                                      HoS.within = sample.HoS.byTaxa.typeMean.within$x,
                                      row.names = NULL))
}

taxaColor = c(Actinobacteria="#66C2A5",
              Alphaproteobacteria="#A6D854",
              Gammaproteobacteria="#E78AC3",
              Deltaproteobacteria="#FC8D62",
              Firmicutes="#8DA0CB",
              Planctomycetia="#FAC92F",
              Others="#B3B3B3")

TaxaSum = tapply(HoS.byTaxa.typeMean$HoS,HoS.byTaxa.typeMean$type,sum)
HoS.byTaxa.typeMean$HoS.pro = HoS.byTaxa.typeMean$HoS *100 / (TaxaSum[HoS.byTaxa.typeMean$type])

HoS.byTaxa.typeMean$Taxa = factor(HoS.byTaxa.typeMean$Taxa,levels = names(taxaColor),ordered = T)

# Relative contributions to homogeneous selection(HoS)  of different taxa 
# in Leaf(L),Root(R),and Rhizosphere soil(S)
# weighted by bin abundance

hoS.plot = ggplot(data=HoS.byTaxa.typeMean)+
  geom_bar(aes(x = type,y = HoS.pro,fill = Taxa),stat = "identity",show.legend = F)+
  scale_fill_manual(values = taxaColor)+
  theme_bw()+genelateBWtheme()+
  scale_y_continuous(name="Relative contributions to homogeneous selection (%)",expand = c(0,0))+
  theme(axis.ticks.length.x = unit(0,'pt'),axis.ticks.x = element_blank(),axis.title.x = element_blank())


##丰度绘图
bin.abund = read.table(paste0(PREFIX,".07-0bin.abundPerSample.tsv"))

bin.abundByType = NULL
Taxa.abundByType = NULL
for (type in levels(sampleMeta$type)) {
  type.bin.abund = bin.abund[,sampleMeta[sampleMeta$type==type,"sampleID"]]
  type.bin.abund.mean = apply(type.bin.abund,1,mean)
  
  bin.abundByType = rbind(bin.abundByType,data.frame(
    type=type,
    bin_id = names(type.bin.abund.mean),
    abund = type.bin.abund.mean/1e3,row.names = NULL
  ))
  
  type.taxa.abund = tapply(type.bin.abund.mean,bin.classify[names(type.bin.abund.mean),"Taxa"],FUN = sum)
  Taxa.abundByType = rbind(Taxa.abundByType,data.frame(
    type=type,
    Taxa = names(type.taxa.abund),
    abund = type.taxa.abund/1e3,row.names = NULL
  ))
}

Taxa.abundByType$Taxa = factor(Taxa.abundByType$Taxa,levels = names(taxaColor),ordered = T)

identical(HoS.byTaxa.typeMean[c(1,2)],Taxa.abundByType[c(1,2)])
HoS.byTaxa.typeMean$abund = Taxa.abundByType$abund

abund.plot = ggplot(data=HoS.byTaxa.typeMean)+
  geom_bar(aes(x = type,y = abund,fill = Taxa),stat = "identity",show.legend = F)+
  scale_fill_manual(values = taxaColor)+
  theme_bw()+genelateBWtheme()+
  scale_y_continuous(name="Relative abundance (%)",expand = c(0,0))+
  theme(axis.ticks.length.x = unit(0,'pt'),axis.ticks.x = element_blank(),axis.title.x = element_blank())


write.table(HoS.byTaxa.typeMean,file=paste0(PREFIX,".08-1TaxaAbundHoSPercent.tsv"),quote = F,sep = "\t")

##increase/Decrease of sex hoS
# weighted by bin abundance

# 
##weighted mean by type sex
typesexMean.w.data = NULL
for(type in levels(sampleMeta$type)){
  M.sampleList = sampleMeta[sampleMeta$sex=='M' & sampleMeta$type==type,"sampleID"]
  M.sample.iCAMP.data = iCAMP.data[iCAMP.data$samp1 %in% M.sampleList & iCAMP.data$samp2 %in% M.sampleList,]
  M.contri.iCAMP.data = aggregate(x = M.sample.iCAMP.data$weighted.CbMPD,by=list(bin_id=M.sample.iCAMP.data$bin_id),sum)
  F.sampleList = sampleMeta[sampleMeta$sex=='F' & sampleMeta$type==type,"sampleID"]
  F.sample.iCAMP.data = iCAMP.data[iCAMP.data$samp1 %in% F.sampleList & iCAMP.data$samp2 %in% F.sampleList,]
  F.contri.iCAMP.data = aggregate(x = F.sample.iCAMP.data$weighted.CbMPD,by=list(bin_id=F.sample.iCAMP.data$bin_id),sum)
  if(!identical(M.contri.iCAMP.data$bin_id,F.contri.iCAMP.data$bin_id)){
    stop("not identical")
  }
  typesexMean.w.data = rbind(typesexMean.w.data,
                           data.frame(type=type,
                                      bin_id=M.contri.iCAMP.data$bin_id,
                                      M.hoS = M.contri.iCAMP.data$x,
                                      F.hoS = F.contri.iCAMP.data$x))
}
write.table(typesexMean.w.data,file = paste0(PREFIX,".06-8typesex.bin.weightedhoS.tsv"),sep="\t",quote = F)



typesexMean.HoS = typesexMean.w.data
rownames(typesexMean.HoS)=NULL
typesexMean.HoS$MF.diff = typesexMean.HoS$M.hoS-typesexMean.HoS$F.hoS
typesexMean.HoS$Group = ifelse(typesexMean.HoS$MF.diff>0,'M','F')

identical(bin.abundByType[,c(1,2)],typesexMean.HoS[,c(1,2)])

typesexMean.HoS$abund = bin.abundByType$abund
typesexMean.HoS$weighted.MF.diff = typesexMean.HoS$MF.diff * typesexMean.HoS$abund/1e2
typesexMean.HoS$Taxa = bin.classify[typesexMean.HoS$bin_id,"Taxa"]

Taxa.typesex.HoS = aggregate.data.frame(typesexMean.HoS$weighted.MF.diff,
                                        by=list(type=typesexMean.HoS$type,
                                                Taxa=typesexMean.HoS$Taxa,
                                                Group=typesexMean.HoS$Group),
                                        sum)
colnames(Taxa.typesex.HoS)[4]='weighted.MF.diff'
Taxa.typesex.HoS$Taxa = factor(Taxa.typesex.HoS$Taxa,levels = names(taxaColor),ordered = T)


sexDiff.plot = ggplot(data=Taxa.typesex.HoS)+
  geom_bar(aes(x = type,y = weighted.MF.diff,fill = Taxa),stat = "identity",show.legend = F)+
  facet_wrap(~type,scales = "free")+
  scale_fill_manual(values = taxaColor)+
  theme_bw()+genelateBWtheme()+
  scale_y_continuous(name="Sex-induced difference (%)",expand = c(0,0),limits = function(o){max(abs(o))*c(-1,1)})+
  scale_x_discrete(expand = c(0,0))+
  geom_hline(yintercept = 0,linewidth=.linewidth)+
  theme(axis.ticks.length.x = unit(0,'pt'),axis.ticks.x = element_blank(),axis.title.x = element_blank(),
        strip.background = element_blank(),strip.text = element_blank())

write.table(Taxa.typesex.HoS,file = paste0(PREFIX,".08-2TaxaHosSexDiff.tsv"))

bac.abund.plot = abund.plot;
bac.hos.plot=hoS.plot;
bac.sexDiff=sexDiff.plot
save(bac.abund.plot,bac.hos.plot,bac.sexDiff,
     file=paste0(PREFIX,".29StatPlot.Rdata"))


