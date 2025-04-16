library(vegan)
source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)

getCleanData <- function(rawdata){
  cleandata = rawdata[which(rowSums(rawdata)>0),]
  return(t(cleandata))
}

calculatePERMANOVA<- function(y,myformula,...,method="bray",permutations = 999,by = "terms",seed=2024){
  #PERMANOVA教程  https://blog.csdn.net/qazplm12_3/article/details/120620477
  
  meta.clean = sampleMeta[rownames(y),]
  ## 这里默认进行贯序检验sequential test,变量顺序有意义，依次添加
  set.seed(seed)
  result0 = adonis2(myformula,meta.clean,method=method,permutations=permutations,by=by,...)
  result.frame = as.data.frame(result0)
  colnames(result.frame) = c('Df','SumOfSqs','R2','F','Pr')
  result.frame$signif = with(result.frame,ifelse(Pr<=0.001,"***",ifelse(Pr<=0.01,"**",
                                                                        ifelse(Pr<0.05,"*",ifelse(Pr<0.1,"."," ")))))
  
  ##检验离散度
  # result.frame$disper_P = NA
  # y.dist <- vegdist(y,method = method,binary = F)
  # for(item in row.names(result.frame)[1:(nrow(result.frame)-2)]){
  #   if(item %in% colnames(meta.clean)){
  #     item.dispersion <- betadisper(y.dist,group=meta.clean[[item]])
  #     item.disp_test = permutest(item.dispersion)
  #     item.pvalue = item.disp_test$tab$`Pr(>F)`[1]
  #     result.frame[item,'disper_P']=item.pvalue
  #   }
  # }
  return(result.frame)
}

bac.alldata = read.table(file="../../00data/bacterial/feature_table_20norm.tsv")
fgi.alldata = read.table(file = "../../00data/fungi/feature_table_20norm.tsv")
bac.sampleMeta = sampleMeta[colnames(bac.alldata),]
fgi.sampleMeta = sampleMeta[colnames(fgi.alldata),]
# ##all sample
# alldata.clean = getCleanData(alldata)
# 
# alldata.permanova = calculatePERMANOVA(alldata.clean,alldata.clean~type*period*sex)
# alldata.permanova
# 
# save(bac.alldata.NMDS.plot3,file = "./NMDS_bac/00alldata.plot.Rdata")


res=NULL
for(type in unique(bac.sampleMeta$typeByName)){
  typedata.clean = getCleanData(bac.alldata[,which(bac.sampleMeta$typeByName==type)])
  
  typedata.permanova = calculatePERMANOVA(typedata.clean,typedata.clean~period*sex)
  total = typedata.permanova['Total','SumOfSqs']
  res=rbind(res,data.frame(type=type,
                           Variables=rownames(typedata.permanova),
                           df=typedata.permanova$Df,
                           Fvalue = typedata.permanova$F,
                           SumOfSqs=typedata.permanova$SumOfSqs,
                           contri=typedata.permanova$SumOfSqs/total,P=typedata.permanova$Pr)[1:3,])
}

for(type in unique(fgi.sampleMeta$typeByName)){
  typedata.clean = getCleanData(fgi.alldata[,which(fgi.sampleMeta$typeByName==type)])
  
  typedata.permanova = calculatePERMANOVA(typedata.clean,typedata.clean~period*sex)
  total = typedata.permanova['Total','SumOfSqs']
  res=rbind(res,data.frame(type=type,
                           Variables=rownames(typedata.permanova),
                           df=typedata.permanova$Df,
                           Fvalue = typedata.permanova$F,
                           SumOfSqs=typedata.permanova$SumOfSqs,
                           contri=typedata.permanova$SumOfSqs/total,P=typedata.permanova$Pr)[1:3,])
}
write.table(res,"permanova.bac+fgi.tsv",sep="\t",quote = F,row.names = F)
