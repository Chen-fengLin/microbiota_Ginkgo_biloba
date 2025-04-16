source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",sep = "\t",header = T)
sampleMeta=fmtMeta(sampleMeta)
alpha_diversity=read.table("../32alpha_diversity/alpha_diversity.tsv",sep="\t")

bac.alpha_diversity = alpha_diversity[alpha_diversity$OTUtype=='bacterial',]
bac.total.frame = NULL
for(params in unique(bac.alpha_diversity$params)){
  param.frame = NULL
  bac.params.alpha_diversity=alpha_diversity[alpha_diversity$OTUtype=='bacterial' & alpha_diversity$params==params,c("sampleID",'value')]
  bac.params.alpha_diversity.combind = data.frame(bac.params.alpha_diversity,sampleMeta[bac.params.alpha_diversity$sampleID,])
  total.res=as.data.frame(anova(lm(1+value~type*period*sex,data = bac.params.alpha_diversity.combind)))[1:3,]
  param.frame = data.frame(otu="bacterial",Samples='All samples',Variables=rownames(total.res),
                           Fvalue=total.res$`F value`,df=total.res$Df,P=total.res$`Pr(>F)`)
  for(type in unique(bac.params.alpha_diversity.combind$typeByName)){
    bac.params.alpha_diversity.combind.type = bac.params.alpha_diversity.combind[bac.params.alpha_diversity.combind$typeByName==type,]
    type.res=as.data.frame(anova(lm(1+value~period*sex,data = bac.params.alpha_diversity.combind.type)))[1:2,]
    param.frame = rbind(param.frame,data.frame(otu="bacterial",Samples=type,Variables=rownames(type.res),
                                               Fvalue=type.res$`F value`,df=type.res$Df,P=type.res$`Pr(>F)`))
  }
  colnames(param.frame)[4:6]<-paste(params,colnames(param.frame)[4:6],sep="_")
  if(class(bac.total.frame)=='NULL'){
    bac.total.frame=param.frame
  }else{
    bac.total.frame=cbind(bac.total.frame,param.frame[4:6])
  }
}


fgi.alpha_diversity = alpha_diversity[alpha_diversity$OTUtype=='fungi',]
fgi.total.frame = NULL
for(params in unique(fgi.alpha_diversity$params)){
  param.frame = NULL
  fgi.params.alpha_diversity=alpha_diversity[alpha_diversity$OTUtype=='fungi' & alpha_diversity$params==params,c("sampleID",'value')]
  fgi.params.alpha_diversity.combind = data.frame(fgi.params.alpha_diversity,sampleMeta[fgi.params.alpha_diversity$sampleID,])
  total.res=as.data.frame(anova(lm(1+value~type*period*sex,data = fgi.params.alpha_diversity.combind)))[1:3,]
  param.frame = data.frame(otu="fungi",Samples='All samples',Variables=rownames(total.res),
                           Fvalue=total.res$`F value`,df=total.res$Df,P=total.res$`Pr(>F)`)
  for(type in unique(fgi.params.alpha_diversity.combind$typeByName)){
    fgi.params.alpha_diversity.combind.type = fgi.params.alpha_diversity.combind[fgi.params.alpha_diversity.combind$typeByName==type,]
    type.res=as.data.frame(anova(lm(1+value~period*sex,data = fgi.params.alpha_diversity.combind.type)))[1:2,]
    param.frame = rbind(param.frame,data.frame(otu="fungi",Samples=type,Variables=rownames(type.res),
                                               Fvalue=type.res$`F value`,df=type.res$Df,P=type.res$`Pr(>F)`))
  }
  colnames(param.frame)[4:6]<-paste(params,colnames(param.frame)[4:6],sep="_")
  if(class(fgi.total.frame)=='NULL'){
    fgi.total.frame=param.frame
  }else{
    fgi.total.frame=cbind(fgi.total.frame,param.frame[4:6])
  }
}

write.table(rbind(bac.total.frame,fgi.total.frame),file = "alpha.linear.tsv",sep="\t",quote = F,row.names = F)
