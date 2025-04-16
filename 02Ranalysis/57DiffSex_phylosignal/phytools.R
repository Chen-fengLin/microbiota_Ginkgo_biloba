source("../../public/public.R")
library(phytools)


bac.tree = read.tree("../../00data/bacterial/bac.tree")

bac.signal=NULL
for(type in c('L','R','S')){
  for(stage in c('S1','S2','S3')){
    print(paste(type,stage))
    diffres.total = read.table(paste0("../52Diff_abund_sex/Diff_res_bac/",type,".",stage,".sexDiff.tsv"))
    diffres = diffres.total[apply(diffres.total[,c('F.Mean','M.Mean')],1,max)>=0.1,]
    K.test = phylosig(bac.tree,diffres[bac.tree$tip.label,'LFC'],method ='K' ,test = T)
    plot(K.test)
    lambda.test = phylosig(bac.tree,diffres[bac.tree$tip.label,'LFC'],method ='lambda' ,test = T)
    bac.signal = rbind(bac.signal,data.frame(Niche=type,stage=stage,K=K.test$K,P_K=K.test$P,star_K=star(K.test$P),
                                lambda=lambda.test$lambda,
                                P_lambda=lambda.test$P,
                                star_lambda = star(lambda.test$P)))
  }
}


fgi.tree = read.tree("../../00data/fungi/fgi.tree")

fgi.signal=NULL
for(type in c('L','R','S')){
  for(stage in c('S1','S2','S3')){
    print(paste(type,stage))
    diffres.total = read.table(paste0("../52Diff_abund_sex/Diff_res_fgi/",type,".",stage,".sexDiff.tsv"))
    diffres = diffres.total[apply(diffres.total[,c('F.Mean','M.Mean')],1,max)>=0.1,]
    K.test = phylosig(fgi.tree,diffres[fgi.tree$tip.label,'LFC'],method ='K' ,test = T)
    plot(K.test)
    lambda.test = phylosig(fgi.tree,diffres[fgi.tree$tip.label,'LFC'],method ='lambda' ,test = T)
    fgi.signal = rbind(fgi.signal,data.frame(Niche=type,stage=stage,K=K.test$K,P_K=K.test$P,star_K=star(K.test$P),
                                             lambda=lambda.test$lambda,
                                             P_lambda=lambda.test$P,
                                             star_lambda = star(lambda.test$P)))
  }
}

signal = rbind(bac.signal,fgi.signal)
write.table(signal,file = "bac+fgi.signal.tsv",sep="\t",quote = F,row.names = F)

