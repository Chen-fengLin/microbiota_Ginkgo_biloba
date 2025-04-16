source("../01enrich_analysis/enrich.R")
source("../../public/public.R")
library(phytools)

arrangeEnrich <-function(diff,all,classify,func,notes){
  if(!length(diff))return(NULL)
  res = rbind(
    calculateTaxonEnrich(diff,all,classify,levl='phylum'),
    calculateTaxonEnrich(diff,all,classify,levl='class'),
    calculateTaxonEnrich(diff,all,classify,levl='order'),
    calculateTaxonEnrich(diff,all,classify,levl='family'),
    calculateTaxonEnrich(diff,all,classify,levl='genus')
    # calculateEnrich(diff,all,func)
  )
  if(nrow(res)==0)return(NULL)
  data.frame(notes=notes,res)
}


bacterial_function_table = read.table("../../00data/bacterial/function_table_10clean.tsv")
bacterial_classify_table = read.table("../../00data/bacterial/classified_table_10clean.tsv",sep="\t")
bac.tree = read.tree("../../00data/bacterial/bac.tree")


bac.res=NULL
bac.signal=NULL
for(type in c('leaf','root','soil')){
  for(compare in c('S1S2','S2S3')){
    diffres = read.table(paste0("../41Diff_abund_period/DiffRes/bac.",type,".",compare,"diff.tsv"))
    K.test = phylosig(bac.tree,diffres[bac.tree$tip.label,'LFC'],method ='K' ,test = T)
    lambda.test = phylosig(bac.tree,diffres[bac.tree$tip.label,'LFC'],method ='lambda' ,test = T)
    all.OTU = row.names(diffres)
    diff.up.OTU = all.OTU[which(diffres$p_value<=0.05 & diffres$LFC >= 1)]
    diff.dn.OTU = all.OTU[which(diffres$p_value<=0.05 & diffres$LFC <= -1)]
    diff.OTU = c(diff.up.OTU,diff.dn.OTU)
    # up.enrich = arrangeEnrich(diff.up.OTU,all.OTU,bacterial_classify_table,bacterial_function_table,paste(type,compare,'up',sep = "-"))
    # dn.enrich = arrangeEnrich(diff.dn.OTU,all.OTU,bacterial_classify_table,bacterial_function_table,paste(type,compare,'dn',sep = "-"))
    both.enrich = arrangeEnrich(diff.OTU,all.OTU,bacterial_classify_table,bacterial_function_table,paste(type,compare,'both',sep = "-"))
    bac.res = rbind(bac.res,both.enrich)
    bac.signal = rbind(bac.signal,data.frame(Niche=type,compare=compare,K=K.test$K,P_K=K.test$P,star_K=star(K.test$P),
                                lambda=lambda.test$lambda,
                                P_lambda=lambda.test$P,
                                star_lambda = star(lambda.test$P)))
  }
}


bac.res.signif = bac.res[bac.res$p<=0.05,]
write.table(bac.res.signif,file = "bac.enrich.signif.tsv",sep = "\t",quote = F,row.names = F)
write.table(bac.signal,file = "bac.deve.signal.tsv",sep="\t",quote = F,row.names = F)


fungi_function_table = read.table("../../00data/fungi/function_table_10clean.tsv")
fungi_classify_table = read.table("../../00data/fungi/classified_table_10clean.tsv",sep="\t")
fgi.tree = read.tree("../../00data/fungi/fgi.tree")

fgi.res=NULL
fgi.signal=NULL
for(type in c('leaf','root','soil')){
  for(compare in c('S1S2','S2S3')){
    diffres = read.table(paste0("../41Diff_abund_period/DiffRes/fgi.",type,".",compare,"diff.tsv"))
    K.test = phylosig(fgi.tree,diffres[fgi.tree$tip.label,'LFC'],method ='K' ,test = T)
    lambda.test = phylosig(fgi.tree,diffres[fgi.tree$tip.label,'LFC'],method ='lambda' ,test = T)
    all.OTU = row.names(diffres)
    diff.up.OTU = all.OTU[which(diffres$p_value<=0.05 & diffres$LFC >= 1)]
    diff.dn.OTU = all.OTU[which(diffres$p_value<=0.05 & diffres$LFC <= -1)]
    diff.OTU = c(diff.up.OTU,diff.dn.OTU)
    # up.enrich = arrangeEnrich(diff.up.OTU,all.OTU,fungi_classify_table,fungi_function_table,paste(type,compare,'up',sep = "-"))
    # dn.enrich = arrangeEnrich(diff.dn.OTU,all.OTU,fungi_classify_table,fungi_function_table,paste(type,compare,'dn',sep = "-"))
    both.enrich = arrangeEnrich(diff.OTU,all.OTU,fungi_classify_table,fungi_function_table,paste(type,compare,'both',sep = "-"))
    fgi.res = rbind(fgi.res,both.enrich)
    fgi.signal = rbind(fgi.signal,data.frame(Niche=type,compare=compare,K=K.test$K,P_K=K.test$P,star_K=star(K.test$P),
                                             lambda=lambda.test$lambda,
                                             P_lambda=lambda.test$P,
                                             star_lambda = star(lambda.test$P)))
  }
}

fgi.res.signif = fgi.res[fgi.res$p<=0.05,]
write.table(fgi.res.signif,file = "fgi.enrich.signif.tsv",sep = "\t",quote = F,row.names = F)
write.table(fgi.signal,file = "fgi.deve.signal.tsv",sep="\t",quote = F,row.names = F)

