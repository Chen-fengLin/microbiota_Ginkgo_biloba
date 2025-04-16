source("../../public/public.R")
library(phylosignal)
library(phylobase)
library(ape)

bac.tree = read.tree("../../00data/bacterial/bac.tree")
bac.tree$node.label = NULL

bac.signal=NULL
for(type in c('leaf','root','soil')){
  for(compare in c('S1S2','S2S3')){
    diffres = read.table(paste0("../41Diff_abund_period/DiffRes/bac.",type,".",compare,"diff.tsv"))
    clean.tree = keep.tip(bac.tree,tip = rownames(diffres))
    phyloTree = phylo4d(clean.tree,tip.data=as.data.frame(diffres['LFC']))
    
    result=phyloSignal(phyloTree,reps=999)
    
    bac.signal = rbind(bac.signal,data.frame(Niche=type,compare=compare,
                                             result$stat,result$pvalue,star(result$pvalue)))
  }
}


fgi.tree = read.tree("../../00data/fungi/fgi.tree")
fgi.tree$node.label = NULL

fgi.signal=NULL
for(type in c('leaf','root','soil')){
  for(compare in c('S1S2','S2S3')){
    diffres = read.table(paste0("../41Diff_abund_period/DiffRes/fgi.",type,".",compare,"diff.tsv"))
    clean.tree = keep.tip(fgi.tree,tip = rownames(diffres))
    phyloTree = phylo4d(clean.tree,tip.data=as.data.frame(diffres['LFC']))
    
    result=phyloSignal(phyloTree,reps=999)
    
    fgi.signal = rbind(fgi.signal,data.frame(Niche=type,compare=compare,
                                             result$stat,result$pvalue,star(result$pvalue)))
  }
}

signal.all = rbind(bac.signal,fgi.signal)

write.table(signal.all,"bac_fgi.signal.tsv",sep="\t",row.names = F,quote = F)
