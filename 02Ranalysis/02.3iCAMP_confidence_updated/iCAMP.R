library(iCAMP)
library(ape)
library(vegan)
NWORKER=64

##bacterial
bac.otu = as.data.frame(t(read.table("../../00data/bacterial/feature_table_10clean.tsv")))
bac.otu = rrarefy(bac.otu,sample=min(rowSums(bac.otu)))
bac.tree = read.tree("../../00data/bacterial/bac.tree")

bac.icamp.res = icamp.big(bac.otu,bac.tree,nworker = NWORKER,detail.null = TRUE,bin.size.limit=24,
                          prefix="bac.iCAMP",pd.wd=paste0(tempdir(),"/bac.pdbig.icampbig"))
save(bac.icamp.res,file="bac.icamp.res.Rdata")
write.table(bac.icamp.res$CbMPDiCBraya,file = "bac.icamp.res.tsv",sep="\t",row.names = F,quote = F)

##fungi
fgi.otu = as.data.frame(t(read.table("../../00data/fungi/feature_table_10clean.tsv")))
fgi.otu = rrarefy(fgi.otu,sample=min(rowSums(fgi.otu)))
fgi.tree = read.tree("../../00data/fungi/fgi.tree")

fgi.icamp.res = icamp.big(fgi.otu,fgi.tree,nworker = NWORKER,detail.null = TRUE,bin.size.limit=12,
                          prefix="fgi.iCAMP",pd.wd=paste0(tempdir(),"/fgi.pdbig.icampbig"))
save(fgi.icamp.res,file="fgi.icamp.res.Rdata")
write.table(fgi.icamp.res$CbMPDiCBraya,file = "fgi.icamp.res.tsv",sep="\t",row.names = F,quote = F)


