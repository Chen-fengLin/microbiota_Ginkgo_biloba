library(vegan)


##bacterial
bac.alldata = read.table(file="../../00data/bacterial/feature_table_10clean.tsv")
bac.normlized.data = bac.alldata
for (sampleName in colnames(bac.alldata)) {
  sample.sum = sum(bac.alldata[sampleName])
  bac.normlized.data[sampleName] = bac.alldata[sampleName]*1e5/sample.sum
}
write.table(bac.normlized.data,file="../../00data/bacterial/feature_table_20norm.tsv",
            sep="\t",quote = FALSE)
bac.distmatrix <- vegdist(t(bac.normlized.data),method="bray")

save(bac.distmatrix,file="bacterial_BC_dismatrix.Rdata")

##fungi

fgi.alldata = read.table(file="../../00data/fungi/feature_table_10clean.tsv")
fgi.normlized.data = fgi.alldata
for (sampleName in colnames(fgi.alldata)) {
  sample.sum = sum(fgi.alldata[sampleName])
  fgi.normlized.data[sampleName] = fgi.alldata[sampleName]*1e5/sample.sum
}
write.table(fgi.normlized.data,file="../../00data/fungi/feature_table_20norm.tsv",
            sep="\t",quote = FALSE)
fgi.distmatrix <- vegdist(t(fgi.normlized.data),method="bray")

save(fgi.distmatrix,file="fungi_BCdismatrix.Rdata")
