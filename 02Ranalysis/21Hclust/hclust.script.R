library(ape)

##Bray-Curtis
#bacterial
loadS = load("../20Distance_matrix/bacterial_BC_dismatrix.Rdata")
bac.dismatrix = eval(as.symbol(loadS))
bac.cluster = hclust(bac.dismatrix,method = 'ward.D2')
bac.treedata <- as.phylo(bac.cluster)# 将聚类结果转成系统发育格式
write.tree(phy=bac.treedata, file="bacterial_BC_tree.nwk") # 输出newick格式文件

#fungi
loadS = load("../20Distance_matrix/fungi_BCdismatrix.Rdata")
fgi.dismatrix = eval(as.symbol(loadS))
fgi.cluster = hclust(fgi.dismatrix,method = 'ward.D2')
fgi.treedata <- as.phylo(fgi.cluster)# 将聚类结果转成系统发育格式
write.tree(phy=fgi.treedata, file="fungi_BC_tree.nwk") # 输出newick格式文件


