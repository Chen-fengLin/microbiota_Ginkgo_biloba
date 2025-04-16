library(ggplot2)
library(ggtree)
# ggtree v3.6.2 For help: https://yulab-smu.top/treedata-book/
#   
#   If you use the ggtree package suite in published research, please cite
# the appropriate paper(s):
#   
#   Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
# ggtree: an R package for visualization and annotation of phylogenetic
# trees with their covariates and other associated data. Methods in
# Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628
# 
# Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods
# for mapping and visualizing associated data on phylogeny using ggtree.
# Molecular Biology and Evolution. 2018, 35(12):3041-3043.
# doi:10.1093/molbev/msy194
# 
# Guangchuang Yu. Using ggtree to visualize data on tree-like structures.
# Current Protocols in Bioinformatics. 2020, 69:e96. doi:10.1002/cpbi.96
source("../../public/public.R")


sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)
library(ape)

##bacterial
bac.treedata = read.tree("bacterial_BC_tree.nwk")
typeByNameGroup = split(row.names(sampleMeta),sampleMeta$typeByName)
bac.treedata = groupOTU(bac.treedata,typeByNameGroup)
colnames(sampleMeta)[1]='taxa'
bac.tree_plot = ggtree(bac.treedata,layout = "circular",linewidth=.linewidth,ladderize = F,
       aes(color=group)) %<+% sampleMeta +
  geom_tippoint(aes(shape=periodSexByName,color=typeByName),
                fill="white",size=.symbolSize*0.8)+
  scale_color_manual(name="",values=c('transparent',colorList),
                     labels=c('','Leaf','Root','Rhizosphere soil'),
                     guide=guide_legend(override.aes = list(alpha = 1,size=.symbolSize,linetype=0)))+
  scale_shape_manual(name="",values =shapeList2,
                     guide=guide_legend(override.aes = list(alpha = 1,size=.symbolSize)))

ggsave(filename="bacterial.hclustTree.pdf",bac.tree_plot,width=16,height=10,unit="cm")

ggsave(filename="bacterial.hclustTree_nolegend.pdf",
       bac.tree_plot+theme(legend.position = "none"),width=8,height=8,unit="cm")



##fungi

fgi.treedata = read.tree("fungi_BC_tree_rotated.nwk")
typeByNameGroup = split(row.names(sampleMeta),sampleMeta$typeByName)
fgi.treedata = groupOTU(fgi.treedata,typeByNameGroup)

fgi.tree_plot = ggtree(fgi.treedata,layout = "circular",linewidth=.linewidth,ladderize = F,
                       aes(color=group)) %<+% sampleMeta +
  geom_tippoint(aes(shape=periodSexByName,color=typeByName),
                fill="white",size=.symbolSize)+
  scale_color_manual(name="",values=c(colorList),
                     labels=c('Leaf','Root','Rhizosphere soil'),
                     guide=guide_legend(override.aes = list(alpha = 1,size=.symbolSize,linewidth=.linewidth,linetype=0)))+
  scale_shape_manual(name="",values =shapeList2,
                     guide=guide_legend(override.aes = list(alpha = 1,size=.symbolSize,linewidth=.linewidth,fill="white")))

ggsave(filename="fungi.hclustTree.pdf",fgi.tree_plot,width=18,height=18,unit="cm")
ggsave(filename="fungi.hclustTree_nolegend.pdf",
       fgi.tree_plot+theme(legend.position = "none"),width=8,height=8,unit="cm")

