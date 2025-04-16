library(ggplot2)
source("../../public/public.R")
typeName = data.frame(name=c('Leaf','Root','Rhizo.'),row.names = c('L','R','S'))
bac.lefse = read.table("lefse.result.bac.tsv",sep="\t")
bac.lefse$type = factor(typeName[bac.lefse$type,'name'],levels = c('Leaf','Root','Rhizo.'),ordered = T)
bac.lefse = bac.lefse[order(bac.lefse$Taxa),]
bac.lefse$TaxaShort = gsub(".*\\|","",bac.lefse$Taxa)
bac.lefse$Group = factor(ifelse(bac.lefse$Group=='F','Female','Male'),levels=c('Female','Male'),ordered = T)



#筛选出目和科的级别

bac.lefse.of = bac.lefse[grep("[o]__",bac.lefse$TaxaShort),]

bac.lefse.of.clean = bac.lefse.of[c('Group',"LDA",'type',"TaxaShort")]

count.LDA = as.data.frame(tapply(bac.lefse.of.clean$LDA, bac.lefse.of.clean$TaxaShort, length))
mean.LDA = as.data.frame(tapply(bac.lefse.of.clean$LDA, bac.lefse.of.clean$TaxaShort, mean))
identical(rownames(count.LDA),rownames(mean.LDA))
im.LDA = data.frame(TaxaShort=rownames(mean.LDA),mean.LDA=mean.LDA[,1],count.LDA=count.LDA[,1])
im.LDA.order = im.LDA[order(im.LDA$mean.LDA,im.LDA$TaxaShort,decreasing = T),]

important.Taxa = im.LDA.order$TaxaShort[1:10]

important.lefse = bac.lefse[bac.lefse$TaxaShort %in% important.Taxa,c("Group","LDA","Significance","type","TaxaShort")]
for (taxa in important.Taxa) {
  for(type in unique(important.lefse$type)){
    if(sum(important.lefse$TaxaShort==taxa & important.lefse$type==type)==0){
      important.lefse = rbind(important.lefse,
                              data.frame(Group=NA,LDA=0,Significance='ns',type=type,TaxaShort=taxa))
    }
  }
}

important.lefse$TaxaShort = factor(important.lefse$TaxaShort,levels = rev(important.Taxa),ordered = T)
important.lefse = important.lefse[order(important.lefse$TaxaShort),]

ggplot(data = important.lefse)+
  geom_tile(aes(x=type,y=TaxaShort,fill = LDA),color="grey",linewidth = .linewidth)+
  theme_bw()+genelateBWtheme()+
  scale_x_discrete(expand = c(0,0),labels=c("L","R","S"))+
  scale_y_discrete(expand = c(0,0),labels=function(s){sub('.__','',s)})+
  scale_fill_gradientn(colors=c("#fff8f8","#ff0000"),name="LDA",na.value = "#ffffff",
                       values = scales::rescale(c(2,5)),
                       limits=c(2,5),breaks=c(2,3,4,5),)+
  coord_fixed(ratio = 1)+
  theme(axis.title = element_blank(),axis.ticks = element_blank(),
        legend.key.height = unit(.symbolSizePT*3,'pt'),legend.key.width = unit(.symbolSizePT*3,'pt'),
        axis.ticks.length = unit(0,'pt'),panel.border = element_blank(),
        legend.frame = element_rect(colour = "#ddd",linewidth = .linewidth),
        legend.ticks = element_line(colour = "#ddd",linewidth = .linewidth),
        legend.title = element_text(size = .fontsizeTitlePT))
ggsave("LEfse.LDA.heat.bacOrder10.pdf",width=60,height =45,units = "mm")



#筛选出目和科的级别

bac.lefse.of = bac.lefse[grep("[f]__",bac.lefse$TaxaShort),]

bac.lefse.of.clean = bac.lefse.of[c('Group',"LDA",'type',"TaxaShort")]

count.LDA = as.data.frame(tapply(bac.lefse.of.clean$LDA, bac.lefse.of.clean$TaxaShort, length))
mean.LDA = as.data.frame(tapply(bac.lefse.of.clean$LDA, bac.lefse.of.clean$TaxaShort, mean))
identical(rownames(count.LDA),rownames(mean.LDA))
im.LDA = data.frame(TaxaShort=rownames(mean.LDA),mean.LDA=mean.LDA[,1],count.LDA=count.LDA[,1])
im.LDA.order = im.LDA[order(im.LDA$mean.LDA,im.LDA$TaxaShort,decreasing = T),]

important.Taxa = im.LDA.order$TaxaShort[1:10]

important.lefse = bac.lefse[bac.lefse$TaxaShort %in% important.Taxa,c("Group","LDA","Significance","type","TaxaShort")]
for (taxa in important.Taxa) {
  for(type in unique(important.lefse$type)){
    if(sum(important.lefse$TaxaShort==taxa & important.lefse$type==type)==0){
      important.lefse = rbind(important.lefse,
                              data.frame(Group=NA,LDA=0,Significance='ns',type=type,TaxaShort=taxa))
    }
  }
}

important.lefse$TaxaShort = factor(important.lefse$TaxaShort,levels = rev(important.Taxa),ordered = T)
important.lefse = important.lefse[order(important.lefse$TaxaShort),]

ggplot(data = important.lefse)+
  geom_tile(aes(x=type,y=TaxaShort,fill = LDA),color="grey",linewidth = .linewidth)+
  theme_bw()+genelateBWtheme()+
  scale_x_discrete(expand = c(0,0),labels=c("L","R","S"))+
  scale_y_discrete(expand = c(0,0),labels=function(s){sub('.__','',s)})+
  scale_fill_gradientn(colors=c("#fff8f8","#ff0000"),name="LDA",na.value = "#ffffff",
                       values = scales::rescale(c(2,5)),
                       limits=c(2,5),breaks=c(2,3,4,5),)+
  coord_fixed(ratio = 1)+
  theme(axis.title = element_blank(),axis.ticks = element_blank(),
        legend.key.height = unit(.symbolSizePT*3,'pt'),legend.key.width = unit(.symbolSizePT*3,'pt'),
        axis.ticks.length = unit(0,'pt'),panel.border = element_blank(),
        legend.frame = element_rect(colour = "#ddd",linewidth = .linewidth),
        legend.ticks = element_line(colour = "#ddd",linewidth = .linewidth),
        legend.title = element_text(size = .fontsizeTitlePT))
ggsave("LEfse.LDA.heat.bacFamily10.pdf",width=60,height =45,units = "mm")
