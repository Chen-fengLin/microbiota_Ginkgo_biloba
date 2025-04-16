source("../01enrich_analysis/enrich.R")


source("../../public/public.R")

enrichSexByType = function(diffRes,type,classifyTable,functionTable,threshold=0.05){
  diffRes.type = rbind(data.frame(diffRes[[type]]$S1,type=type,period='S1'),
                       data.frame(diffRes[[type]]$S2,type=type,period='S2'),
                       data.frame(diffRes[[type]]$S3,type=type,period='S3'))
  rownames(diffRes.type)=NULL
  diffRes.type.diff = diffRes.type[diffRes.type$p.value<=0.05,]
  all_OTU_list = unique(diffRes.type$OTU_id)
  diff_OTU_list = unique(diffRes.type.diff$OTU_id)
  
  diffRes.enrich = rbind(
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='phylum'),
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='class'),
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='order'),
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='family'),
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='genus'),
    calculateEnrich(diff_OTU_list,all_OTU_list,functionTable)
  )
  plots = plotEnrichBarSep(diffRes.enrich,threshold = threshold)
  list(diff = diffRes.type,enrich = diffRes.enrich,
       enrich.sig = diffRes.enrich[diffRes.enrich$p<=threshold,],
       plot = plots$plot,plot.type=plots$plot.type)
}

enrichSexByTypePeriod = function(diffRes,type,period,classifyTable,functionTable,threshold=0.05){
  diffRes.type = data.frame(diffRes[[type]][[period]],type=type,period=period)
  rownames(diffRes.type)=NULL
  diffRes.type.diff = diffRes.type[diffRes.type$p.value<=0.05,]
  all_OTU_list = unique(diffRes.type$OTU_id)
  diff_OTU_list = unique(diffRes.type.diff$OTU_id)
  
  diffRes.enrich = rbind(
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='phylum'),
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='class'),
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='order'),
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='family'),
    calculateTaxonEnrich(diff_OTU_list,all_OTU_list,classifyTable,levl='genus'),
    calculateEnrich(diff_OTU_list,all_OTU_list,functionTable)
  )
  diffRes.enrich = diffRes.enrich[diffRes.enrich$hub_in_Term>1,]
  plots = plotEnrichBarSep(diffRes.enrich,threshold = threshold)
  print(paste(type,period,sum(diffRes.enrich$p<=threshold)))
  list(diff = diffRes.type,enrich = diffRes.enrich,
       enrich.sig = diffRes.enrich[diffRes.enrich$p<=threshold,],
       plot = plots$plot,plot.type=plots$plot.type)
}


bacterial_function_table = read.table("../../00data/bacterial/function_table_10clean.tsv")
bacterial_classify_table = read.table("../../00data/bacterial/classified_table_10clean.tsv",sep="\t")
loadS = load("../52Diff_abund_sex/bac.diffSex.sex.Rdata")
bac.DiffRes = eval(as.symbol(loadS))
# bac
bac.leaf.res = enrichSexByType(bac.DiffRes,'L',bacterial_classify_table,bacterial_function_table)
bac.leaf.res$plot
bac.root.res = enrichSexByType(bac.DiffRes,'R',bacterial_classify_table,bacterial_function_table)
bac.root.res$plot
bac.soil.res = enrichSexByType(bac.DiffRes,'S',bacterial_classify_table,bacterial_function_table)
bac.soil.res$plot
# ggsave(filename = "bac-L-S123-sexDiff.Enrich.plot.pdf",width=10,height=12,units="cm")
p = function(a){a+theme(legend.key.width = unit(.symbolSizePT*3,'pt'),legend.key.height = unit(.symbolSizePT*3,'pt'))}
ggarrange(bac.leaf.res$plot.type,p(bac.leaf.res$plot),
          bac.root.res$plot.type,p(bac.root.res$plot),
          bac.soil.res$plot.type,p(bac.soil.res$plot),
          ncol = 2,nrow = 3,widths = c(3,3),align = "hv",heights = c(262,172,132))

ggsave(filename = "bac-LSR-sexDiff.Enrich.plot.pdf",width=11,height=20,units="cm")

bac.diffData = rbind(bac.leaf.res$diff,bac.root.res$diff,bac.soil.res$diff)
bac.EnrichData = rbind(bac.leaf.res$enrich,bac.root.res$enrich,bac.soil.res$enrich)
write.table(bac.diffData,"bac.diffData.LRS.tsv",sep="\t",quote = F)
write.table(bac.EnrichData,"bac.enrichData.LRS.tsv",sep="\t",quote = F)


bac.leafS1.res = enrichSexByTypePeriod(bac.DiffRes,'L','S1',bacterial_classify_table,bacterial_function_table)
bac.leafS2.res = enrichSexByTypePeriod(bac.DiffRes,'L','S2',bacterial_classify_table,bacterial_function_table)
bac.leafS3.res = enrichSexByTypePeriod(bac.DiffRes,'L','S3',bacterial_classify_table,bacterial_function_table)
bac.rootS1.res = enrichSexByTypePeriod(bac.DiffRes,'R','S1',bacterial_classify_table,bacterial_function_table)
bac.rootS2.res = enrichSexByTypePeriod(bac.DiffRes,'R','S2',bacterial_classify_table,bacterial_function_table)
bac.rootS3.res = enrichSexByTypePeriod(bac.DiffRes,'R','S3',bacterial_classify_table,bacterial_function_table)
bac.soilS1.res = enrichSexByTypePeriod(bac.DiffRes,'S','S1',bacterial_classify_table,bacterial_function_table)
bac.soilS2.res = enrichSexByTypePeriod(bac.DiffRes,'S','S2',bacterial_classify_table,bacterial_function_table)
bac.soilS3.res = enrichSexByTypePeriod(bac.DiffRes,'S','S3',bacterial_classify_table,bacterial_function_table)

bac1 = ggarrange(bac.leafS1.res$plot.type,p(bac.leafS1.res$plot)+ggtitle("Leaf S1")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          bac.leafS2.res$plot.type,p(bac.leafS2.res$plot)+ggtitle("Leaf S2")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          bac.leafS3.res$plot.type,p(bac.leafS3.res$plot)+ggtitle("Leaf S3")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          nrow = 3,ncol=2,align = "hv",heights = c(349,100,117),common.legend = T)


bac2 = ggarrange(bac.rootS1.res$plot.type,p(bac.rootS1.res$plot)+ggtitle("Root S1")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          bac.rootS2.res$plot.type,p(bac.rootS2.res$plot)+ggtitle("Root S2")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          bac.rootS3.res$plot.type,p(bac.rootS3.res$plot)+ggtitle("Root S3")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          bac.rootS3.res$plot.type,p(bac.rootS3.res$plot)+ggtitle("Root S3")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          nrow = 4,ncol=2,align = "hv",heights = c(64,155,173,173),common.legend = T)

bac3 = ggarrange(bac.soilS1.res$plot.type,p(bac.soilS1.res$plot)+ggtitle("Soil S1")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          bac.soilS2.res$plot.type,p(bac.soilS2.res$plot)+ggtitle("Soil S2")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          bac.soilS3.res$plot.type,p(bac.soilS3.res$plot)+ggtitle("Soil S3")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
          nrow = 3,ncol=2,align = "hv",heights = c(173,315,78),common.legend = T)

ggarrange(bac1,bac2,bac3,nrow = 1,ncol = 3,align = "hv")
ggsave(filename = "bac-LSR123-sexDiff.Enrich.plot.pdf",width=40,height=20,units="cm")

bac.diffData2 = rbind(bac.leafS1.res$diff,bac.leafS2.res$diff,bac.leafS3.res$diff,bac.rootS1.res$diff,bac.rootS2.res$diff,bac.rootS3.res$diff,bac.soilS1.res$diff,bac.soilS2.res$diff,bac.soilS3.res$diff)
bac.EnrichData2 = rbind(bac.leafS1.res$enrich,bac.leafS2.res$enrich,bac.leafS3.res$enrich,bac.rootS1.res$enrich,bac.rootS2.res$enrich,bac.rootS3.res$enrich,bac.soilS1.res$enrich,bac.soilS2.res$enrich,bac.soilS3.res$enrich)
write.table(bac.diffData2,"bac.diffData.LRS123.tsv",sep="\t",quote = F)
write.table(bac.EnrichData2,"bac.enrichData.LRS123.tsv",sep="\t",quote = F)




fungi_function_table = read.table("../../00data/fungi/function_table_10clean.tsv")
fungi_classify_table = read.table("../../00data/fungi/classified_table_10clean.tsv",sep="\t")
loadS = load("../52Diff_abund_sex/fgi.diffSex.sex.Rdata")
fgi.DiffRes = eval(as.symbol(loadS))
# fgi
fgi.leaf.res = enrichSexByType(fgi.DiffRes,'L',fungi_classify_table,fungi_function_table)
fgi.leaf.res$plot
fgi.root.res = enrichSexByType(fgi.DiffRes,'R',fungi_classify_table,fungi_function_table)
fgi.root.res$plot
fgi.soil.res = enrichSexByType(fgi.DiffRes,'S',fungi_classify_table,fungi_function_table)
fgi.soil.res$plot
# ggsave(filename = "bac-L-S123-sexDiff.Enrich.plot.pdf",width=10,height=12,units="cm")
p = function(a){a+theme(legend.key.width = unit(.symbolSizePT*3,'pt'),legend.key.height = unit(.symbolSizePT*3,'pt'))}
ggarrange(fgi.leaf.res$plot.type,p(fgi.leaf.res$plot),
          fgi.root.res$plot.type,p(fgi.root.res$plot),
          fgi.soil.res$plot.type,p(fgi.soil.res$plot),
          ncol = 2,nrow = 3,widths = c(3,3),align = "hv",heights = c(112,122,82),common.legend = T,legend = "right")

ggsave(filename = "fgi-LSR-sexDiff.Enrich.plot.pdf",width=13,height=11.17,units="cm")

fgi.diffData = rbind(fgi.leaf.res$diff,fgi.root.res$diff,fgi.soil.res$diff)
fgi.EnrichData = rbind(fgi.leaf.res$enrich,fgi.root.res$enrich,fgi.soil.res$enrich)
write.table(fgi.diffData,"fgi.diffData.LRS.tsv",sep="\t",quote = F)
write.table(fgi.EnrichData,"fgi.enrichData.LRS.tsv",sep="\t",quote = F)


fgi.leafS1.res = enrichSexByTypePeriod(fgi.DiffRes,'L','S1',fungi_classify_table,fungi_function_table)
fgi.leafS2.res = enrichSexByTypePeriod(fgi.DiffRes,'L','S2',fungi_classify_table,fungi_function_table)
fgi.leafS3.res = enrichSexByTypePeriod(fgi.DiffRes,'L','S3',fungi_classify_table,fungi_function_table)
fgi.rootS1.res = enrichSexByTypePeriod(fgi.DiffRes,'R','S1',fungi_classify_table,fungi_function_table)
fgi.rootS2.res = enrichSexByTypePeriod(fgi.DiffRes,'R','S2',fungi_classify_table,fungi_function_table)
fgi.rootS3.res = enrichSexByTypePeriod(fgi.DiffRes,'R','S3',fungi_classify_table,fungi_function_table)
fgi.soilS1.res = enrichSexByTypePeriod(fgi.DiffRes,'S','S1',fungi_classify_table,fungi_function_table)
fgi.soilS2.res = enrichSexByTypePeriod(fgi.DiffRes,'S','S2',fungi_classify_table,fungi_function_table)
fgi.soilS3.res = enrichSexByTypePeriod(fgi.DiffRes,'S','S3',fungi_classify_table,fungi_function_table)

p2 = function(a){
  
  a+scale_fill_gradientn(colors=c("#fff","#fff","#fff8f8","#f00"),name="p value",
                         values = scales::rescale(c(0,-log10(0.05)-1e-10,-log10(0.05),4)),
                         limits=c(0,4),breaks=c(0,-log10(0.05),4),labels=c(1,0.05,paste0("1e-",4)))+
    theme(legend.key.width = unit(.symbolSizePT*3,'pt'),legend.key.height = unit(.symbolSizePT*3,'pt'))
}

fgi1 = ggarrange(fgi.leafS1.res$plot.type,p2(fgi.leafS1.res$plot)+ggtitle("Leaf S1")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
                 fgi.leafS2.res$plot.type,p2(fgi.leafS2.res$plot)+ggtitle("Leaf S2")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT))+scale_x_continuous(expand = c(0,0,0.05,0),name = "Number of OTUs",limits = c(0,5)),
                 fgi.leafS3.res$plot.type,p2(fgi.leafS3.res$plot)+ggtitle("Leaf S3")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
                 nrow = 3,ncol=2,align = "hv",heights = c(144,125,125),common.legend = T)


fgi2 = ggarrange(fgi.rootS1.res$plot.type,p2(fgi.rootS1.res$plot)+ggtitle("Root S1")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
                 fgi.rootS2.res$plot.type,p2(fgi.rootS2.res$plot)+ggtitle("Root S2")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT))+scale_x_continuous(expand = c(0,0,0.05,0),name = "Number of OTUs",limits = c(0,5)),
                 fgi.rootS3.res$plot.type,p2(fgi.rootS3.res$plot)+ggtitle("Root S3")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
                 nrow = 3,ncol=2,align = "hv",heights = c(144,125,125),common.legend = T)

fgi3 = ggarrange(fgi.soilS1.res$plot.type,p2(fgi.soilS1.res$plot)+ggtitle("Soil S1")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT))+scale_x_continuous(expand = c(0,0,0.05,0),name = "Number of OTUs",limits = c(0,10)),
                 fgi.soilS2.res$plot.type,p2(fgi.soilS2.res$plot)+ggtitle("Soil S2")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT)),
                 fgi.soilS3.res$plot.type,p2(fgi.soilS3.res$plot)+ggtitle("Soil S3")+theme(title = element_text(face = "bold",size = .fontsizeTitlePT))+scale_x_continuous(expand = c(0,0,0.05,0),name = "Number of OTUs",limits = c(0,5)),
                 nrow = 3,ncol=2,align = "hv",heights = c(160,131,103),common.legend = T)

ggarrange(fgi1,fgi2,fgi3,nrow = 1,ncol = 3,align = "hv",common.legend = T)
ggsave(filename = "fgi-LSR123-sexDiff.Enrich.plot.pdf",width=40,height=14,units="cm")

fgi.diffData2 = rbind(fgi.leafS1.res$diff,fgi.leafS2.res$diff,fgi.leafS3.res$diff,fgi.rootS1.res$diff,fgi.rootS2.res$diff,fgi.rootS3.res$diff,fgi.soilS1.res$diff,fgi.soilS2.res$diff,fgi.soilS3.res$diff)
fgi.EnrichData2 = rbind(fgi.leafS1.res$enrich,fgi.leafS2.res$enrich,fgi.leafS3.res$enrich,fgi.rootS1.res$enrich,fgi.rootS2.res$enrich,fgi.rootS3.res$enrich,fgi.soilS1.res$enrich,fgi.soilS2.res$enrich,fgi.soilS3.res$enrich)
write.table(fgi.diffData2,"fgi.diffData.LRS123.tsv",sep="\t",quote = F)
write.table(fgi.EnrichData2,"fgi.enrichData.LRS123.tsv",sep="\t",quote = F)
