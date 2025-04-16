myconcat = function(vec,OTU_list){
  paste(OTU_list[vec>0],sep="|",collapse = "|")
}

calculateEnrich = function(hub_list,all_list,function_table,type='function'){
  hub_function = function_table[hub_list,]
  all_function = function_table[all_list,]
  hub_count = length(hub_list)
  all_count = length(all_list)
  hub_function_count = colSums(hub_function)
  all_function_count = colSums(all_function)
  hub_function_count.filtered = hub_function_count[hub_function_count>0]
  all_function_count.filtered = all_function_count[hub_function_count>0]
  p_value = 1-phyper(hub_function_count.filtered-1,all_function_count.filtered,all_count-all_function_count.filtered,hub_count)
  rich_factor = (hub_function_count.filtered/hub_count)/(all_function_count.filtered/all_count)
  hub_function.filtered = as.data.frame(hub_function[hub_function_count>0])
  # print(hub_function.filtered)
  hub_function.OTUname = apply(hub_function.filtered,2,OTU_list = rownames(hub_function.filtered),myconcat)
  if(length(hub_function_count.filtered)==0){
    return(data.frame(Term=c(),
                      count_in_Term=c(),
                      hub_in_Term=c(),
                      rich_factor=c(),
                      OTUname=c(),
                      p=c(),
                      type=c()
    ))
  }
  result = data.frame(Term = names(hub_function_count.filtered),
                      count_in_Term = all_function_count.filtered,
                      hub_in_Term = hub_function_count.filtered,
                      rich_factor = rich_factor,
                      OTUname = hub_function.OTUname,
                      p=p_value,
                      type=type
                      )
  result[order(result$p),]
}

sort=function(vec){
  vec[order(vec)]
}

calculateTaxonEnrich = function(hub_list,all_list,classify_table,levl){
  row.names(classify_table) = classify_table$OTU_id
  classify_table.using = classify_table[levl]
  Taxon.list = sort(unique(classify_table.using[,1]))
  Taxon.list = Taxon.list[!is.na(Taxon.list)]
  Taxon.frame = as.data.frame(matrix(0,nrow = nrow(classify_table.using),ncol=length(Taxon.list)),
                              row.names = rownames(classify_table.using))
  colnames(Taxon.frame) = Taxon.list
  for(OTU_id in row.names(classify_table.using)){
    if(classify_table.using[OTU_id,1]!='' & !is.na(classify_table.using[OTU_id,1])){
      Taxon.frame[OTU_id,classify_table.using[OTU_id,1]]=1
    }
  }
  
  calculateEnrich(hub_list,all_list,Taxon.frame,type=levl)
}


library(ggplot2)
library(ggpubr)
plotEnrichBarSep = function(EnrichRes,threshold=0.05){
  plotData = EnrichRes[EnrichRes$p<=threshold,]
  plotData$termTag = paste(plotData$type,plotData$Term,sep = "-")
  plotData$type = factor(plotData$type,levels = c('domain','phylum','class','order','family','genus','function'),ordered = T)
  plotData$termTag = factor(plotData$termTag,ordered = TRUE,
                            levels = plotData$termTag[order(plotData$type,plotData$p,decreasing = T)])
  plotData = plotData[order(plotData$termTag),]
  plotData$logp = -log10(plotData$p)
  maxP = ceiling(max(plotData$logp)/2)*2
  # maxP = 8 
  plot.enrich = ggplot(data=plotData)+
    geom_bar(aes(x=hub_in_Term,y=termTag,fill=logp),color="#ddd",stat = "identity",linewidth=.linewidth,)+
    theme_bw()+genelateBWtheme()+
    scale_x_continuous(expand = c(0,0,0.05,0),name = "Number of OTUs")+
    scale_y_discrete(breaks = plotData$termTag,labels = plotData$Term,name = NULL,expand = c(0,0))+
    scale_fill_gradientn(colors=c("#fff","#fff","#fff8f8","#f00"),name="p value",
                         values = scales::rescale(c(0,-log10(0.05)-1e-10,-log10(0.05),maxP)),
                         limits=c(0,maxP),breaks=c(0,-log10(0.05),maxP),labels=c(1,0.05,paste0("1e-",maxP)))+
    # scale_color_manual(values = c("#333","#ddd","#333","#ddd","#333","#ddd"),
    #                    breaks = c('domain','phylum','class','order','family','genus','function'),
    #                    guide=FALSE)+
    theme(legend.frame = element_rect(colour = "#ddd",linewidth = .linewidth),
          legend.ticks = element_line(colour = "#ddd",linewidth = .linewidth),
          legend.title = element_text(colour = "#000",size = .fontsizeTitlePT),
          legend.box.background = element_rect(colour = "#000",linewidth = .linewidth/2,fill = 'transparent'),
          legend.position = "inside",legend.position.inside = c(1,0),legend.justification = c(1,0),
          axis.ticks.y = element_blank(),axis.ticks.length.y = unit(0,'pt'))
  plot.type=ggplot(data=plotData)+
    geom_bar(aes(x=1,y=termTag,fill = type),stat="identity",width = 1,show.legend = F,color="transparent")+
    theme_bw()+genelateBWtheme()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0),
                     labels=plotData$type[tapply(c(1:nrow(plotData)),as.character(plotData$type),median)],
                     breaks=plotData$termTag[tapply(c(1:nrow(plotData)),as.character(plotData$type),median)])+
    scale_fill_manual(values = c("#4ddb4d","#66E066","#7FE57F","#99EB99","#B3f0b3","#CCF5CC","#E5FAE5"),
                      breaks = c('domain','phylum','class','order','family','genus','function'))+
    theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),
          axis.text.y = element_text(size = .fontsizeTitlePT,face = "bold"),
          plot.background = element_blank(),panel.background = element_blank(),panel.border = element_blank())
  list(plot=plot.enrich,plot.type=plot.type)
}
plotEnrichBar = function(EnrichRes,threshold=0.05,plot.par=c(1,6)){
  plot = plotEnrichBarSep(EnrichRes,threshold=threshold)
  ggarrange(plot$plot.type,plot$plot,ncol = 2,nrow=1,align = "h",widths = plot.par)
}
