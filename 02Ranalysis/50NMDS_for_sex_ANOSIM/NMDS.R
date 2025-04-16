library(vegan)
library(ggplot2)
library(rlang)
source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)

getCleanData <- function(rawdata){
  cleandata = rawdata[which(rowSums(rawdata)>0),]
  return(as.data.frame(t(cleandata)))
}

calculateNMDS <- function(type,period,feature,distance="bray",dim=2,try=10,trymax=50,maxit=200,seed=2024,...){
  set.seed(seed)
  msd_data = metaMDS(feature,distance=distance,k=dim,try=try,trymax=trymax,maxit=maxit,...)
  points_data = as.data.frame(scores(msd_data,display = "sites",tidy=TRUE))
  stress = msd_data$stress
  result = list(msd_data=msd_data,points = points_data, stress= stress)
  return(result)
}

plotNMDS <- function(type,period,NMDSdata,aesList,show.stressplot=FALSE,show.goodnessplot=FALSE,...,
                     ellipse.level=0.95,ellipseAES=aes(),ellipse.params=list()){
  # if(!is.null(ellipse.params$linewidth)){
  #   LINEPT = ggplot2::.pt*72.27/96
  #   ellipse.params$linewidth = ellipse.params$linewidth / LINEPT
  # }
  color.dir = list(S=soilColor,L=leafColor,R=rootColor)
  color=color.dir[[type]]
  shape=list(S1=shapeList[1],S2=shapeList[2],S3=shapeList[3])[[period]]
  stress = NMDSdata$stress
  points_data = as.data.frame(scores(NMDSdata,display = "sites",tidy=TRUE))
  if(show.stressplot){
    stressplot(NMDSdata)
  }
  if(show.goodnessplot){
    gof <- goodness(NMDSdata)
    plot(points_data$NMDS1,points_data$NMDS2,cex = gof*100,xlab="NMDS1",ylab="NMDS2")
  }
  meta.clean = sampleMeta[points_data$label,]
  plotdata = cbind(points_data,meta.clean)
  plot1 <- ggplot(data = plotdata,mapping = aes(x=NMDS1,y=NMDS2))+
    geom_point(mapping = aes(fill=sexByName),size=.symbolSize,color=color,shape=shape,show.legend = FALSE)+
    scale_fill_manual(values = c("transparent",color))+
    labs(x="NMDS1",y="NMDS2")+theme_bw()+
    theme(panel.grid = element_blank())
  
  plot2 <- plot1+
    layer(mapping = aes(linetype=sexByName), stat = StatEllipse, 
          geom = "path", position = "identity", show.legend = FALSE, 
          inherit.aes = TRUE, 
          params = list2(type = "t",level = ellipse.level, 
          segments = 51, na.rm = FALSE,linewidth=.linewidth,color=color))+
    scale_linetype_manual(values = c(2,1))
  result = list(p=plot1,ellipse = plot2)
  return(result)
  
}
formatANOSIM <- function(ANOSIM,Round){
  Pr=ANOSIM$signif
  if(Pr<0.001){
    resList = "sex: p<0.001"
  }else{
    resList = paste0("sex: p=",round(Pr,Round))
  }
  return(resList)
}
addtext <- function(plot0,stress=NA,ANOSIM=NA,Round = 3,position=c(-1,1)){
  textList = c()
  if(class(stress)!=class(NA)){
    textList = c(textList,paste0("stress: ",round(stress,Round)))
  }
  if(class(ANOSIM)!=class(NA)){
    ANOSIMList = formatANOSIM(ANOSIM,Round)
    textList = c(textList,ANOSIMList)
  }
  resplot = plot0+
    annotate(geom="text",x=position[1]*Inf,y=position[2]*Inf,
             label=paste(textList,collapse = "\n"),size = .fontsizeTitle/.pt,
             hjust=position[1]/2+0.5,vjust=position[2]/2+0.5)
  return(resplot)
}

calculateanosim<- function(y,groupingName,distance ="bray",permutations = 9999,seed=2025){
  meta.clean = factor(sampleMeta[rownames(y),groupingName])
  ## 这里默认进行贯序检验sequential test,变量顺序有意义，依次添加
  set.seed(seed)
  anosim(y,meta.clean,distance=distance,permutations=permutations)
}
