library(vegan)
library(ggplot2)
library(rlang)
source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",header = T,sep="\t")
sampleMeta = fmtMeta(sampleMeta)

getCleanData <- function(rawdata){
  cleandata = rawdata[which(rowSums(rawdata)>0),]
  return(t(cleandata))
}

calculateNMDS <- function(feature,distance="bray",dim=2,try=10,trymax=50,maxit=200,seed=2024,...){
  set.seed(seed)
  msd_data = metaMDS(feature,distance=distance,k=dim,try=try,trymax=trymax,maxit=maxit,...)
  points_data = as.data.frame(scores(msd_data,display = "sites",tidy=TRUE))
  stress = msd_data$stress
  result = list(msd_data=msd_data,points = points_data, stress= stress)
  return(result)
}
plotNMDS <- function(NMDSdata,aesList,show.stressplot=FALSE,show.goodnessplot=FALSE,...,
                     ellipse.level=0.95,ellipseAES=aes(),ellipse.params=list()){
  if(!is.null(ellipse.params$linewidth)){
    LINEPT = ggplot2::.pt*72.27/96
    ellipse.params$linewidth = ellipse.params$linewidth / LINEPT
  }
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
    geom_point(mapping = aesList,...)+
    labs(x="NMDS1",y="NMDS2")+theme_bw()+
    theme(panel.grid = element_blank())
  
  plot2 <- plot1+
    layer(mapping = ellipseAES, stat = StatEllipse, 
          geom = "path", position = "identity", show.legend = NA, 
          inherit.aes = TRUE, 
          params = list2(type = "t",level = ellipse.level, 
          segments = 51, na.rm = FALSE,!!!ellipse.params))
  result = list(p=plot1,ellipse = plot2)
  return(result)
  
}
formatPREM <- function(PERMANOVA){
  signifList = PERMANOVA[grep("[*]",PERMANOVA$signif),]
  nameList = row.names(signifList)
  nameList = gsub(":"," × ",nameList)
  resList = paste0(nameList,": ",signifList$signif)
  return(resList)
}
addtext <- function(plot0,stress=NA,PERMANOVA=NA,fontsize=.fontsize,stressRound = 3,position=c(-1,1)){
  textList = c()
  if(class(stress)!=class(NA)){
    textList = c(textList,paste0("stress: ",round(stress,stressRound)))
  }
  if(class(PERMANOVA)!=class(NA)){
    PERMANOVAList = formatPREM(PERMANOVA)
    textList = c(textList,PERMANOVAList)
  }
  resplot = plot0+
    annotate(geom="text",x=position[1]*Inf,y=position[2]*Inf,
             label=paste(textList,collapse = "\n"),size = fontsize,
             hjust=position[1]/2+0.5,vjust=position[2]/2+0.5)
  return(resplot)
}

calculatePERMANOVA<- function(y,myformula,...,method="bray",permutations = 999,by = "terms",seed=2024){
  
  meta.clean = sampleMeta[rownames(y),]
  ## 这里默认进行贯序检验sequential test,变量顺序有意义，依次添加
  set.seed(seed)
  result0 = adonis2(myformula,meta.clean,method=method,permutations=permutations,by=by,...)
  result.frame = as.data.frame(result0)
  colnames(result.frame) = c('Df','SumOfSqs','R2','F','Pr')
  result.frame$signif = with(result.frame,ifelse(Pr<=0.001,"***",ifelse(Pr<=0.01,"**",
                             ifelse(Pr<0.05,"*",ifelse(Pr<0.1,"."," ")))))
  
  ##检验离散度
  result.frame$disper_P = NA
  y.dist <- vegdist(y,method = method,binary = F)
  for(item in row.names(result.frame)[1:(nrow(result.frame)-2)]){
    if(item %in% colnames(meta.clean)){
      item.dispersion <- betadisper(y.dist,group=meta.clean[[item]])
      item.disp_test = permutest(item.dispersion)
      item.pvalue = item.disp_test$tab$`Pr(>F)`[1]
      result.frame[item,'disper_P']=item.pvalue
    }
  }
  return(result.frame)
}