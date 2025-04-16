sourceMean = function(sourceTable){
  Mean=apply(sourceTable[,1:(ncol(sourceTable)/2)],2, mean)*100
  Sd=apply(sourceTable[,1:(ncol(sourceTable)/2)],2, sd)*100
  N=apply(sourceTable[,1:(ncol(sourceTable)/2)],2, length)
  Se = Sd/sqrt(N)
  as.data.frame(t(rbind(Mean,Sd,N,Se)))
}

sourceSummary = function(file){
  loadS = load(file)
  sourceRes = eval(as.symbol(loadS))
  lapply(sourceRes,sourceMean)
}


bac.fixedAlpha = sourceSummary("bac.sourceTrack.byType.fixedAlpha.Rdata")

fgi.fixedAlpha = sourceSummary("fgi.sourceTrack.byType.fixedAlpha.Rdata")
