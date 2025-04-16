source("../../public/sourcetracker-master/src/SourceTracker.r")
mysource <-function(feature,env,source,sink,alpha1,alpha2,...){
  source.index = which(source)
  sink.index = which(sink)
  if((missing(alpha1) || missing(alpha2))){  ##未指定alpha1和alpha2
    if(length(unique(env[source.index]))>1){ ##且可以计算最佳alpha1和alpha2
      message("alpha is missing, calculating best alpha...")
      tune.res = tune.st(feature[source.index,], env[source.index],
                         alpha1=10**(-3:0), alpha2=10**(-3:0),)
      alpha1 = tune.res$best.alpha1
      alpha2 = tune.res$best.alpha2
    }else{
      alpha1=1e-3
      alpha2=1e-1
    }
  }
  st <- sourcetracker(feature[source.index,], env[source.index])
  results <- predict.sourcetracker(st,feature[sink.index,],alpha1 = alpha1,alpha2 = alpha2,...)
  colnames(results$proportions_sd) <- paste0(colnames(results$proportions_sd),"_sd")
  cbind(results$proportions,results$proportions_sd)
}




