RowColApply <- function(X,rowIndex,rowFun,colIndex,colFun){
  X1 = NULL
  for (i in rownames(X)) {
    rowData = X[i,]
    rowData.group = as.data.frame(t(tapply(t(rowData),rowIndex,rowFun)),
                                  row.names = ifelse(is.numeric(i),paste0('V',i),i))
    X1 = rbind(X1,rowData.group)
  }
  
  X2=NULL
  for (j in colnames(X1)) {
    colData = X1[,j]
    colData.group = as.data.frame(tapply(colData,colIndex,colFun))
    colnames(colData.group) = ifelse(is.numeric(j),paste0('V',j),j)
    if(class(X2)==class(NULL))
      X2 = colData.group
    else
      X2 = cbind(X2,colData.group)
  }
  X2
}

star = function(pvalue){
  if(is.na(pvalue))return(NA)
  if(pvalue<0.001)return('***');
  if(pvalue<0.01)return('**');
  if(pvalue<0.05)return('*');
  return(NA)
}

checkSignif <- function(data,groupBy,group.in.rowname = FALSE,paired=FALSE){
  if(group.in.rowname)
    data = as.data.frame(t(data))
  groupNames = unique(groupBy)
  signifRes = NULL
  for(rowname in rownames(data)){
    for(name1 in groupNames){
      for(name2 in groupNames){
        if(name1<name2){
          vec1 = as.numeric(data[rowname,groupBy==name1])
          vec2 = as.numeric(data[rowname,groupBy==name2])
          mean1 = mean(vec1);mean2=mean(vec2);
          t_pvalue = t.test(vec1,vec2,paired = paired)$p.value
          w_pvalue = wilcox.test(vec1,vec2,paired = paired)$p.value
          signifRes = rbind(signifRes,data.frame(
            name = rowname,
            compare1 = name1,
            compare2 = name2,
            mean1 = mean1,mean2 = mean2,
            greater = ifelse(mean1>mean2,name1,ifelse(mean1<mean2,name2,NA)),
            t.pvalue = t_pvalue,t.star = star(t_pvalue),
            wilcox.pvalue = w_pvalue,wilcox.star = star(w_pvalue)
          ))
        }
      }
    }
  }
  signifRes
}
