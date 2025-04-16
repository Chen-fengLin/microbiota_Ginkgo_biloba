source("../../public/public.R")
feature_table = read.table(file="../../00data/fungi/feature_table_10clean.tsv")
sampleMeta = read.table("../../00data/sample_meta.tsv",sep="\t",header = T)
sampleMeta = fmtMeta(sampleMeta)
sampleMeta = sampleMeta[colnames(feature_table),]

classified_table = read.table('../../00data/fungi/classified_table_10clean.tsv',header=TRUE, sep="\t",na.strings = "NA",comment.char = "")
tax_table = data.frame(
  Kingdom = paste0('k__',ifelse(is.na(classified_table$domain),"",classified_table$domain)),
  Phylum = paste0('p__',ifelse(is.na(classified_table$phylum),"",classified_table$phylum)),
  Class = paste0('c__',ifelse(is.na(classified_table$class),"",classified_table$class)),
  Order = paste0('o__',ifelse(is.na(classified_table$order),"",classified_table$order)),
  Family = paste0('f__',ifelse(is.na(classified_table$family),"",classified_table$family)),
  Genus = paste0('g__',ifelse(is.na(classified_table$genus),"",classified_table$genus))
)
row.names(tax_table)=classified_table$OTU_id

library(microeco)

lefse.result.frame = NULL
lefse.abund.frame=NULL
for(type in c('L','R','S')){
  # for(stage in c('S1','S2','S3')){
    sampleMeta.cur = sampleMeta[sampleMeta$type==type,]
    feature_table.cur = feature_table[,sampleMeta$type==type]
    tax_table.cur = tax_table
    dataset.cur <- microtable$new(sample_table = sampleMeta.cur,
                              otu_table = feature_table.cur, 
                              tax_table = tax_table.cur)
    lefse.cur = trans_diff$new(dataset = dataset.cur, 
                               method = "lefse", 
                               group = "sex", 
                               filter_thres = 0.001,
                               alpha = 10, 
                               p_adjust_method = "none",
                               lefse_subgroup = NULL)
    lefse.result.cur = lefse.cur$res_diff
    lefse.result.cur$taxa = row.names(lefse.result.cur)
    lefse.result.cur$type=type
    row.names(lefse.result.cur)=NULL
    lefse.result.frame = rbind(lefse.result.frame,lefse.result.cur)
    
    lefse.abund.cur = lefse.cur$res_abund
    lefse.abund.cur$type=type
    row.names(lefse.abund.cur)=NULL
    lefse.abund.frame = rbind(lefse.abund.frame,lefse.abund.cur)
 
  # }
}

write.table(lefse.result.frame,file = "lefse.result.fgi.tsv",sep="\t",quote = F,col.names = T,row.names = T)
write.table(lefse.abund.frame,file = "lefse.abund.fgi.tsv",sep="\t",quote = F,col.names = T,row.names = T)
