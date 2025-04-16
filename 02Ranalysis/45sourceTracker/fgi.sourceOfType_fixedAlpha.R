source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",header = TRUE,sep="\t")
sampleMeta = fmtMeta(sampleMeta)

fgi.alldata = read.table("../../00data/fungi/feature_table_10clean.tsv")
fgi.sampleMeta = sampleMeta[colnames(fgi.alldata),]
fgi.alldata = t(as.matrix(fgi.alldata))

Envs = paste0(fgi.sampleMeta$type,fgi.sampleMeta$period)
TypePeriod = paste0(fgi.sampleMeta$type,fgi.sampleMeta$period)

source("mysource.R")

results = list()

results[['LS1']]=mysource(fgi.alldata,Envs,sink = TypePeriod=='LS1',alpha1 = 0.001,alpha2 = 0.1,
                        source = TypePeriod=='RS1' | TypePeriod=='SS1')


results[['RS1']]=mysource(fgi.alldata,Envs,sink = TypePeriod=='RS1',alpha1 = 0.001,alpha2 = 0.1,
                        source = TypePeriod=='SS1')


results[['LS2']]=mysource(fgi.alldata,Envs,sink = TypePeriod=='LS2',alpha1 = 0.001,alpha2 = 0.1,
                        source = TypePeriod=='RS2' | TypePeriod=='SS2' | TypePeriod=='LS1')


results[['RS2']]=mysource(fgi.alldata,Envs,sink = TypePeriod=='RS2',alpha1 = 0.001,alpha2 = 0.1,
                        source = TypePeriod=='SS2' | TypePeriod=='RS1')

results[['SS2']]=mysource(fgi.alldata,Envs,sink = TypePeriod=='SS2',alpha1 = 0.001,alpha2 = 0.1,
                          source = TypePeriod=='SS1')

results[['LS3']]=mysource(fgi.alldata,Envs,sink = TypePeriod=='LS3',alpha1 = 0.001,alpha2 = 0.1,
                        source = TypePeriod=='RS3' | TypePeriod=='SS3' | TypePeriod=='LS2')


results[['RS3']]=mysource(fgi.alldata,Envs,sink = TypePeriod=='RS3',alpha1 = 0.001,alpha2 = 0.1,
                        source = TypePeriod=='SS3' | TypePeriod=='RS2')

results[['SS3']]=mysource(fgi.alldata,Envs,sink = TypePeriod=='SS3',alpha1 = 0.001,alpha2 = 0.1,
                          source = TypePeriod=='SS2')

save(results,file = "fgi.sourceTrack.byType.fixedAlpha.Rdata")