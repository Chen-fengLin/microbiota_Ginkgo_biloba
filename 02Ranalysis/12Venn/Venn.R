source("../../public/public.R")
sampleMeta = read.table("../../00data/sample_meta.tsv",sep="\t",header = T)
sampleMeta = fmtMeta(sampleMeta)

bac.alldata = read.table("../../00data/bacterial/feature_table_20norm.tsv")
bac.sampleMeta = sampleMeta[colnames(bac.alldata),]

  
bac.leaf.alldata = bac.alldata[,bac.sampleMeta$type=='L']
bac.root.alldata = bac.alldata[,bac.sampleMeta$type=='R']
bac.soil.alldata = bac.alldata[,bac.sampleMeta$type=='S']

bac.leaf.ID = rownames(bac.leaf.alldata)[rowSums(bac.leaf.alldata>0)>1|rowSums(bac.leaf.alldata)/1e3>0.1]
bac.root.ID = rownames(bac.leaf.alldata)[rowSums(bac.root.alldata>0)>1|rowSums(bac.leaf.alldata)/1e3>0.1]
bac.soil.ID = rownames(bac.leaf.alldata)[rowSums(bac.soil.alldata>0)>1|rowSums(bac.leaf.alldata)/1e3>0.1]


length(bac.leaf.ID)
length(bac.root.ID)
length(bac.soil.ID)
sum(bac.leaf.ID %in% bac.root.ID)
sum(bac.root.ID %in% bac.soil.ID)
sum(bac.leaf.ID %in% bac.soil.ID)

sum(bac.leaf.ID %in% bac.soil.ID & bac.leaf.ID %in% bac.root.ID)


fgi.alldata = read.table("../../00data/fungi/feature_table_10clean.tsv")
fgi.sampleMeta = sampleMeta[colnames(fgi.alldata),]

fgi.leaf.alldata = fgi.alldata[,fgi.sampleMeta$type=='L']
fgi.root.alldata = fgi.alldata[,fgi.sampleMeta$type=='R']
fgi.soil.alldata = fgi.alldata[,fgi.sampleMeta$type=='S']

fgi.leaf.ID = rownames(fgi.leaf.alldata)[rowSums(fgi.leaf.alldata>0)>1|rowSums(fgi.leaf.alldata)/1e3>0.1]
fgi.root.ID = rownames(fgi.leaf.alldata)[rowSums(fgi.root.alldata>0)>1|rowSums(fgi.root.alldata)/1e3>0.1]
fgi.soil.ID = rownames(fgi.leaf.alldata)[rowSums(fgi.soil.alldata>0)>1|rowSums(fgi.soil.alldata)/1e3>0.1]


length(fgi.leaf.ID)
length(fgi.root.ID)
length(fgi.soil.ID)
sum(fgi.leaf.ID %in% fgi.root.ID)
sum(fgi.root.ID %in% fgi.soil.ID)
sum(fgi.leaf.ID %in% fgi.soil.ID)

sum(fgi.leaf.ID %in% fgi.soil.ID & fgi.leaf.ID %in% fgi.root.ID)
