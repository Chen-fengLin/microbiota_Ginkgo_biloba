cut -f 2,6 03classified_table_10clean.tsv |grep -v 'NA' |grep -v 'Incertae sedis' > 04OTU2order
cut -f2 04OTU2order |sed 1d|sort|uniq >05uniqOrder
perl TreeToConstraints.pl 08exactedOrderTree.nwk > 09constraint_alignment
R
	constraint = as.data.frame(matrix(read.table("09constraint_alignment")$V1,ncol=2,byrow=T))
	colnames(constraint)<-c("name","value")
	name2order = read.csv("06used_list.csv")[,c(1,8)]
	rownames(name2order) = paste0(">",name2order$old_taxonID_linked_genome_sequence)
	constraint$orderName = name2order[constraint$name,]$order
	order2OTU=read.table("04OTU2order",header=T)
	rownames(constraint) = constraint$orderName
	order2OTU$value = constraint[ order2OTU$order, ]$value
	clean.OTU.constraint = order2OTU[ !is.na(order2OTU$value), c(1,3) ]
	clean.OTU.constraint$OTU_id = paste0('>',clean.OTU.constraint$OTU_id)
	as.vector(t(as.matrix(clean.OTU.constraint)))->res
	write.table(res,"10my_OTU_constraint",sep="\t",quote=F,col.names=F,row.names=F)

mafft --auto 02OTU.filtered.fasta > 11OTU.filtered.mafft
fasttree -nt -gtr -constraints 10my_OTU_constraint < 11OTU.filtered.mafft > fgi.tree
