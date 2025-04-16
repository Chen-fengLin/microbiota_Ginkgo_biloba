mafft --auto OTU.filtered.fasta > OTU.filtered.mafft
nohup fasttree -nt -gtr OTU.filtered.mafft > bac.tree
