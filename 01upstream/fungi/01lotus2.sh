~/main/software/lotus2-master/lotus2 -i ~/main/fungi_data/02unzipdata/ \
	-o ./lotusOut \
	-map mapping.txt \
	-t 16  \
	-verbosity 3 \
	-CL uparse\
	-taxAligner blast \
	-refDB UNITE \
	-amplicon_type ITS \
	-tax_group fungi \
	-rdp_thr 0.5 \
	-buildPhylo 0 

