~/main/software/lotus2-master/lotus2 -i ~/main/bacteria_data/02unzipdata/ \
	-o ./lotusOut \
	-map mapping.txt \
	-t 16  \
	-verbosity 3 \
	-CL uparse\
	-taxAligner 0 \
	-amplicon_type SSU \
	-tax_group bacteria \
	-rdp_thr 0.5 \
	-buildPhylo 0 

