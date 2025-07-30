# make 1 Mb regions
bedtools makewindows -g path_to_file/GRCh38.p13/STAR/chrNameLength.txt -w 1000000 | grep '^chr' > ../intermediate/GRCh38.p13_1Mb.bed
