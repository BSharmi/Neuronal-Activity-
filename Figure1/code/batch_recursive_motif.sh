#!/bin/bash

## get directory
directory=/home/bsharmi6/NA_TF_project/scMethylome/Clustering_in_ETRM/subset_DMS/

## count number of dms clusters
num_DMS=$(find $directory -mindepth 1 -maxdepth 1 -name '*' -type d  | wc -l)

## get names of DMS folder using IFS. If this does not work  the two lines below
IFS=' ' read -r -a name_DMS <<< $(find $directory -mindepth 1  -maxdepth 1 -type d -name '*')

## do for each dms
for ((idms=0;idms<=num_DMS-1;idms++)); do
	## get name
	outname=${name_DMS[idms]##*/}
	## run clustering
	sbatch --export=idir=$directory/$outname/,Hpath=/home/bsharmi6/HOMER_custom/,refpath=/home/bsharmi6/mm10bowtie2/mm10.fa,seqextractpath=/home/bsharmi6/MEME/sequence_extractor.pl,Rpath=/home/bsharmi6/NA_TF_project/R_scripts/ recursive_top_pval_mask_general.sbatch
done