#!/bin/bash

## get directory
directory=/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/ChIP_seq_analysis/MEF2A/

## count number of dms clusters
num_sc=$(find $directory -mindepth 1 -maxdepth 1 -name 'm*' -type d  | wc -l)

## get names of DMS folder using IFS. If this does not work  the two lines below
IFS=' ' read -r -a name_sc <<< $(find $directory -mindepth 1  -maxdepth 1 -type d -name 'm*')

## define topmotif
motif=mef2a
## do for each dms
for ((isc=0;isc<=num_sc-1;isc++)); do
	## get name
	outname=${name_sc[isc]##*/}
	## run TRM
	sbatch --export=idir=$directory/$outname/,Hpath=/home/bsharmi6/HOMER_custom/,refpath=/home/bsharmi6/mm10bowtie2/mm10.fa,seqextractpath=/home/bsharmi6/MEME/sequence_extractor.pl,Rpath=/home/bsharmi6/NA_TF_project/R_scripts/,topmotif=$motif,target=${outname}_dms.txt,background=${outname}_background.txt recursive_ChIPseq_top_motif_general.sbatch
done

