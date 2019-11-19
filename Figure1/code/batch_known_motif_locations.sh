#!/bin/bash

#################################### OLD DO NOT DELETE #################################
# ## get directory
# directory=/home/bsharmi6/NA_TF_project/scMethylome/all_DMS_known_TFs/
# outpath=/home/bsharmi6/NA_TF_project/scMethylome/Clustering_in_ETRM

# ## count number of dms clusters
# num_TFs=$(find $directory -mindepth 1 -maxdepth 1 -name '*minus' -type d  | wc -l)

# ## get names of DMS folder using IFS. If this does not work  the two lines below
# IFS=' ' read -r -a name_TFs <<< $(find $directory -mindepth 1  -maxdepth 1 -type d -name '*minus')

# ## find motif locations for each TF
# for ((iTF=0;iTF<=num_TFs-1;iTF++)); do
# 	## get name
# 	TF_minus=${name_TFs[iTF]##*/}  
# 	TF=${TF_minus%_*}
# 	## call job
# 	sbatch --export=idir=$outpath,Hpath=/home/bsharmi6/HOMER_custom/,motif=$TF top_motif_location_once_general.sbatch
# done

## rename X-b to X-box (bad temporary approach)
#mv $directory/Enrichment_motifs_in_neurons/TRM/X.trm.txt $directory/Enrichment_motifs_in_neurons/TRM/X.box.trm.txt 

#################################### NEW ##################################################
# ## get directory
directory=/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/Enrichment_motifs_in_neurons/TRM/
outpath=/home/bsharmi6/NA_TF_project/scMethylome/Clustering_in_ETRM/subset_DMS/

## count number of dms clusters
num_TFs=$(find $directory -mindepth 1 -maxdepth 1 -name '*'  | wc -l)

## get names of DMS folder using IFS. If this does not work  the two lines below
IFS=' ' read -r -a name_TFs <<< $(find $directory -mindepth 1  -maxdepth 1 -name '*')

## find motif locations for each TF
for ((iTF=0;iTF<=num_TFs-1;iTF++)); do
	## get name
	TF_minus=${name_TFs[iTF]##*/};
	TF_minus=${TF_minus%%.*}  
	TF=${TF_minus,,}
	## call job
	sbatch --export=idir=$outpath,Hpath=/home/bsharmi6/HOMER_custom/,motif=${TF:0:3} top_motif_location_once_general.sbatch
done