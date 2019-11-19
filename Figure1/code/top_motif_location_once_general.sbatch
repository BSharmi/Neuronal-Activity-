#!/bin/bash

#SBATCH --time=100:05:00

#SBATCH --nodes=1 --ntasks-per-node=1

#SBATCH --mem=122280

#SBATCH -p pegasus_q

#SBATCH -A other_research

# Below here enter the commands to start your job
## check if any input parameter is empty 
if [[ -z "$idir" ]] || [[ -z "$Hpath" ]] || [[ -z "$motif" ]] ; then 
	echo -e "Usage: `basename $0` [-h/--help] <-idir/--input-dir-to-DMS string> <-Hpath/--HOMER-motif-path string> <-motif/--motif string> <args>\n";  
	exit 1; 
fi

# cd to HOMER_custom
fpath=$Hpath
cd $fpath/data/knownTFs/motifs/

## add HOMER path
PATH=$PATH:$Hpath/bin/

## get path to dms files
dmspath=/home/bsharmi6/NA_TF_project/scMethylome/all_DMS_known_TFs/

## create directory
mkdir $idir/$motif

## create known motif library
awk "BEGIN{IGNORECASE=1} FILENAME~ /^${motif}/" *.motif | awk '{print $FILENAME}' > $idir/$motif/${motif}.search.motif

# run find motifs to find all known motifs
findMotifs.pl $dmspath/*dms.withSquence.fa  mouse $idir/$motif/ -mset vertebrates  -fastaBg  $dmspath/*background.withSquence.fa  -find $idir/$motif/${motif}.search.motif > $idir/$motif/${motif}.motif.txt
