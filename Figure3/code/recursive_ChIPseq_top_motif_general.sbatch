#!/bin/bash

#SBATCH --time=100:05:00

#SBATCH --nodes=1 --ntasks-per-node=1

#SBATCH --mem=122280

#SBATCH -p pegasus_q

#SBATCH -A other_research


# Add any modules you might require:
module add R/3.5.1-foss-2018b


# Below here enter the commands to start your job
## check if any input parameter is empty 
if [[ -z "$idir" ]] || [[ -z "$Hpath" ]] || [[ -z "$refpath" ]] || [[ -z "$seqextractpath" ]] || [[ -z "$Rpath" ]] || [[ -z "$topmotif" ]]; then 
	echo -e "Usage: `basename $0` [-h/--help] <-idir/--input-dir-to-DMS string> <-Hpath/--HOMER-motif-path string> <-refpath/--reference-genome-fasta-path string> <-seqextractpath/--sequence-extractor-path string> <-Rpath/--Rscript-path string> <-motif/--top-motif string> <args>\n";  
	exit 1; 
fi

# cd to HOMER_custom
fpath=$Hpath
cd $fpath/data/knownTFs/motifs/

## add HOMER path
PATH=$PATH:$Hpath/bin/

## set directory
directory=$idir

# read target bed file and convert to fasta
perl $seqextractpath -r $refpath $directory/$target

# read background bed file and convert to fasta
perl $seqextractpath -r $refpath $directory/$background

## create known motif library
awk "BEGIN{IGNORECASE=1} FILENAME~ /^${topmotif}/" *.motif | awk '{print $FILENAME}' > $directory/${topmotif}.search.motif

## find locations of the motif of interest
findMotifs.pl $directory/${target%*.txt}.withSquence.fa  mouse $directory/ -mset vertebrates  -fastaBg  $directory/${background%*.txt}.withSquence.fa  -find $directory/${topmotif}.search.motif > $directory/${topmotif}.motif.txt

# make directory
mkdir $directory/${topmotif}	

## shift locations
Rscript $Rpath/Recursive_motif_search_dplyr.R  $directory ${topmotif} $directory/${topmotif}.motif.txt $directory/$target

## create fasta sequence 
perl $seqextractpath -r $refpath  $directory/${topmotif}/${topmotif}.txt
## find secondary motifs by masking
#findMotifs.pl ${name_DMS[idms]}/${topmotif:0:3}/${topmotif:0:3}.withSquence.fa  mouse ${name_DMS[idms]}/${topmotif:0:3}/  -mset vertebrates  -fastaBg ${name_DMS[idms]}/${name_DMS[idms]##*/}_background.withSquence.fa  -nomotif -p 10 -maskMotif ${name_DMS[idms]}/${topmotif:0:3}.search.motif
findMotifs.pl $directory/${topmotif}/${topmotif}.withSquence.fa  mouse $directory/${topmotif}/  -mset vertebrates  -fastaBg $directory/*_background.withSquence.fa  -nomotif -p 10 
## delete unwanted motifs sed does not work
ed -s $directory/${topmotif}/knownResults.txt <<< $'g/^SeqB/d\nw' 
ed -s $directory/${topmotif}/knownResults.txt  <<< $'g/^Unkno/d\nw' 