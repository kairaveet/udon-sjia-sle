#!/bin/bash

INPUTFILE=$1
SAMPLE=$(basename $INPUTFILE .txt) #removes the .txt from the file name only (parent folder name)
DIR=$(pwd)
#EXPDIR=$DIR/exp.MergedFiles-filtered-fold_increments_clean_UID.txt
#OUTPUT=$DIR/$INPUTFILE

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 24:00
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -M 128000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE


module load python

python /data/salomonis2/software/AltAnalyze/import_scripts/mergeFiles.py --i /data/salomonis2/LabFiles/Kairavee/Emely/H5files_Schulert/H5Files

EOF
