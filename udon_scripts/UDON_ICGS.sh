#!/bin/bash

DIR=$(pwd)
INPUTFILE=$1
SAMPLE=$(basename $INPUTFILE .txt) #removes the .txt from the file name only (parent folder name)
EXPDIR=$DIR/MergedFiles_fold.txt
OUTPUT=$DIR/$INPUTFILE

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 24:00
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -M 64000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE

cd $DIR
module load R/3.2.2
module load python/2.7.5

python /data/salomonis2/software/AltAnalyze/AltAnalyze.py --runICGS yes --species Hs --platform RNASeq --excludeCellCycle no --removeOutliers no --restrictBy protein_coding --rho 0.3 --markerPearsonCutoff 0.2 --k 15 --FoldDiff 2 --SamplesDiffering 3 --column_method ward --row_method ward --expdir $EXPDIR --output $OUTPUT --expname "k15_GS"

EOF
