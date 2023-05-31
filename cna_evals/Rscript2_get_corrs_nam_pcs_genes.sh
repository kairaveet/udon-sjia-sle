INPUTFILE=$1
SAMPLE=$(basename $INPUTFILE .txt) #removes the .txt from the file name only (parent folder name)
DIR=$(pwd)


cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 50:00
#BSUB -n 2
#BSUB -R "span[ptile=4]"
#BSUB -M 250000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE

cd $DIR
module load R/4.0.0
module load python/2.7.5

Rscript get_corrs_nam_pcs_genes.R


EOF
#for i in a; do cellranger.sh $i | bsub; done
