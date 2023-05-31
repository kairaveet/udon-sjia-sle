INPUTFILE=$1
SAMPLE="CNA" #removes the .txt from the file name only (parent folder name)
DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 24:00
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -M 128000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE

cd $DIR
module load anaconda3
source activate /data/salomonis2/LabFiles/Kairavee/cna_env/



python3 /data/salomonis2/LabFiles/Kairavee/cna_evaluations_paper_revisions/cna_100cellsperpseudobulk.py


EOF
