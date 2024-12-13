#!/bin/bash

#SBATCH -N 1
#SBATCH -t 4:30:00
#SBATCH -n 1
#SBATCH -D '/tigress/SHAEVITZ/dsmendez/diana'
#SBATCH --mem=32000
#SBATCH --output='logs/out-%A_%a.log'
#SBATCH --error='logs/error-%A_%a.err'
#SBATCH --mail-type=end
#SBATCH --mail-user=dsmendez@princeton.edu
#SBATCH --array=1-55
#-55
#-98 +15
#-110
#3-150





#use this command to find out how many jobs:  find parentDir -name \*.tif | wc -l


module load matlab/R2017a
export MATLABPATH=$PWD'/code/Diana/'

#export EXPPATH='20211014_ovChargeSweep/Vuldi'
#export EXPPATH='20211013_newFilterOvnights/AqLS/100msExp'
export EXPPATH='20200730_32vs37/VuldiExp'
export CODEPATH=$MATLABPATH$EXPPATH 


export PARENTDIR=$PWD'/'$EXPPATH
export CELLSPATH=$PARENTDIR'/cellSelec.mat' 
export FILES2FIXPATH=$CODEPATH'/files2Fix.mat'




#fname=$1
#ARGS="$1"
#echo "Script: $1"
let t=$SLURM_ARRAY_TASK_ID

#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); trackingParticlesprobConvFunction('$CODEPATH'); recalculatePositionsConv('$FILES2FIXPATH',1,$t); exit"
#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); trackingParticlesprobConvFunction('$CODEPATH'); recalculatePositionsConvProb('$FILES2FIXPATH',0,$t); exit"
matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); recalculatePositionsConv('$FILES2FIXPATH',1,$t); exit"
#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); recalculatePositionsConvProb('$FILES2FIXPATH',0,$t); exit"
