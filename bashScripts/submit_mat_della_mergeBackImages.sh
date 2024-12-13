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
#SBATCH --array=1-10

#use this command to find out how many jobs:  find parentDir -name \*.tif | wc -l


module load matlab/R2017a
export MATLABPATH=$PWD'/code/Diana/'
#export EXPPATH='20210616_stationaryGrowthConditions/V_Glu'
#export EXPPATH='20210415_ChargeSweepOvernights/P'
#export EXPPATH='20211013_newFilterOvnights/AqLS/30msExp'
#export EXPPATH='20211014_ovChargeSweep/PFV'
#export EXPPATH='20210415_ChargeSweepOvernights/P-18'
#export EXPPATH='20211014_ovChargeSweep/P+9'
#export EXPPATH='20200730_32vs37/VuldiExp'
export EXPPATH='20200730_32vs37/PFVExp'
export CODEPATH=$MATLABPATH$EXPPATH 


parentDir=$(echo $PWD'/'$EXPPATH)



let t=$SLURM_ARRAY_TASK_ID

matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); imgSubDirs = getDirectoryNames('$CODEPATH');mergeBackImages('$CODEPATH',imgSubDirs,$t); exit"