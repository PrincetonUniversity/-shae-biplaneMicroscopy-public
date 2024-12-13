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
#SBATCH --array=1-13

#use this command to find out how many jobs:  find parentDir -name \*.tif | wc -l


module load matlab/R2017a
export MATLABPATH=$PWD'/code/Diana/'
#export EXPPATH='20210616_stationaryGrowthConditions/V_Gly'
#export EXPPATH='20210921_inductionTestingOvernights/AnohistInduced1ulwga2ulVybrant'
export EXPPATH='20200730_32vs37/PFVExp'
export CODEPATH=$MATLABPATH$EXPPATH 


#currdir=$(echo $PWD'/VK/Data/WT/20181012/Run3')

export PARENTDIR=$PWD'/'$EXPPATH
export FILESTRUCTPATHS=$PARENTDIR'/filePaths.mat'
export EXPONPATH=$PARENTDIR'/expCoeff.mat'



#fname=$1
#ARGS="$1"
#echo "Script: $1"
let t=$SLURM_ARRAY_TASK_ID

matlab -nodisplay -nodesktop -r "cd('$MATLABPATH');  divideImages('$PARENTDIR', $t); exit"
