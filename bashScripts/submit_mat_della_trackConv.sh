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
#-176
#-153
#-389
#-15 p-30
#-176 P-18
#-46 +9
#-30 +15
#14 +48




#use this command to find out how many jobs:  find parentDir -name \*.tif | wc -l


module load matlab/R2017a
export MATLABPATH=$PWD'/code/Diana/'
#export EXPPATH='20210415_ChargeSweepOvernights/p+15'
#export EXPPATH='20200911_charge_exponential/+15/reorganized'
#export EXPPATH='20200908_charge_exponential/-18/reorganized'
#export EXPPATH='20210415_ChargeSweepOvernights/V'
#export EXPPATH='20211014_ovChargeSweep/PFV'
#export EXPPATH='20211013_newFilterOvnights/Vuldi'
#export EXPPATH='20210415_ChargeSweepOvernights/P+9'
#export EXPPATH='20211014_ovChargeSweep/P+9'
export EXPPATH='20200730_32vs37/VuldiExp'



export CODEPATH=$MATLABPATH$EXPPATH 

#currdir=$(echo $PWD'/VK/Data/WT/20181012/Run3')


export PARENTDIR=$PWD'/'$EXPPATH
export CELLSPATH=$PARENTDIR'/cellSelec.mat' 
export FILESTRUCTPATHS=$PARENTDIR'/filePaths.mat'




#fname=$1
#ARGS="$1"
#echo "Script: $1"
let t=$SLURM_ARRAY_TASK_ID

#calibrations20210416.mat
#calibrations20200930.mat
#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH');  extractFilesMulticolor('$PARENTDIR'); [jFile,iiCell]=whichCell('$CELLSPATH',$t);  convolutionTrackingNucFirstExp('$CELLSPATH','$FILESTRUCTPATHS', jFile,iiCell); exit"

#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH');  [jFile,iiCell]=whichCell('$CELLSPATH',$t);  convolutionTrackingNucFirstExp('$CELLSPATH','$FILESTRUCTPATHS', jFile,iiCell); exit"
matlab -nodisplay -nodesktop -r "cd('$MATLABPATH');  [jFile,iiCell]=whichCell('$CELLSPATH',$t);  convolutionTrackingNucFirstDrawnMask('$CELLSPATH','$FILESTRUCTPATHS', jFile,iiCell); exit"

#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH');  [jFile,iiCell]=whichCell('$CELLSPATH',$t);  convolutionTrackingNucFirst('$CELLSPATH','$FILESTRUCTPATHS', jFile,iiCell); exit"
