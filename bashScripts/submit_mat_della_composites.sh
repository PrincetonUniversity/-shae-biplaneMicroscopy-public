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
#SBATCH --array=1-45
#use this command to find out how many jobs:  find parentDir -name \*.tif | wc -l


module load matlab/R2017a
export MATLABPATH=$PWD'/code/Diana/'
#export EXPPATH='20210616_stationaryGrowthConditions/V_Gly'
#export EXPPATH='20210921_inductionTestingOvernights/AnohistInduced1ulwga2ulVybrant'
#export EXPPATH='20211014_ovChargeSweep/Vuldi'
#export EXPPATH='20210415_ChargeSweepOvernights/p+48'
#export EXPPATH='20211110_growthCurve/OD3.14'
#export EXPPATH='20211014_ovChargeSweep/P+36'
#export EXPPATH='20200730_32vs37/PFVExp'
export EXPPATH='20201103_cccp_Everywhere/PFVexp'

export CODEPATH=$MATLABPATH$EXPPATH 


#currdir=$(echo $PWD'/VK/Data/WT/20181012/Run3')

export PARENTDIR=$PWD'/'$EXPPATH
export FILESTRUCTPATHS=$PARENTDIR'/filePaths.mat'
export EXPONPATH=$PARENTDIR'/expCoeff.mat'

#parentDir=$(echo $PWD'/'$EXPPATH)
#parentDir=$(echo $PWD'/20201126_cccp/A/cropped')
#parentDir=$(echo $PWD'/20200704_37vs32/GLY_V/cropped')
#parentDir=$(echo $PWD'/20200730_32vs37/37at37/cropped')
#parentDir=$(echo $PWD'/20201103_cccp_Everywhere/Cropped/vcccp')
#parentDir=$(echo $PWD'/20201008_exponentialAqls/Cropped2')
#parentDir=$(echo $PWD'/20201023_AqLS/A/Cropped')
#parentDir=$(echo $PWD'/20191112_growthOvernight12HoursOnPad')
#parentDir=$(echo $PWD'/20191213_labeledVsUnlabeledOvernights/glucoseLabeled')
#parentDir=$(echo $PWD'/20191217_labeledVsUnlabeled2')
#parentDir =$(echo $PWD'/20181025_DVM001_experiment')
#parentDir=$(echo $PWD'/20200730_32vs37/VuldiExp')


#fname=$1
#ARGS="$1"
#echo "Script: $1"
let t=$SLURM_ARRAY_TASK_ID
#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); infile = fullfile('$currdir', '$fname'); calcgaussCurv(infile, 0.689, 50, 2000, tiffList($t)); exit"
#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); load('calibrations20200930.mat'); tiffList =chooseTiffs('$parentDir'); createComposites(f,st,stprime,tiffList,$t); exit"

#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); load('calibrations20210416.mat'); load('$EXPONPATH'); tiffList =chooseTiffsMulticolor('$PARENTDIR'); createCompositesMulticolorBckSub(f,st,stprime,tiffList,$t,avgCoff,'$FILESTRUCTPATHS'); exit"
#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); load('calibrations20210416.mat'); load('$EXPONPATH'); extractFilesMulticolor('$PARENTDIR'); tiffList =chooseTiffsMulticolor('$PARENTDIR'); createCompositesMulticolorBckSub(f,st,stprime,tiffList,$t,avgCoff,'$FILESTRUCTPATHS'); exit"
matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); load('calibrations20210622.mat'); extractFilesMulticolor('$PARENTDIR'); tiffList =chooseTiffsMulticolor('$PARENTDIR'); createCompositesMulticolor(f,st,stprime,tiffList,$t); exit"
#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); load('calibrations20210622.mat');  tiffList =chooseTiffs('$PARENTDIR'); createComposites(f,st,stprime,tiffList,$t); exit"
#matlab -nodisplay -nodesktop -r "cd('$MATLABPATH'); load('calibrations20210416.mat'); tiffList =chooseTiffsMulticolor('$PARENTDIR'); createCompositesMulticolor(f,st,stprime,tiffList,$t); exit"