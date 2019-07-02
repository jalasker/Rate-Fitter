#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a.out
#SBATCH --error=DiagFiles/RateFitter%a.out
#SBATCH --time=00:60:00
#SBATCH --partition=broadwl-lc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
###SBATCH --qos debug
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1

### 21 for same beta, same k, same contam
### 71 for '' '', 2x contam
###

cd /project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/

if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData.out
fi
