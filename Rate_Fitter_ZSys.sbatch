#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a_ZSys.out
#SBATCH --error=DiagFiles/RateFitter%a_ZSys.out
#SBATCH --time=00:90:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
###SBATCH --qos debug
#SBATCH --mem-per-cpu=32000
#SBATCH --array=1-2

### 21 for same beta, same k, same contam
### 71 for '' '', 2x contam
###

cd /project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/


if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_ZSys' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_ZSys.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_ZSys.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_ZSys' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_ZSys.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_ZSys.out
fi

