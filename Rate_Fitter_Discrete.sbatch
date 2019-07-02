#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a.out
#SBATCH --error=DiagFiles/RateFitter%a.out
#SBATCH --time=00:15:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
#SBATCH --qos debug
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1



cd /project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/


if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_Discrete' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_Discrete.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_Discrete.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_Discrete' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_Discrete.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_Discrete.out
fi

