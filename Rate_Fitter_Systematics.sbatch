#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a_WeakCuts.out
#SBATCH --error=DiagFiles/RateFitter%a_WeakCuts.out
#SBATCH --time=06:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
###SBATCH --qos debug
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1-4


cd /project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/




if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/SystematicSimsAndTest/outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam_Sys-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/SystematicSimsAndTest/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_WeakCuts_Sys' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_WeakCuts_Sys' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys.DUMP' 0 0 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_WeakCuts_Sys.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/SystematicSimsAndTest/outFit_datasize_2Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam_Sys-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/SystematicSimsAndTest/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_WeakCuts_Sys' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_WeakCuts_Sys' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys.DUMP' 0 0 2.11 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_WeakCuts_Sys.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/SystematicSimsAndTest/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam_Sys-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/SystematicSimsAndTest/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_Sys' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Sys' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys.DUMP' 0 0 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_Sys.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/SystematicSimsAndTest/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam_Sys-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/SystematicSimsAndTest/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_Sys' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Sys' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_Sys.DUMP' 0 0 2.11 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_Sys.out
fi


