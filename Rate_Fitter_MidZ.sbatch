#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a.out
#SBATCH --error=DiagFiles/RateFitter%a.out
#SBATCH --time=00:60:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
##SBATCH --qos debug
#SBATCH --mem-per-cpu=24000
#SBATCH --array=1-9

### 21 for same beta, same k, same contam
### 71 for '' '', 2x contam
###

cd /project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/

if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZ' 'BlindDataMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindData.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMCBeta18MidZ' 'BlindDataMCBeta18MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindDataMCBeta18.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_Beta05/JLDES_R0_1-7E-5_Beta_0-5_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMCBeta05MidZ' 'BlindDataMCBeta05' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_0-5_Datalike_Contam/JLDES_R0_1-7E-5_Beta_0-5_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindDataMCBeta05.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsMidZ' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsSameBetaSameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombinedMidZ' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombinedMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsCombinedSameBetaSameKContamNNCut.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 6 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100SimsMidZ' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.91 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsBeta19SameKContamNNCut.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 7 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100SimsMidZ' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.01 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsBeta20SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 8 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100SimsMidZ' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.21 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsBeta22SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 9 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23S/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsMidZ' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.31 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsBeta23SameKContamNNCut.out
fi