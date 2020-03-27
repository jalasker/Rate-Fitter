#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a.out
#SBATCH --error=DiagFiles/RateFitter%a.out
#SBATCH --time=00:15:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
#SBATCH --qos debug
#SBATCH --mem-per-cpu=24000
#SBATCH --array=1

### 21 for same beta, same k, same contam
### 71 for '' '', 2x contam
###

cd /project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/

if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData' 'BlindData' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMCBeta18' 'BlindDataMCBeta18' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataMCBeta18.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_Beta05/JLDES_R0_1-7E-5_Beta_0-5_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMCBeta05' 'BlindDataMCBeta05' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_0-5_Datalike_Contam/JLDES_R0_1-7E-5_Beta_0-5_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataMCBeta05.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData_Broken' 'BlindData_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py 'outFit_data_NoAGN/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataAGNCut' 'BlindDataAGNCut' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataAGNCut.out
fi