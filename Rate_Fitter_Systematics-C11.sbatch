#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a.out
#SBATCH --error=DiagFiles/RateFitter%a.out
#SBATCH --time=00:60:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
###SBATCH --qos debug
#SBATCH --mem-per-cpu=24000
#SBATCH --array=1-43



cd /project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/

if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData-C11' 'BlindData-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-K10-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData-K10-C11' 'BlindData-K10-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData-K10-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-J17b-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData-J17b-C11' 'BlindData-J17b-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b-C11.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData-J17b-C11.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-K10-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-K10-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-K10-C11Data' 'SameBetaSameK-K10-C11Data' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-K10-C11Data.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-J17b-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-J17b-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-J17b-C11Data' 'SameBetaSameK-J17b-C11Data' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-J17b-C11Data.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 6 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-K10-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-K10-C11MC' 'SameBetaSameK-K10-C11MC' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-K10-C11MC.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 7 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-J17b-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-J17b-C11MC' 'SameBetaSameK-J17b-C11MC' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b-C11.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-J17b-C11MC.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 8 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCut-C11s/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKWeakCuts-C11' 'SameBetaSameKWeakCuts-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameKWeakCuts-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 12 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKCosmologyCuts-C11' 'SameBetaSameKCosmologyCuts-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams.txt CosmologyCuts.txt  >& DiagFiles/RateOutputSameBetaSameKCosmologyCuts-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 13 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataWeakCuts-C11' 'BlindDataWeakCuts-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataWeakCuts-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 17 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataCosmologyCuts-C11' 'BlindDataCosmologyCuts-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_Blind.txt CosmologyCuts.txt  >& DiagFiles/RateOutputBlindDataCosmologyCuts-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 18 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightNEpoch-C11' 'BlindDataTightNEpoch-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_Blind.txt TightNEpoch.txt  >& DiagFiles/RateOutputBlindDataTightNEpoch-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 19 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightFITPROB-C11' 'BlindDataTightFITPROB-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_Blind.txt TightFITPROB.txt  >& DiagFiles/RateOutputBlindDataTightFITPROB-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 20 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightSNRMAX-C11' 'BlindDataTightSNRMAX-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_Blind.txt TightSNRMAX.txt  >& DiagFiles/RateOutputBlindDataTightSNRMAX-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 21 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightNEpoch-C11' 'SameBetaSameKTightNEpoch-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams.txt TightNEpoch.txt  >& DiagFiles/RateOutputSameBetaSameKTightNEpoch-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 22 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightFITPROB-C11' 'SameBetaSameKTightFITPROB-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams.txt TightFITPROB.txt  >& DiagFiles/RateOutputSameBetaSameKTightFITPROB-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 23 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightSNRMAX-C11' 'SameBetaSameKTightSNRMAX-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams.txt TightSNRMAX.txt  >& DiagFiles/RateOutputSameBetaSameKTightSNRMAX-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 24 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZ-C11' 'BlindDataMidZ-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindData-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 25 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsMidZ-C11' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 2.11 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsSameBetaSameKContamNNCut-C11.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 26 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZ_WeakCuts-C11' 'BlindDataMidZ_WeakCuts-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindData_WeakCuts-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 27 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_WeakCuts-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCut-C11s/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsMidZ_WeakCuts-C11' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ_WeakCuts-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsSameBetaSameKContamNNCut_WeakCuts-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 28 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataZSys-C11' 'BlindDataZSys-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_ZSys_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataZSys.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 29 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKZSys-C11' 'SameBetaSameKZSys-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_ZSys.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameKZSys-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 30 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightNEpochMidZ' 'BlindDataTightNEpochMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt TightNEpoch.txt  >& DiagFiles/RateOutputBlindDataTightNEpochMidZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 31 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightFITPROBMidZ-C11' 'BlindDataTightFITPROBMidZ-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt TightFITPROB.txt  >& DiagFiles/RateOutputBlindDataTightFITPROBMidZ-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 32 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightSNRMAXMidZ-C11' 'BlindDataTightSNRMAXMidZ-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt TightSNRMAX.txt  >& DiagFiles/RateOutputBlindDataTightSNRMAXMidZ-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 33 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightNEpochMidZ-C11' 'SameBetaSameKTightNEpochMidZ-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_MidZ.txt TightNEpoch.txt  >& DiagFiles/RateOutputSameBetaSameKTightNEpochMidZ-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 34 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightFITPROBMidZ-C11' 'SameBetaSameKTightFITPROBMidZ-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_MidZ.txt TightFITPROB.txt  >& DiagFiles/RateOutputSameBetaSameKTightFITPROBMidZ-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 35 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightSNRMAXMidZ-C11' 'SameBetaSameKTightSNRMAXMidZ-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_MidZ.txt TightSNRMAX.txt  >& DiagFiles/RateOutputSameBetaSameKTightSNRMAXMidZ-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 36 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData1CCZBin-C11' 'BlindData1CCZBin-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_Blind_1CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData1CCZBin-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 37 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK1CCZBin-C11' 'SameBetaSameK1CCZBin-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_1CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK1CCZBin-C11.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 38 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts-C11/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts-K10-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataWeakCuts-K10-C11' 'BlindDataWeakCuts-K10-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataWeakCuts-K10-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 39 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts-K10-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKWeakCuts-K10-C11MC' 'SameBetaSameKWeakCuts-K10-C11MC' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10-C11.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameKWeakCuts-K10-C11MC.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 40 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData2CCZBin-C11' 'BlindData2CCZBin-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_Blind_2CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData2CCZBin-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 41 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2CCZBin-C11' 'SameBetaSameK2CCZBin-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_2CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK2CCZBin-C11.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 42 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData4CCZBin-C11' 'BlindData4CCZBin-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_Blind_4CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData4CCZBin-C11.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 43 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-C11-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK4CCZBin-C11' 'SameBetaSameK4CCZBin-C11' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-C11.DUMP' 0 0 0.00 GlobalParams_4CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK4CCZBin-C11.out
fi

