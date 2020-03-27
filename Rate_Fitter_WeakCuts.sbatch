#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a_WeakCuts.out
#SBATCH --error=DiagFiles/RateFitter%a_WeakCuts.out
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
###SBATCH --qos debug
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1-10

cd /project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/




if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_WeakCuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_WeakCuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_WeakCuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_WeakCuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19_WeakCuts/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims_WeakCuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.91 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut_WeakCuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22_WeakCuts/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims_WeakCuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.21 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut_WeakCuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23_WeakCuts/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims_WeakCuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut_WeakCuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 6 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20_WeakCuts/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims_WeakCuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.01 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut_WeakCuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 7 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_L1Cuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_L1Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams.txt Level1Cuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_L1Cuts.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 8 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_L2Cuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_L2Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams.txt Level2Cuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_L2Cuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 9 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_L3Cuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_L3Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams.txt Level3Cuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_L3Cuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 10 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_CosmoCuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_CosmoCuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams.txt CosmologyCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_CosmoCuts.out
fi