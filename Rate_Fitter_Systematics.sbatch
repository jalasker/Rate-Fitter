#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a.out
#SBATCH --error=DiagFiles/RateFitter%a.out
#SBATCH --time=00:60:00
#SBATCH --partition=broadwl-lc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
###SBATCH --qos debug
#SBATCH --mem-per-cpu=24000
#SBATCH --array=1-47
####SBATCH --array=1-43



cd /project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/

if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData' 'BlindData' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData-K10' 'BlindData-K10' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData-K10.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData-J17b' 'BlindData-J17b' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData-J17b.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-K10-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-K10Data' 'SameBetaSameK-K10Data' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-K10Data.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-J17b-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-J17bData' 'SameBetaSameK-J17bData' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-J17bData.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 6 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-K10MC' 'SameBetaSameK-K10MC' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-K10MC.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 7 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-J17bMC' 'SameBetaSameK-J17bMC' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-J17bMC.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 8 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKWeakCuts' 'SameBetaSameKWeakCuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameKWeakCuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 9 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKLevel1Cuts' 'SameBetaSameKLevel1Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt Level1Cuts.txt  >& DiagFiles/RateOutputSameBetaSameKLevel1Cuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 10 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKLevel2Cuts' 'SameBetaSameKLevel2Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt Level2Cuts.txt  >& DiagFiles/RateOutputSameBetaSameKLevel2Cuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 11 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKLevel3Cuts' 'SameBetaSameKLevel3Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt Level3Cuts.txt  >& DiagFiles/RateOutputSameBetaSameKLevel3Cuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 12 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKCosmologyCuts' 'SameBetaSameKCosmologyCuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt CosmologyCuts.txt  >& DiagFiles/RateOutputSameBetaSameKCosmologyCuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 13 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataWeakCuts' 'BlindDataWeakCuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataWeakCuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 14 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataLevel1Cuts' 'BlindDataLevel1Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt Level1Cuts.txt  >& DiagFiles/RateOutputBlindDataLevel1Cuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 15 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataLevel2Cuts' 'BlindDataLevel2Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt Level2Cuts.txt  >& DiagFiles/RateOutputBlindDataLevel2Cuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 16 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataLevel3Cuts' 'BlindDataLevel3Cuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt Level3Cuts.txt  >& DiagFiles/RateOutputBlindDataLevel3Cuts.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 17 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataCosmologyCuts' 'BlindDataCosmologyCuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt CosmologyCuts.txt  >& DiagFiles/RateOutputBlindDataCosmologyCuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 18 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightNEpoch' 'BlindDataTightNEpoch' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt TightNEpoch.txt  >& DiagFiles/RateOutputBlindDataTightNEpoch.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 19 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightFITPROB' 'BlindDataTightFITPROB' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt TightFITPROB.txt  >& DiagFiles/RateOutputBlindDataTightFITPROB.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 20 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightSNRMAX' 'BlindDataTightSNRMAX' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt TightSNRMAX.txt  >& DiagFiles/RateOutputBlindDataTightSNRMAX.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 21 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightNEpoch' 'SameBetaSameKTightNEpoch' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt TightNEpoch.txt  >& DiagFiles/RateOutputSameBetaSameKTightNEpoch.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 22 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightFITPROB' 'SameBetaSameKTightFITPROB' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt TightFITPROB.txt  >& DiagFiles/RateOutputSameBetaSameKTightFITPROB.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 23 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightSNRMAX' 'SameBetaSameKTightSNRMAX' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt TightSNRMAX.txt  >& DiagFiles/RateOutputSameBetaSameKTightSNRMAX.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 24 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZ' 'BlindDataMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindData.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 25 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsMidZ' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsSameBetaSameKContamNNCut.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 26 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZ_WeakCuts' 'BlindDataMidZ_WeakCuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindData_WeakCuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 27 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsMidZ_WeakCuts' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZ_WeakCuts' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_MidZ.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsSameBetaSameKContamNNCut_WeakCuts.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 28 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataZSys' 'BlindDataZSys' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_ZSys_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataZSys.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 29 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKZSys' 'SameBetaSameKZSys' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_ZSys.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameKZSys.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 30 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightNEpochMidZ' 'BlindDataTightNEpochMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt TightNEpoch.txt  >& DiagFiles/RateOutputBlindDataTightNEpochMidZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 31 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightFITPROBMidZ' 'BlindDataTightFITPROBMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt TightFITPROB.txt  >& DiagFiles/RateOutputBlindDataTightFITPROBMidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 32 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightSNRMAXMidZ' 'BlindDataTightSNRMAXMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt TightSNRMAX.txt  >& DiagFiles/RateOutputBlindDataTightSNRMAXMidZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 33 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightNEpochMidZ' 'SameBetaSameKTightNEpochMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt TightNEpoch.txt  >& DiagFiles/RateOutputSameBetaSameKTightNEpochMidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 34 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightFITPROBMidZ' 'SameBetaSameKTightFITPROBMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt TightFITPROB.txt  >& DiagFiles/RateOutputSameBetaSameKTightFITPROBMidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 35 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightSNRMAXMidZ' 'SameBetaSameKTightSNRMAXMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt TightSNRMAX.txt  >& DiagFiles/RateOutputSameBetaSameKTightSNRMAXMidZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 36 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData1CCZBin' 'BlindData1CCZBin' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind_1CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData1CCZBin.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 37 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK1CCZBin' 'SameBetaSameK1CCZBin' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_1CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK1CCZBin.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 38 ]
then
python -u  Rate_Fitter.py 'outFit_data_WeakCuts/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataWeakCuts-K10' 'BlindDataWeakCuts-K10' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataWeakCuts-K10.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 39 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam_WeakCuts/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_WeakCuts-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKWeakCuts-K10MC' 'SameBetaSameKWeakCuts-K10MC' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10.DUMP' 0 0 0.00 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameKWeakCuts-K10MC.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 40 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData2CCZBin' 'BlindData2CCZBin' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind_2CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData2CCZBin.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 41 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2CCZBin' 'SameBetaSameK2CCZBin' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_2CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK2CCZBin.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 42 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData4CCZBin' 'BlindData4CCZBin' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind_4CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindData4CCZBin.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 43 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK4CCZBin' 'SameBetaSameK4CCZBin' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_4CCZBin.txt StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK4CCZBin.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 44 ]
then
python -u  Rate_Fitter.py 'NoFudge/outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/NoFudge/outFit_datalike_Contam_NoFudge/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_NoFudge/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataNoFudge' 'BlindDataFudge' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_NoFudge/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataNoFudge.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 45 ]
then
python -u  Rate_Fitter.py 'NoFudge/outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/NoFudge/outFit_datalike_Contam_NoFudge/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_NoFudge/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZNoFudge' 'BlindDataMidZNoFudge' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_NoFudge/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindDataNoFudge.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 46 ]
then
python -u  Rate_Fitter.py 'DES3YRFudge/outFit_data_DES3YR/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/DES3YRFudge/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataDES3YRFudge' 'BlindDataDES3YRFudge' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR.DUMP' 0 0 0.00 GlobalParams_Blind.txt StandardCuts.txt  >& DiagFiles/RateOutputBlindDataDES3YRFudge.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 47 ]
then
python -u  Rate_Fitter.py 'DES3YRFudge/outFit_data_DES3YR/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/DES3YRFudge/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZDES3YRFudge' 'BlindDataMidZDES3YRFudge' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindDataDES3YRFudge.out
fi



