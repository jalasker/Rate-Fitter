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
#SBATCH --array=1-25
####SBATCH --array=1-43



cd /project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/

if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataCCScaleMod' 'BlindDataCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModBlind.txt  StandardCuts.txt  >& DiagFiles/RateOutputBlindDataCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData-K10CCScaleMod' 'BlindData-K10CCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10.DUMP' 0 0 0.00 GlobalParams_CCScaleModBlind.txt  StandardCuts.txt  >& DiagFiles/RateOutputBlindData-K10CCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindData-J17bCCScaleMod' 'BlindData-J17bCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b.DUMP' 0 0 0.00 GlobalParams_CCScaleModBlind.txt  StandardCuts.txt  >& DiagFiles/RateOutputBlindData-J17bCCScaleMod.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-K10-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-K10DataCCScaleMod' 'SameBetaSameK-K10DataCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleMod.txt  StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-K10DataCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-J17b-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-J17bDataCCScaleMod' 'SameBetaSameK-J17bDataCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleMod.txt  StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-J17bDataCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 6 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-K10MCCCScaleMod' 'SameBetaSameK-K10MCCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-K10.DUMP' 0 0 0.00 GlobalParams_CCScaleMod.txt  StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-K10MCCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 7 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK-J17bMCCCScaleMod' 'SameBetaSameK-J17bMCCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-J17b.DUMP' 0 0 0.00 GlobalParams_CCScaleMod.txt  StandardCuts.txt  >& DiagFiles/RateOutputSameBetaSameK-J17bMCCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 8 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightNEpochCCScaleMod' 'BlindDataTightNEpochCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModBlind.txt  TightNEpoch.txt  >& DiagFiles/RateOutputBlindDataTightNEpochCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 9 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightFITPROBCCScaleMod' 'BlindDataTightFITPROBCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModBlind.txt  TightFITPROB.txt  >& DiagFiles/RateOutputBlindDataTightFITPROBCCScaleMod.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 10 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightSNRMAXCCScaleMod' 'BlindDataTightSNRMAXCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModBlind.txt  TightSNRMAX.txt  >& DiagFiles/RateOutputBlindDataTightSNRMAXCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 11 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightNEpochCCScaleMod' 'SameBetaSameKTightNEpochCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleMod.txt  TightNEpoch.txt  >& DiagFiles/RateOutputSameBetaSameKTightNEpochCCScaleMod.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 12 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightFITPROBCCScaleMod' 'SameBetaSameKTightFITPROBCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleMod.txt  TightFITPROB.txt  >& DiagFiles/RateOutputSameBetaSameKTightFITPROBCCScaleMod.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 13 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightSNRMAXCCScaleMod' 'SameBetaSameKTightSNRMAXCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleMod.txt  TightSNRMAX.txt  >& DiagFiles/RateOutputSameBetaSameKTightSNRMAXCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 14 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZCCScaleMod' 'BlindDataMidZCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZBlind.txt  StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindDataCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 15 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsMidZCCScaleMod' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_ContamMidZCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_CCScaleModMidZ.txt  StandardCuts.txt  >& DiagFiles/RateOutputMidZ100SimsSameBetaSameKContamNNCutCCScaleMod.out
fi






if [ $SLURM_ARRAY_TASK_ID  -eq 16 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightNEpochMidZCCSCaleMod' 'BlindDataTightNEpochMidZCCSCaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZBlind.txt  TightNEpoch.txt  >& DiagFiles/RateOutputBlindDataTightNEpochMidZCCSCaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 17 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightFITPROBMidZCCScaleMod' 'BlindDataTightFITPROBMidZCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZBlind.txt  TightFITPROB.txt  >& DiagFiles/RateOutputBlindDataTightFITPROBMidZCCScaleMod.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 18 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightSNRMAXMidZCCScaleMod' 'BlindDataTightSNRMAXMidZCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZBlind.txt  TightSNRMAX.txt  >& DiagFiles/RateOutputBlindDataTightSNRMAXMidZCCScaleMod.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 19 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightNEpochMidZCCScaleMod' 'SameBetaSameKTightNEpochMidZCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZ.txt  TightNEpoch.txt  >& DiagFiles/RateOutputSameBetaSameKTightNEpochMidZCCScaleMod.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 20 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightFITPROBMidZCCScaleMod' 'SameBetaSameKTightFITPROBMidZCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZ.txt  TightFITPROB.txt  >& DiagFiles/RateOutputSameBetaSameKTightFITPROBMidZCCScaleMod.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 21 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightSNRMAXMidZCCScaleMod' 'SameBetaSameKTightSNRMAXMidZCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZ.txt  TightSNRMAX.txt  >& DiagFiles/RateOutputSameBetaSameKTightSNRMAXMidZCCScaleMod.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 22 ]
then
python -u  Rate_Fitter.py 'NoFudge/outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/NoFudge/outFit_datalike_Contam_NoFudge/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_NoFudge/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataNoFudgeCCScaleMod' 'BlindDataFudgeCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_NoFudge/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModBlind.txt  StandardCuts.txt  >& DiagFiles/RateOutputBlindDataNoFudgeCCScaleMod.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 23 ]
then
python -u  Rate_Fitter.py 'NoFudge/outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/NoFudge/outFit_datalike_Contam_NoFudge/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_NoFudge/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZNoFudgeCCScaleMod' 'BlindDataMidZNoFudgeCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam_NoFudge/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZBlind.txt  StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindDataNoFudgeCCScaleMod.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 24 ]
then
python -u  Rate_Fitter.py 'DES3YRFudge/outFit_data_DES3YR/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/DES3YRFudge/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataDES3YRFudgeCCScaleMod' 'BlindDataDES3YRFudgeCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR.DUMP' 0 0 0.00 GlobalParams_CCScaleModBlind.txt  StandardCuts.txt  >& DiagFiles/RateOutputBlindDataDES3YRFudgeCCScaleMod.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 25 ]
then
python -u  Rate_Fitter.py 'DES3YRFudge/outFit_data_DES3YR/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/DES3YRFudge/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataMidZDES3YRFudgeCCScaleMod' 'BlindDataMidZDES3YRFudgeCCScaleMod' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam-DES3YR.DUMP' 0 0 0.00 GlobalParams_CCScaleModMidZBlind.txt  StandardCuts.txt  >& DiagFiles/RateOutputMidZBlindDataDES3YRFudgeCCScaleMod.out
fi



