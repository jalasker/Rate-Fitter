#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a.out
#SBATCH --error=DiagFiles/RateFitter%a.out
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
###SBATCH --qos debug
#SBATCH --mem-per-cpu=24000
#SBATCH --array=13-24
####SBATCH --array=7-12



cd /project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/


if [ $SLURM_ARRAY_TASK_ID  -eq 1 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightSNRMAX' 'BlindDataTightSNRMAX' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt TightSNRMAX.txt  >& DiagFiles/RateOutputBlindDataTightSNRMAX.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX55' 'BlindDataSNRMAX55' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt SNRMAX55.txt  >& DiagFiles/RateOutputBlindDataSNRMAX55.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX6' 'BlindDataSNRMAX6' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt SNRMAX6.txt  >& DiagFiles/RateOutputBlindDataSNRMAX6.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX65' 'BlindDataSNRMAX65' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt SNRMAX65.txt  >& DiagFiles/RateOutputBlindDataSNRMAX65.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX7' 'BlindDataSNRMAX7' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt SNRMAX7.txt  >& DiagFiles/RateOutputBlindDataSNRMAX7.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 6 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX75' 'BlindDataSNRMAX75' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_Blind.txt SNRMAX75.txt  >& DiagFiles/RateOutputBlindDataSNRMAX75.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 7 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightSNRMAX' 'SameBetaSameKTightSNRMAX' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt TightSNRMAX.txt  >& DiagFiles/RateOutputSameBetaSameKTightSNRMAX.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 8 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX55' 'SameBetaSameKSNRMAX55' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt SNRMAX55.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX55.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 9 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX6' 'SameBetaSameKSNRMAX6' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt SNRMAX6.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX6.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 10 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX65' 'SameBetaSameKSNRMAX65' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt SNRMAX65.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX65.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 11 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX7' 'SameBetaSameKSNRMAX7' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt SNRMAX7.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX7.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 12 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX75' 'SameBetaSameKSNRMAX75' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams.txt SNRMAX75.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX75.out
fi










if [ $SLURM_ARRAY_TASK_ID  -eq 13 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataTightSNRMAXMidZ' 'BlindDataTightSNRMAXMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt TightSNRMAX.txt  >& DiagFiles/RateOutputBlindDataTightSNRMAXMidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 14 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX55MidZ' 'BlindDataSNRMAX55MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt SNRMAX55.txt  >& DiagFiles/RateOutputBlindDataSNRMAX55MidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 15 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX6MidZ' 'BlindDataSNRMAX6MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt SNRMAX6.txt  >& DiagFiles/RateOutputBlindDataSNRMAX6MidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 16 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX65MidZ' 'BlindDataSNRMAX65MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt SNRMAX65.txt  >& DiagFiles/RateOutputBlindDataSNRMAX65MidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 17 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX7MidZ' 'BlindDataSNRMAX7MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt SNRMAX7.txt  >& DiagFiles/RateOutputBlindDataSNRMAX7MidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 18 ]
then
python -u  Rate_Fitter.py 'outFit_data/DESALL_forcePhoto_real_snana_fits/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/BlindDataSNRMAX75MidZ' 'BlindDataSNRMAX75MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZBlind.txt SNRMAX75.txt  >& DiagFiles/RateOutputBlindDataSNRMAX75MidZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 19 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKTightSNRMAXMidZ' 'SameBetaSameKTightSNRMAXMidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt TightSNRMAX.txt  >& DiagFiles/RateOutputSameBetaSameKTightSNRMAXMidZ.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 20 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX55MidZ' 'SameBetaSameKSNRMAX55MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt SNRMAX55.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX55MidZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 21 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX6MidZ' 'SameBetaSameKSNRMAX6MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt SNRMAX6.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX6MidZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 22 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX65MidZ' 'SameBetaSameKSNRMAX65MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt SNRMAX65.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX65MidZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 23 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX7MidZ' 'SameBetaSameKSNRMAX7MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt SNRMAX7.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX7MidZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 24 ]
then
python -u  Rate_Fitter.py 'outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-0{0:03d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameKSNRMAX75MidZ' 'SameBetaSameKSNRMAX75MidZ' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 0.00 GlobalParams_MidZ.txt SNRMAX75.txt  >& DiagFiles/RateOutputSameBetaSameKSNRMAX75MidZ.out
fi