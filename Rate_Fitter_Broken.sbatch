#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a_Broken.out
#SBATCH --error=DiagFiles/RateFitter%a_Broken.out
#SBATCH --time=00:60:00
#SBATCH --partition=broadwl-lc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
###SBATCH --qos debug
#SBATCH --mem-per-cpu=24000
#SBATCH --array=21,71,55,67,94-97

###SBATCH --array=96-97
###SBATCH --array=21,71,55,67,94,95

### 21 for same beta, same k, same contam
### 71 for '' '', 2x contam
###

cd /project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/

if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta26.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta26SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.1 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameKContamNNCut_CheatTypeZ_Broken.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.1 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameKContamNNCut_CheatType_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 6 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.1 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameKContamNNCut_CheatZ_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 7 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta22.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.21 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta22SameKContamNNCut_CheatTypeZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 8 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta22.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.21 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta22SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 9 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta22.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.21 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta22SameKContamNNCut_CheatType_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 10 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta23.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta23SameKContamNNCut_CheatZ_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 11 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta23.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta23SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 12 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta23.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta23SameKContamNNCut_CheatType_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 13 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta23.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsCombined_Broken' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta23SameKContamNNCut_CheatZ_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 15 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeDataCat_001.txt' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeSimCat.txt' 'TestSameKSameBeta/FakeCatalog' 'FakeCatalogData' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeGenCat.txt' 0 0 >& RateOutputFakeCatalogs_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 16 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeDataCat_001.txt' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeSimCat.txt' 'TestSameKSameBeta/FakeCatalog' 'FakeCatalogData' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeGenCat.txt' 0 1 >& RateOutputFakeCatalogs_CheatZ_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 17 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeDataCat_001.txt' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeSimCat.txt' 'TestSameKSameBeta/FakeCatalog' 'FakeCatalogData' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeGenCat.txt' 1 0 >& RateOutputFakeCatalogs_CheatType_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 18 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeDataCat_001.txt' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeSimCat.txt' 'TestSameKSameBeta/FakeCatalog' 'FakeCatalogData' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeGenCat.txt' 1 1 >& RateOutputFakeCatalogs_CheatTypeZ_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 19 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.1 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_CheatTypeZ_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 20 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.1 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_CheatType_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 21 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 22 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.1 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 23 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.21 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 24 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.21 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 25 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.21 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 26 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.21 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut_CheatTypeZ_Broken.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 27 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 28 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 29 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 30 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut_CheatTypeZ_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 31 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta216/JLDES_R0_1-7E-5_Beta_2-16_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta216SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-16_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.16 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta216SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 32 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta216/JLDES_R0_1-7E-5_Beta_2-16_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta216SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-16_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.16 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta216SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 33 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta216/JLDES_R0_1-7E-5_Beta_2-16_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta216SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-16_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.16 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta216SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 34 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta216/JLDES_R0_1-7E-5_Beta_2-16_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta216SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.16 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta216SameKContamNNCut_CheatTypeZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 35 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta24/JLDES_R0_1-7E-5_Beta_2-4_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta24SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-4_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.41 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta24SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 36 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta24/JLDES_R0_1-7E-5_Beta_2-4_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta24SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-4_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.41 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta24SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 37 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta24/JLDES_R0_1-7E-5_Beta_2-4_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta24SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-4_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.41 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta24SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 38 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta24/JLDES_R0_1-7E-5_Beta_2-4_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta24SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-4_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.41 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta24SameKContamNNCut_CheatTypeZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 39 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta25/JLDES_R0_1-7E-5_Beta_2-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta25SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-5_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta25SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 40 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta25/JLDES_R0_1-7E-5_Beta_2-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta25SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-5_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta25SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 41 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta25/JLDES_R0_1-7E-5_Beta_2-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta25SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-5_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta25SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 42 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta25/JLDES_R0_1-7E-5_Beta_2-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta25SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-5_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta25SameKContamNNCut_CheatTypeZ_Broken.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 43 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 44 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 45 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 46 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameKContamNNCut_CheatTypeZ_Broken.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 47 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-0_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.01 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 48 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-0_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.01 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 49 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-0_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.01 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 50 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-0_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.01 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut_CheatTypeZ_Broken.out
fi






if [ $SLURM_ARRAY_TASK_ID  -eq 51 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-9_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.91 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 52 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-9_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.91 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 53 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-9_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.91 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 54 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-9_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.91 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut_CheatTypeZ_Broken.out
fi









if [ $SLURM_ARRAY_TASK_ID  -eq 55 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.81 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 56 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.81 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 57 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.81 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 58 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.81 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameKContamNNCut_CheatTypeZ_Broken.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 59 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta17/JLDES_R0_1-7E-5_Beta_1-7_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta17SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-7_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.71 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta17SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 60 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta17/JLDES_R0_1-7E-5_Beta_1-7_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta17SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-7_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.71 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta17SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 61 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta17/JLDES_R0_1-7E-5_Beta_1-7_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta17SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-7_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.71 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta17SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 62 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta17/JLDES_R0_1-7E-5_Beta_1-7_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta17SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-7_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.71 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta17SameKContamNNCut_CheatTypeZ_Broken.out
fi









if [ $SLURM_ARRAY_TASK_ID  -eq 63 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 64 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 65 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 66 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameKContamNNCut_CheatTypeZ_Broken.out
fi
























if [ $SLURM_ARRAY_TASK_ID  -eq 67 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta15/JLDES_R0_1-7E-5_Beta_1-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta15SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 68 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta15/JLDES_R0_1-7E-5_Beta_1-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta15SameKContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 69 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta15/JLDES_R0_1-7E-5_Beta_1-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta15SameKContamNNCut_CheatType_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 70 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta15/JLDES_R0_1-7E-5_Beta_1-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta15SameKContamNNCut_CheatTypeZ_Broken.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 71 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 72 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.11 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_CheatType_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 73 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.11 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 74 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.11 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_CheatTypeZ_Broken.out
fi





if [ $SLURM_ARRAY_TASK_ID  -eq 75 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameK2ContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 76 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameK2ContamNNCut_CheatType_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 77 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameK2ContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 78 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameK2ContamNNCut_CheatTypeZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 79 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameK2ContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 80 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameK2ContamNNCut_CheatType_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 81 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameK2ContamNNCut_CheatZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 82 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameK2ContamNNCut_CheatTypeZ_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 83 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta26_2Contam.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100SimsCombined_2Contam_Broken' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_2Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta26SameK2ContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 84 ]
then
python -u  Rate_Fitter.py 'Combined100Sims_2Contam.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined_2Contam_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameK2ContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 85 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta16_2Contam.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100SimsCombined_2Contam_Broken' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_2Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.61 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta16SameK2ContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 86 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta203/JLDES_R0_1-7E-5_Beta_2-03_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta203SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.03 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta203SameK2ContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 87 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta207/JLDES_R0_1-7E-5_Beta_2-07_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta207SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.07 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta207SameK2ContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 88 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta215/JLDES_R0_1-7E-5_Beta_2-15_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta215SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.15 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta215SameK2ContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 89 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta219/JLDES_R0_1-7E-5_Beta_2-19_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta219SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.19 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta219SameK2ContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 90 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.21 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameK2ContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 91 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.31 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameK2ContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 92 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta218/JLDES_R0_1-7E-5_Beta_1-8_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.81 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameK2ContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 93 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK2Contam100Sims_Broken' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.91 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameK2ContamNNCut_Broken.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 94 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta15_Contam.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100SimsCombined_Contam_Broken' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.51 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta15SameKContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 95 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta18_Contam.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100SimsCombined_Contam_Broken' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam100SimsCombined_Broken' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.81 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta18SameKContamNNCut_Broken.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 96 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta18_Contam.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/MCBeta18DataBeta18SameK100SimsCombined_Contam_Broken' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam100SimsCombinedMCBeta18' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam.DUMP' 0 0 1.81 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedMCBeta18DataBeta18SameKContamNNCut_Broken.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 97 ]
then
python -u  Rate_Fitter.py '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project2/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/MCBeta18DataBeta21SameK100Sims_Contam_Broken' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam100SimsMCBeta18' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam/JLDES_R0_1-7E-5_Beta_1-8_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams_Broken.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedMCBeta18DataBeta21SameKContamNNCut_Broken.out
fi
