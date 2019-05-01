#!/bin/bash
#SBATCH --job-name=RateFitter
#SBATCH --output=DiagFiles/RateFitter%a.out
#SBATCH --error=DiagFiles/RateFitter%a.out
#SBATCH --time=00:10:00
#SBATCH --partition=debug
##SBATCH --account=kicp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
#SBATCH --mem-per-cpu=8000
##SBATCH --exclude=midway172
#SBATCH --array=71-85
###SBATCH --array=3-14,19-46
###SBATCH --array=22,23
cd /project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/

if [ $SLURM_ARRAY_TASK_ID  -eq 2 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta26.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta26SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 3 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 4 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameKContamNNCut_CheatTypeZ.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 5 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameKContamNNCut_CheatType.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 6 ]
then
python -u  Rate_Fitter.py 'Combined100Sims.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameKContamNNCut_CheatZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 7 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta22.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.21 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta22SameKContamNNCut_CheatTypeZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 8 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta22.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.21 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta22SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 9 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta22.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.21 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta22SameKContamNNCut_CheatType.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 10 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta23.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta23SameKContamNNCut_CheatZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 11 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta23.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta23SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 12 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta23.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta23SameKContamNNCut_CheatType.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 13 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta23.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100SimsCombined' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta23SameKContamNNCut_CheatZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 15 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeDataCat_001.txt' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeSimCat.txt' 'TestSameKSameBeta/FakeCatalog' 'FakeCatalogData' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeGenCat.txt' 0 0 >& RateOutputFakeCatalogs.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 16 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeDataCat_001.txt' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeSimCat.txt' 'TestSameKSameBeta/FakeCatalog' 'FakeCatalogData' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeGenCat.txt' 0 1 >& RateOutputFakeCatalogs_CheatZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 17 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeDataCat_001.txt' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeSimCat.txt' 'TestSameKSameBeta/FakeCatalog' 'FakeCatalogData' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeGenCat.txt' 1 0 >& RateOutputFakeCatalogs_CheatType.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 18 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeDataCat_001.txt' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeSimCat.txt' 'TestSameKSameBeta/FakeCatalog' 'FakeCatalogData' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/FakeGenCat.txt' 1 1 >& RateOutputFakeCatalogs_CheatTypeZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 19 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_CheatTypeZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 20 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_CheatType.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 21 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 22 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.1 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 23 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.21 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 24 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.21 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 25 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.21 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 26 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta22/JLDES_R0_1-7E-5_Beta_2-2_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta22SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.21 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta22SameKContamNNCut_CheatTypeZ.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 27 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 28 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 29 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 30 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta23/JLDES_R0_1-7E-5_Beta_2-3_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta23SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-3_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.31 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta23SameKContamNNCut_CheatTypeZ.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 31 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta216/JLDES_R0_1-7E-5_Beta_2-16_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta216SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-16_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.16 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta216SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 32 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta216/JLDES_R0_1-7E-5_Beta_2-16_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta216SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-16_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.16 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta216SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 33 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta216/JLDES_R0_1-7E-5_Beta_2-16_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta216SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-16_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.16 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta216SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 34 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta216/JLDES_R0_1-7E-5_Beta_2-16_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta216SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-2_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.16 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta216SameKContamNNCut_CheatTypeZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 35 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta24/JLDES_R0_1-7E-5_Beta_2-4_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta24SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-4_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.41 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta24SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 36 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta24/JLDES_R0_1-7E-5_Beta_2-4_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta24SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-4_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.41 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta24SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 37 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta24/JLDES_R0_1-7E-5_Beta_2-4_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta24SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-4_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.41 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta24SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 38 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta24/JLDES_R0_1-7E-5_Beta_2-4_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta24SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-4_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.41 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta24SameKContamNNCut_CheatTypeZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 39 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta25/JLDES_R0_1-7E-5_Beta_2-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta25SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-5_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.51 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta25SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 40 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta25/JLDES_R0_1-7E-5_Beta_2-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta25SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-5_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.51 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta25SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 41 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta25/JLDES_R0_1-7E-5_Beta_2-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta25SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-5_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.51 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta25SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 42 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta25/JLDES_R0_1-7E-5_Beta_2-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta25SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-5_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.51 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta25SameKContamNNCut_CheatTypeZ.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 43 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 44 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 45 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 46 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameKContamNNCut_CheatTypeZ.out
fi



if [ $SLURM_ARRAY_TASK_ID  -eq 47 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-0_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.01 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 48 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-0_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.01 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 49 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-0_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.01 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 50 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta20/JLDES_R0_1-7E-5_Beta_2-0_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta20SameK100Sims' 'JLDES_R0_1-7E-5_Beta_2-0_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.01 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta20SameKContamNNCut_CheatTypeZ.out
fi






if [ $SLURM_ARRAY_TASK_ID  -eq 51 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-9_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.91 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 52 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-9_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.91 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 53 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-9_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.91 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 54 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta19/JLDES_R0_1-7E-5_Beta_1-9_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta19SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-9_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.91 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta19SameKContamNNCut_CheatTypeZ.out
fi









if [ $SLURM_ARRAY_TASK_ID  -eq 55 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.81 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 56 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.81 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 57 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.81 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 58 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta18/JLDES_R0_1-7E-5_Beta_1-8_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta18SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-8_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.81 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta18SameKContamNNCut_CheatTypeZ.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 59 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta17/JLDES_R0_1-7E-5_Beta_1-7_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta17SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-7_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.71 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta17SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 60 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta17/JLDES_R0_1-7E-5_Beta_1-7_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta17SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-7_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.71 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta17SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 61 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta17/JLDES_R0_1-7E-5_Beta_1-7_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta17SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-7_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.71 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta17SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 62 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta17/JLDES_R0_1-7E-5_Beta_1-7_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta17SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-7_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.71 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta17SameKContamNNCut_CheatTypeZ.out
fi









if [ $SLURM_ARRAY_TASK_ID  -eq 63 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 64 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 65 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 66 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameKContamNNCut_CheatTypeZ.out
fi
























if [ $SLURM_ARRAY_TASK_ID  -eq 67 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta15/JLDES_R0_1-7E-5_Beta_1-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.51 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta15SameKContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 68 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta15/JLDES_R0_1-7E-5_Beta_1-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.51 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta15SameKContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 69 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta15/JLDES_R0_1-7E-5_Beta_1-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.51 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta15SameKContamNNCut_CheatType.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 70 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_Contam_Beta15/JLDES_R0_1-7E-5_Beta_1-5_Datasize_Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta15SameK100Sims' 'JLDES_R0_1-7E-5_Beta_1-5_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.51 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta15SameKContamNNCut_CheatTypeZ.out
fi




if [ $SLURM_ARRAY_TASK_ID  -eq 71 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 72 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.11 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_CheatType.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 73 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.11 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 74 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam/JLDES_R0_1-7E-5_Beta_2-1_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.11 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsSameBetaSameK2ContamNNCut_CheatTypeZ.out
fi





if [ $SLURM_ARRAY_TASK_ID  -eq 75 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameK2ContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 76 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameK2ContamNNCut_CheatType.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 77 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameK2ContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 78 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta26/JLDES_R0_1-7E-5_Beta_2-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta26SameK2ContamNNCut_CheatTypeZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 79 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameK2ContamNNCut.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 80 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 0 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameK2ContamNNCut_CheatType.out
fi

if [ $SLURM_ARRAY_TASK_ID  -eq 81 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 1 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameK2ContamNNCut_CheatZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 82 ]
then
python -u  Rate_Fitter.py '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datasize_2Contam_Beta16/JLDES_R0_1-7E-5_Beta_1-6_Datasize_2Contam-00{0:02d}/FITOPT000+SALT2mu.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta16SameK2Contam100Sims' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 1 1 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsBeta16SameK2ContamNNCut_CheatTypeZ.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 83 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta26_2Contam.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100SimsCombined_2Contam' 'JLDES_R0_1-7E-5_Beta_2-6_DataLikePhotZ_2Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta26SameK2ContamNNCut.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 84 ]
then
python -u  Rate_Fitter.py 'Combined100Sims_2Contam.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/SameBetaSameK100SimsCombined_2Contam' 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ_2Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 2.11 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedSameBetaSameK2ContamNNCut.out
fi


if [ $SLURM_ARRAY_TASK_ID  -eq 85 ]
then
python -u  Rate_Fitter.py 'Combined100SimsBeta16_2Contam.fitres' '/project/rkessler/SURVEYS/DES/USERS/jlasker/Rate_Analysis/outFit_datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/FITOPT000+SALT2mu.fitres' 'DiagFiles/Beta26SameK100SimsCombined_2Contam' 'JLDES_R0_1-7E-5_Beta_1-6_DataLikePhotZ_2Contam100SimsCombined' '/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam/JLDES_R0_1-7E-5_Beta_2-1_Datalike_Contam.DUMP' 0 0 1.61 GlobalParams.txt StandardCuts.txt  >& DiagFiles/RateOutput100SimsCombinedBeta16SameK2ContamNNCut.out
fi


