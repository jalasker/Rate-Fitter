#!/software/python-2.7-2014q3-el6-x86_64/bin/python
import SNANA_Reader as simread
import REAL_Reader as dataread
#import astropy.cosmology as cosmo
import traceback
import scipy
import scipy.stats as stats
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import Cosmology
#import emcee as MC
#import corner
import scipy.stats.mstats as mstats
import scipy.stats as stats
from sys import argv
import glob
import time
import os
import gzip
import shutil
import numpy.ma as ma
import subprocess
import iminuit as iM
from iminuit import Minuit as M
from discreteChi2Func import discreteChi2Func as chi2func



class Rate_Fitter:
    def __init__(self, realfilename, realName, simfilename, simName, simgenfilename, MCBeta, zmin=0.1, zmax=1.20 , simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95,  Rate_Model = 'powerlaw',  cheatType = False, cheatZ = False, cheatCCSub = False, cheatCCScale = False, cuts = None, nprint = 5):
        self.zmin = zmin
        self.zmax = zmax
        self.MCBeta = MCBeta
        self.Rate_Model = Rate_Model
        self.cheatType = cheatType
        self.cheatZ = cheatZ
        self.cheatCCSub = cheatCCSub
        self.cheatCCScale = cheatCCScale
        self.cuts = cuts
        self.nprint = nprint

        self.shiftFlagData = False
        self.shiftFlagSim = False

        self.globalChi2Storage = []
        self.globalNDataStorage = []
        self.globalZPhotBinStorage = []
        self.globalNDataIaPhotBinStorage = []
        self.globalNDataCCPhotBinStorage = []
        self.globalZTrueBinStorage = []
        self.globalNDataIaTrueBinStorage = []
        self.globalNDataCCTrueBinStorage = []
        print 'a'
        try: 
            self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
        except:
            self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95, skip_header = 5)
        print 'b'   
        self.simName = simName
        self.simgencat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
        print 'c'   
        try:
            SIMGEN = np.load(simgenfilename + '.npz')['a']
        except:
       
            SIMGEN = np.genfromtxt(simgenfilename, dtype=None, names = True, skip_footer=3, invalid_raise=False)
            np.savez_compressed(simgenfilename+'.npz', a = SIMGEN)

            print "WHY DO YOU HATE ME WHEN I SHOW YOU NOTHING BUT LOVE"
            print simgenfilename

        print 'd'
        SIMGEN = SIMGEN[SIMGEN['GENZ'] != 'GENZ']

        self.simgencat.params = {'flat':True, 'H0': simH0, 'Om0':simOmegaM, 'Ob0': simOb0, 'sigma8': simSigma8, 'ns': simNs}
        self.simgencat.cosmo = Cosmology.setCosmology('simCosmo', self.simcat.params)
        self.simgencat.OrigCatalog = np.copy(SIMGEN)
        self.simgencat.Catalog = np.copy(SIMGEN)
        self.simgencat.Catalog = self.simgencat.Catalog[self.simgencat.Catalog['GENZ'] != 'GENZ']
        self.simgencat.simname = simName
        self.simgencat.NSN = int(len(self.simgencat.Catalog['GENZ']))

        print "SIMGEN NUMBER"
        print self.simgencat.NSN
        print "SIMGENCAT FILE"
        print simfilename

        self.realName = realName
        try:
            self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95, skip_header = 5)
        except:
            self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)


        if self.cheatType:
            print "WARNING, THE FITTER IS CHEATING AND ELIMINATED NON-IAs USING SIM INFO"
            self.realcat.Catalog = self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1]
            self.simcat.Catalog = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1]


        for cut in cuts:
            self.realcat.Catalog = self.realcat.Catalog[(self.realcat.Catalog[cut[0]].astype(type(cut[1])) > cut[1]) & (self.realcat.Catalog[cut[0]].astype(type(cut[2])) < cut[2])]
            self.simcat.Catalog = self.simcat.Catalog[(self.simcat.Catalog[cut[0]].astype(type(cut[1])) > cut[1]) & (self.simcat.Catalog[cut[0]].astype(type(cut[2])) < cut[2])]

        
        
    def newData(self, realfilename, realName, simIndex = 100):
        self.realName = realName
        self.shitFlagData = False
        try:
            self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
        except:
            self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95, skip_header = 5 )
        if self.cheatType:
            print "WARNING, THE FITTER IS CHEATING AND ELIMINATED NON-IAs USING SIM INFO"
            self.realcat.Catalog = self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1]
        if simIndex < self.nprint:
            print 'N precuts'
            print self.realcat.Catalog['FITPROB'].shape
        for cut in cuts:
            self.realcat.Catalog = self.realcat.Catalog[(self.realcat.Catalog[cut[0]].astype(type(cut[1])) > cut[1]) & (self.realcat.Catalog[cut[0]].astype(type(cut[2])) < cut[2])]
            
        if simIndex < self.nprint:
            print "Minimum Fitprob"
            print np.min(self.realcat.Catalog['FITPROB'])
            print 'N postcuts'
            print self.realcat.Catalog['FITPROB'].shape

    def zSystematic(self):
        if self.shiftFlagData:
            print "DONT DOUBLE SHIFT"
            return 0
        if not self.shiftFlagSim:
        
            oldsimz = self.simcat.Catalog['zPHOT']
            oldsimtruez = self.simcat.Catalog['SIM_Z_CMB']
            stat, bins, binnum = stats.binnedstatistic(oldsimz, oldsimz - oldsimtruez, bins = self.bins, statistic = 'mean')
            self.zBiasShifts = stat
            newsimz = oldsimz - stat[binnum]
            assert(np.sum(np.abs(newsimz - oldsimz)) > 0)
            assert((oldzshape - np.arange(0, oldz.shape[0]).shape[0])< 1)
            self.shiftFlagSim = True
        oldz = self.realcat.Catalog['zPHOT']
        _,_, binnum = stats.binnedstatistic(oldz, oldz , bins = self.bins, statistic = 'mean')
        newz = oldz - self.zBiasShifts[binnum]
        oldzshape = oldz.shape[0]
        self.realcat.Catalog['zPHOT'].put(np.arange(0, oldz.shape[0]), newz)
        assert(np.sum(np.abs(newz - oldz)) > 0)
        assert((oldzshape - np.arange(0, oldz.shape[0]).shape[0])< 1)
        self.simFlagData = True
        
    def effCalc(self, fracContamCut = 0.0, nbins = 10, simIndex = 100):
        #### Do we want SNIas or all SN for efficiency?
        self.nbins = nbins
        self.typeString = ''

        if self.cheatZ:
            ztype = 'SIM_ZCMB'
        else:
            ztype = 'zPHOT'

        '''
        if (fracContamCut > 0.000000001) & (fracContamCut < 1.0):
            print " Cutting based on Frac Contam"
            histTot, binsX, binsY = np.histogram2d(self.simcat.Catalog[ztype], self.simcat.Catalog['MURES'], bins = nbins)
           
            histCC, binsX, binsY = np.histogram2d(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) != 1][ztype], self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) != 1]['MURES'], bins = (binsX, binsY))

            fracContam = histCC.astype(np.float)/histTot.astype(np.float)

            for fcRow, i in zip(fracContam, xrange(binsX.shape[0])):
                for fc, j in zip(fcRow, xrange(binsY.shape[0])):
                    if fc < fracContamCut:
                        continue
                    else:
                        simInBin = (self.simcat.Catalog[ztype] > binsX[i]) & (self.simcat.Catalog[ztype] < binsX[i+1]) & (self.simcat.Catalog['MURES'] > binsY[j]) & (self.simcat.Catalog['MURES'] < binsY[j+1])
                        realInBin = (self.realcat.Catalog[ztype] > binsX[i]) & (self.realcat.Catalog[ztype] < binsX[i+1]) & (self.realcat.Catalog['MURES'] > binsY[j]) & (self.realcat.Catalog['MURES'] < binsY[j+1])
                        self.simcat.Catalog = self.simcat.Catalog[np.invert(simInBin)]
                        self.realcat.Catalog = self.realcat.Catalog[np.invert(realInBin)]
        
        '''
        zPHOTs = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1][ztype].astype(float)

        zTRUEs = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1]['SIM_ZCMB'].astype(float)

        self.typeString = self.typeString + 'A1'
        
        
        binList = np.linspace(self.zmin, self.zmax, nbins+1)
        if simIndex < self.nprint:
            print "Type Location A"
            print "Choice A1"
            print zPHOTs.shape
            print zTRUEs.shape
            print binList
        
        self.binList = binList
        counts, zPhotEdges, zTrueEdges, binnumber = scipy.stats.binned_statistic_2d(zPHOTs, zTRUEs, zTRUEs, statistic = 'count', bins =  self.binList)
        assert(zPhotEdges.shape[0] == (self.nbins + 1))
        if simIndex < self.nprint:
            print "Type Location B"
            print "Choice B1"
        
        self.typeString = self.typeString + 'B1'
        zGenHist, zGenBins = np.histogram(self.simgencat.Catalog[self.simgencat.Catalog['GENTYPE'].astype(int) == 1]['GENZ'].astype(float), bins = self.binList)
        zSim1Hist, zSim1Bins = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) ==1]['SIM_ZCMB'].astype(float), bins = self.binList)
        
       
        
        if simIndex < self.nprint:
            print "counts of zTrue in each zPhot vs zTrue bin"
            print counts.astype(int)
            print "zGen Bins"
            print zGenBins
            print 'zGen Histogram'
            print zGenHist
            print "sum zGen events"
            print np.sum(zGenHist)
            print "sum zPhot events"
            print np.sum(counts)
            #print "DEBUG HERE"
            #assert(0)
        self.effmat = np.zeros((self.nbins,self.nbins))
        xMax = zPhotEdges.shape[0] - 2
        yMax = zTrueEdges.shape[0] - 2
        if simIndex < self.nprint:
            print zGenHist
            print counts.astype(int)

        for zPhotLedge, zPhotRedge, row, i in zip(zPhotEdges[:-1], zPhotEdges[1:], counts, xrange(xMax + 1)):
            zPhotCenter = (zPhotLedge + zPhotRedge)/2.0
            for zTrueLedge, zTrueRedge, count, j in zip(zTrueEdges[:-1], zTrueEdges[1:], row, xrange(yMax + 1)):
                zTrueCenter = (zTrueLedge + zTrueRedge)/2.0
                inCell = (zPHOTs > zPhotLedge) & (zPHOTs < zPhotRedge) & (zTRUEs > zTrueLedge)& (zTRUEs < zTrueRedge)
                zPhotCell = zPHOTs[inCell];zTrueCell = zTRUEs[inCell]
                self.effmat[i][j] = np.sum(inCell)
                assert(np.abs(np.sum(inCell) - count < 2))
                
        for row, i in zip(self.effmat, xrange(self.effmat.shape[0])):
            for j in xrange(row.shape[0]):
                self.effmat[i][j] /= zGenHist[j]

        if simIndex < self.nprint:
            print 'effmat'
            print self.effmat

        if simIndex == 0:
            extent = [zPhotEdges[0], zPhotEdges[-1], zTrueEdges[0], zTrueEdges[-1]]
            plt.figure()
            plt.imshow(np.flipud(counts), extent = extent, cmap = 'Blues')
            plt.colorbar()
            plt.savefig(self.realName + 'redshiftDistro.png')
            plt.clf()
            plt.close()
            plt.figure()
            plt.imshow(np.flipud(self.effmat), extent = extent, cmap = 'Blues', norm=mpl.colors.LogNorm())
            plt.colorbar()
            plt.savefig(self.realName + 'efficiencyMatrixLog.png')
            plt.clf()
            plt.close()
            plt.figure()
            plt.imshow(np.flipud(self.effmat), extent = extent, cmap = 'Blues')
            plt.colorbar()
            plt.savefig(self.realName + 'efficiencyMatrix.png')
            plt.clf()
            plt.close()
            
    def fit_rate(self, fixK = False, fixBeta = False, simIndex = 100, trueBeta = 0, CCScale = 1.0, TrueCCScale = 1.0, BetaInit = 0.0, kInit = 1.0, BetaErr = 1, kErr = 1, f_Js = None):
        #import iminuit as iM
        #from iminuit import Minuit as M
        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        if self.cheatZ:
            ztype = 'SIM_ZCMB'
        else:
            ztype = 'zPHOT'
        plt.switch_backend('Agg')

        if simIndex < self.nprint:
            print "Type Location C"
            print "Choice C1"

        if len(self.typeString) <= 4:
            self.typeString = self.typeString + 'C1'
        nSim, simBins = np.histogram(self.simgencat.Catalog[self.simgencat.Catalog['GENTYPE'].astype(int) == 1]['GENZ'].astype(float), bins=self.binList)
        nSim2, simBins2 = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) ==1][ztype].astype(float), bins=self.binList)
        
        
       
        nSim3, simBins3 = np.histogram(self.simcat.Catalog[ztype].astype(float), bins=self.binList)
        

        NCC , _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1][ztype].astype(float), bins=self.binList)

        OrigNCC = np.copy(NCC)
        if self.cheatCCSub:
            if self.cheatCCScale:
                print "WARNING: Only cheating on CC Subtraction not scale"
            print "Setting NCC to infinity to make sure that cheating correctly"
            print "Diagnostics after this point may be nonsense"
            print "NCC BeforeFck"
            print NCC
            NCC = NCC*1E100
            print "NCC AfterFck"
            print NCC    
        elif self.cheatCCScale:
            print "NCC Before1"
            print NCC
            print TrueCCScale
            NCC = applyCCScale(NCC, TrueCCScale)
            print "NCC After1"
            print NCC
        else: 
            print "NCC Before2"
            print NCC
            print CCScale
            NCC = applyCCScale(NCC, CCScale)
            print "NCC After2"
            print NCC
        #assert(0)

        
        NIa , _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1][ztype].astype(float), bins=self.binList)

        DebugNIaPhot, _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1]['zPHOT'].astype(float), bins=self.binList)
        DebugNCCPhot, _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1]['zPHOT'].astype(float), bins=self.binList)
        DebugNCCPhot = applyCCScale(DebugNCCPhot, CCScale)
        DebugNIaTrue, _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1]['SIM_ZCMB'].astype(float), bins=self.binList)
        DebugNCCTrue, _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1]['SIM_ZCMB'].astype(float), bins=self.binList)
        DebugNCCTrue = applyCCScale(DebugNCCTrue, CCScale)

        uselessCtr = 0
        for niap, nccp, niat, ncct, zb in zip(DebugNIaPhot, DebugNCCPhot, DebugNIaTrue, DebugNCCTrue,(self.binList[1:] + self.binList[:-1])/2.0 ):
            uselessCtr +=1
            self.globalZTrueBinStorage.append(zb)
            self.globalZPhotBinStorage.append(zb)
            self.globalNDataIaPhotBinStorage.append(niap)
            self.globalNDataCCPhotBinStorage.append(nccp)
            self.globalNDataIaTrueBinStorage.append(niat)
            self.globalNDataCCTrueBinStorage.append(ncct)
        print "UselessCtr"
        print uselessCtr
        


        try:
            TrueNCC, _ = np.histogram(self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'] !=1][ztype].astype(float), bins=self.binList)
            if simIndex < self.nprint:

                print "True NCC Data"
                print TrueNCC
        except:
            print "Using real data"

            TrueNCC = None

        nData, dataBins = np.histogram(self.realcat.Catalog[ztype].astype(float), bins=self.binList)
        if not(self.cheatCCSub):
            FracBad = NCC*1.0/(1.0*(NCC+NIa))
            nCCData = nData*FracBad
        else: 
            nCCData = TrueNCC*1.0
            FracBad = TrueNCC*1.0/nData

        print "PreScale NCC/nSim"
        print OrigNCC*1.0/(OrigNCC+NIa)
        
        print "PreScale Pred NCC Data"
        print OrigNCC*1.0/(OrigNCC+NIa)*nData

        print "PreScale Pred NCC Data if 2NCC"
        print OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData

        print "PreScale PredNCCData - TrueNCCData"
        print OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData - TrueNCC

        print "PreScale PredNCCData - TrueNCCData/ PredNCCData"
        print (OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData - TrueNCC)/(OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData)
        
        print "Mean of PreScale PredNCCData - TrueNCCData/ PredNCCData"
        print np.mean((OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData - TrueNCC)/(OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData))

        print "PostScale NCC/nData"
        print NCC*1.0/(NCC+NIa)

        if True or simIndex < self.nprint:
            print "Fraction of CCs in each bin"
            print FracBad

            print 'NCC'
            print NCC

            print 'nSim2'
            print nSim2
            print "nData, dataBins, realcat shape pre contam correction"
            print nData
            print dataBins
            print np.sum(self.realcat.Catalog[ztype].astype(float) > zmax)
            print np.sum(self.realcat.Catalog[ztype].astype(float) < zmin)
            print self.realcat.Catalog[ztype].shape
            
            print "Ratio nData/nSim"
            print 1.0*nData/(1.0*nSim)

            print "Ratio nData/nSim2"
            print 1.0*nData/(1.0*nSim2)
    
            print "Ratio nSim/nData"
            print 1.0*nSim2/(1.0*nData)

            print "Ratio nSim2/nData"
            print 1.0*nSim2/(1.0*nData)

            print "FracBad"
            print FracBad
            print 'NCCData'
            print nCCData

        if simIndex < self.nprint:

            print "overall Contam"
            print np.sum(NCC)*1.0/(np.sum(nSim3)*1.0)
        
        def chi2func(nData, nSim, effmat, fnorm, zCenters, k = 1.0, Beta = 0.0, zBreak = 1.0, dump = False, complexdump = False, modelError = False, nIA = None, nCC = None, Rate_Model = 'powerlaw', zbins = None, simIndex = 100, BetaPrior = (-3, 3), KPrior = (0.0, 50.0), TrueNCCData = None, f_1 = 1.0, f_2 = 1.0, f_3 = 1.0, f_4 = 1.0, f_5 = 1.0, f_6 = 1.0, f_7 = 1.0, f_8 = 1.0, f_9 = 1.0, f_10 = 1.0):
            Chi2Temp = 0.0
            if Rate_Model == 'powerlaw':
                f_Js = k*(1+zCenters)**Beta
            elif Rate_Model == 'discrete':
                f_Js = np.array([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10])
            elif Rate_Model == 'brokenpowerlaw':
                zCenters = (zbins[1:]+zbins[:-1])/2.0
                for zC in zCenters:
                    if zC < zBreak:
                        f_Js.append(k*(1+zC)**Beta)
                    elif not(temp is None):
                        f_Js.append(temp)
                    else:
                        temp = f_Js[-1]
                        f_Js.append(temp)
            else: 
                assert(0)
            print "f_Js init"
            print f_Js
            chi2Mat = np.zeros((self.nbins))
            adjNMC = np.zeros((self.nbins))
            if Rate_Model == 'discrete':
                kprior = 0
                betaprior = 0
            else:
                kprior  = weakPrior(k, KPrior)
                betaprior = weakPrior(Beta, BetaPrior)

            if dump and (self.nprint < simInd):
                print "kprior"
                print kprior
                print "betaprior"
                print betaprior
            if (nIA is None) or (nCC is None):
                print "No CC Cut"
                fracCCData = np.zeros(nData.shape)
            elif self.cheatCCSub:
                fracCCData = TrueNCC*1.0/nData 

            else:
                if Rate_Model == 'discrete':
                    print 'f_J adjusted CC Cut'
                    print Rate_Model
                    print nCC
                    print nIA
                    print np.array(f_Js)
                    fracCCData = (nCC*1.0)/((1.0*nCC + nIA*np.array(f_Js)))
                    print fracCCData
                else:
                    print "Beta Adjusted CC Cut"
                    print Rate_Model
                    #BetaRatio = k*(1+zCenters)**(Beta)#/(1+zCenters)**MCBeta
                    BetaRatio = (1+zCenters)**(Beta)#/(1+zCenters)**MCBeta
                    print "BadFracCCData"
                    print (nCC*1.0)/((1.0*nCC + nIA*BetaRatio))
                    print "bad NCCData"
                    print (nCC*1.0)/((1.0*nCC + nIA*BetaRatio))*nData
                    fracCCData = (nCC*1.0)/((1.0*nCC + nIA*BetaRatio))
            
            

            if dump and (self.nprint < simInd):
                print "fracCCData2"
                print fracCCData
                print "unscaled fracCCData"
                print (1.0*nCC)/(1.0*(nCC+nIA))
            if self.cheatCCSub:
                print "Cheating CC Sub"
                assert(not(TrueNCCData is None))
                nCCData = TrueNCCData
            else:
                print "Normal CC Sub"
                nCCData = nData*fracCCData
            if dump and (self.nprint < simInd):
                print "nCCData2"
                print nCCData
                if not(TrueNCCData is None):
                    print "TrueNCCData"
                    print TrueNCCData
                
            #print f_Js
            #Check if I am scaling errors down with increasing MC size. Make MC twice as large as "Data" to test.
            if dump: chi2Storage = []
            if dump: scaledNSimStor = []
            if dump:
                print "actually used NCC"
                #print nCC
                print nCCData
            for row, nDataI, nCCDataI, i in zip(effmat, nData, nCCData, xrange(self.nbins)):
                
                if dump and (self.nprint < simInd):
                    print 'effmat row'
                    print row
                    print 'nDataI'
                    print nDataI
                    print 'nCCDataI'
                    print nCCDataI
                    scaledNSimTemp = 0.0
                
                JSumTempNum = 0.0
                JSumTempDen = 0.0
                for eff, nSimJ, f_J, j in zip(row, nSim, f_Js, xrange(self.nbins)):
                    if dump and (i == j) and (self.nprint < simInd):
                        print 'NGen J'
                        print nSimJ
                        print 'JSumTempNum contr'
                        print nSimJ*f_J*eff*fnorm
                        print 'JSumTempDen contr'
                        print nSimJ*f_J*eff*fnorm*f_J*fnorm
                    if dump and (i != j) and self.cheatZ and (self.nprint < simInd):
                        if nSimJ*f_J*eff*fnorm > 0:
                            print " This should be zero but isnt "
                            print nSimJ*f_J*eff*fnorm
                            assert(0)
                    JSumTempNum += nSimJ*f_J*eff*fnorm
                    JSumTempDen += nSimJ*f_J*eff*fnorm*f_J*fnorm
                dataFunc = np.maximum(nDataI ,1)
                CCFunc = np.ceil(np.maximum(nCCDataI, 1))
                c2t = (nDataI - nCCDataI - JSumTempNum)**2/( dataFunc + CCFunc + JSumTempDen) + kprior + betaprior

                if dump and (self.nprint < simInd):
                    print i
                    print 'nDataI'
                    print nDataI
                    print 'fnCCDataI'
                    print nCCDataI
                    print 'fnorm'
                    print fnorm
                    print "JSumTempNum tot"
                    print JSumTempNum
                    print "JSumTempDen tot"
                    print JSumTempDen
                    print "Chi2Bin"

                    print c2t
                if dump:
                    chi2Storage.append(c2t)
                    if c2t > 5:
                        print 'INSANITY CHECK ABOVE'
                
                #    Chi2Temp += ((nDataI - nCCDataI - JSumTempNum)**2/(JSumTempNum + JSumTempDen))#*fnorm**2
                if nDataI > 1E-11 or JSumTempDen > 1E-11:
                    Chi2Temp += c2t
            if dump:
                return Chi2Temp, chi2Storage
            else:

                return Chi2Temp
        
        zCenters = (simBins[1:] + simBins[:-1])/2.0
 
        #Is this right? Everything else in the other side of the chi2 function should be Ia only
        if self.cheatCCSub:
            self.fracCCData = TrueNCC*1.0/nData
        else:
            self.fracCCData = (NCC*1.0)/(1.0*(NCC + NIa))
        fnorm = float(np.sum(nData*(1-self.fracCCData)))/float(np.sum(nSim))

        if self.Rate_Model == 'powerlaw':
            lamChi2 = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, nIA = NIa, nCC = NCC, simIndex = simIndex, TrueNCCData = TrueNCC)
            lamChi2Dump = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, dump = True, nIA = NIa, nCC = NCC, simIndex = simIndex, TrueNCCData = TrueNCC)
            MinObj = M(lamChi2, k = kInit, error_k = kErr , Beta = BetaInit, error_Beta = BetaErr, limit_k = (0.0, None), fix_k = fixK, fix_Beta = fixBeta)
            c2i, _ = lamChi2Dump(1.0, 0.0)

            print "Chi2 init = {0}".format(round(c2i, 4))
        elif self.Rate_Model == 'brokenpowerlaw':
            lamChi2 = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, 1.0, nCCData = NCCData, simIndex = simIndex, TrueNCCData = TrueNCC, Rate_Model = 'brokenpowerlaw')
            lamChi2Dump = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, 1.0, dump = True, nCCData = nCCData, simIndex = simIndex, TrueNCCData = TrueNCC, Rate_Model = 'brokenpowerlaw')

            MinObj = M(lamChi2, k = kInit, error_k = kErr , Beta = BetaInit, error_Beta = BetaErr, limit_k = (0.0, None), fix_k = fixK, fix_Beta = fixBeta)
            c2i, _ = lamChi2Dump(1.0, 0.0)

            print "Chi2 init = {0}".format(round(c2i, 4))
        elif self.Rate_Model == 'discrete':
            
            lamChi2 = lambda f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10: chi2func(nData, nSim, self.effmat, fnorm, zCenters, 1.0, nIA = NIa, nCC = NCC, simIndex = simIndex, TrueNCCData = TrueNCC, f_1 = f_1, f_2 = f_2,f_3 = f_3, f_4 = f_4,f_5 = f_5, f_6 = f_6,f_7 = f_7, f_8 = f_8,f_9 = f_9, f_10 = f_10, Rate_Model = 'discrete' )
            lamChi2Dump = lambda f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10: chi2func(nData, nSim, self.effmat, fnorm, zCenters, 1.0, nIA = NIa, nCC = NCC, simIndex = simIndex, TrueNCCData = TrueNCC, f_1 = f_1, f_2 = f_2,f_3 = f_3, f_4 = f_4,f_5 = f_5, f_6 = f_6,f_7 = f_7, f_8 = f_8,f_9 = f_9, f_10 = f_10, dump = True, Rate_Model = 'discrete')

            c2i, _ = lamChi2Dump(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

            print "Chi2 init = {0}".format(round(c2i, 4))

            MinObj = M(lamChi2, f_1 = 1.0, error_f_1 = 1.0, limit_f_1 = (0.0, None), f_2 = 1.0, error_f_2 = 1.0, limit_f_2 = (0.0, None), f_3 = 1.0, error_f_3 = 1.0, limit_f_3 = (0.0, None), f_4 = 1.0, error_f_4 = 1.0, limit_f_4 = (0.0, None), f_5 = 1.0, error_f_5 = 1.0, limit_f_5 = (0.0, None), f_6 = 1.0, error_f_6 = 1.0, limit_f_6 = (0.0, None), f_7 = 1.0, error_f_7 = 1.0, limit_f_7 = (0.0, None), f_8 = 1.0, error_f_8 = 1.0, limit_f_8 = (0.0, None), f_9 = 1.0, error_f_9 = 1.0, limit_f_9 = (0.0, None), f_10 = 1.0, error_f_10 = 1.0, limit_f_10 = (0.0, None))

        if self.Rate_Model == 'discrete':
            c2f, c2stor = lamChi2Dump(MinObj.values['f_1'],MinObj.values['f_2'],MinObj.values['f_3'],MinObj.values['f_4'],MinObj.values['f_5'],MinObj.values['f_6'],MinObj.values['f_7'],MinObj.values['f_8'],MinObj.values['f_9'],MinObj.values['f_10'])
        else: 

            c2f, c2stor = lamChi2Dump(MinObj.values['k'], MinObj.values['Beta'])

        

        
        #MinObj = M(lamChi2, k = 1.0, fix_k = True, Beta = 0.0, error_Beta = 0.1)
        

        MinObj.set_strategy(2)
        fmin, param = MinObj.migrad(nsplit= 10)

        

        plt.scatter(nData, c2stor)
        plt.xlabel('nData')
        plt.ylabel('chi2 in bin')
        plt.savefig(self.realName + 'Chi2VsnData.png')
        plt.clf()
        print "Shapes of things"
        print len(c2stor)
        print nData.shape

        print dataBins.shape

        print self.binList.shape
        print DebugNIaPhot.shape
        print DebugNCCPhot.shape
        print DebugNIaTrue.shape
        print DebugNCCTrue.shape

        for c2, nd in zip(c2stor, nData):
            self.globalChi2Storage.append(c2)
            self.globalNDataStorage.append(nd)

        if self.Rate_Model == 'discrete':
            fJList = [MinObj.values['f_1'],MinObj.values['f_2'],MinObj.values['f_3'],MinObj.values['f_4'],MinObj.values['f_5'],MinObj.values['f_6'],MinObj.values['f_7'],MinObj.values['f_8'],MinObj.values['f_9'],MinObj.values['f_10']]
            fJErrList = [MinObj.errors['f_1'],MinObj.errors['f_2'],MinObj.errors['f_3'],MinObj.errors['f_4'],MinObj.errors['f_5'],MinObj.errors['f_6'],MinObj.errors['f_7'],MinObj.errors['f_8'],MinObj.errors['f_9'],MinObj.errors['f_10']]

            
            self.fJList = fJList
            self.fJErrList = fJErrList
            self.Beta = None
            self.k = None
            self.kErr = None
            self.BetaErr = None
            print fJList
            print fJErrList
        else:
            k = MinObj.values['k']
            kErr = MinObj.errors['k']
            Beta = MinObj.values['Beta']
            BetaErr = MinObj.errors['Beta']
            self.k = k
            self.Beta = Beta
            self.kErr = kErr
            self.BetaErr = BetaErr
            #/(self.nbins - 2)
            self.BetaRatio = (1+zCenters)**(Beta)
            self.fJList = None
            self.fracCCData = (NCC*1.0)/(1.0*(1.0*NCC + NIa*self.BetaRatio))
            print "Chi2 final = {0}".format(round(lamChi2Dump(self.k, self.Beta)[0], 4))
        self.chi2 = fmin.fval
        print "Chi2final? = {0}".format(round(fmin.fval, 4))
        

        #fJs = np.ones(zCenters.shape)
        '''
        xgrid,ygrid, sigma, rawdata = MinObj.mncontour_grid('k', 'Beta', numpoints=400, sigma_res = 4, nsigma = 2.0)
        plt.figure()
        plt.clf()
        CS = plt.contour(xgrid, ygrid + self.MCBeta, sigma, levels = [0.5, 1.0, 1.5, 2.0])
        plt.clabel(CS, fontsize=7, inline=1)
        plt.xlabel('k')
        plt.ylabel('Beta')
        plt.savefig('{0}_{1}_k_beta_contour.png'.format(self.realName, self.simName))
        plt.close()
        '''

    
        #plt.axhline(y = self.MCBeta, c = 'k', label = 'True Beta')
        #plt.axhline(y = Beta + self.MCBeta, c = 'g', label= 'Best Fit Beta')
        #plt.axvline(x = k, label = 'Best Fit k')
        

    def chi2V2(self, fJs, fJErrs, zCenters, k, Beta):
            fitfJs = k*(1+zCenters)**Beta
            Chi2Temp = 0
            for fJ, fitfJ, fJErr in zip(fJs, fitfJs, fJErrs):
                Chi2Temp += (fJ - fitfJ)**2/(fJ + fJErr)
            return Chi2Temp

def weakPrior(value, priorTuple):
    if value < priorTuple[1]:
        if value > priorTuple[0]:
            return 1
        else: 
            return (value - priorTuple[0])**4
    else:
        return (value - priorTuple[1])**4



def getCCScale(simCat, dataCat, MURESWindow = (-1, 1), zbins = [0.0, 0.3, 0.6, 0.9, 1.2], muresBins = np.linspace(-2, -1, 4), Beta = None, binList = None, fracCCData = None, outfilePrefix = 'Test', Rate_Model = 'powerlaw', f_Js = None):
    import iminuit as iM
    from iminuit import Minuit as M
    CCScales = []
    CCScaleErrs = []
    if not(f_Js is None):
        f_Js = np.array(f_Js)

    tempSimCC = simCat[simCat['SIM_TYPE_INDEX'].astype(int) != 1]
    tempSimIa = simCat[simCat['SIM_TYPE_INDEX'].astype(int) == 1]

    allSimCC = np.copy(tempSimCC)
    allSimIa = np.copy(tempSimIa)
    allData = np.copy(dataCat)

    simHistCC, simBinsCC  = np.histogram(tempSimCC['MURES'], bins = muresBins)
    simHistIa, simBinsIa  = np.histogram(tempSimIa['MURES'], bins = muresBins)

    simZHistCC, simZBinsCC = np.histogram(tempSimCC['zPHOT'], bins = binList)
    simZHistIa, simZBinsIa = np.histogram(tempSimIa['zPHOT'], bins = binList)

    binCent = (simZBinsIa[1:] + simZBinsIa[:-1])/2.0


    print "components of simHist"
    print simZHistIa
    print binCent
    print simHistCC

    if Rate_Model == 'discrete':
        simHist = simZHistIa*np.array(f_Js) + simZHistCC
    else:
        simHist = (simZHistIa*(1+binCent)**Beta) + simZHistCC
    print 'num and denom of fnorm2'
    print float(dataCat.shape[0])
    print float(np.sum(simHist))

    fnorm2 = float(dataCat.shape[0])/float(np.sum(simHist))
    print "precut"
    print simCat.shape
    print dataCat.shape
    simCat = simCat[(simCat['MURES'] < MURESWindow[0]) | (simCat['MURES'] > MURESWindow[1]) ]
    dataCat = dataCat[(dataCat['MURES'] < MURESWindow[0]) | (dataCat['MURES'] > MURESWindow[1]) ]
    


    print "postcut"
    print simCat.shape
    print dataCat.shape

  
    print 'post MURES Cut Ndata and Nsim'
    print dataCat.shape
    print simCat.shape

    for zl, zh in zip(zbins[:-1], zbins[1:]):

        tempSim = simCat[(simCat['zPHOT'] < zh) & (simCat['zPHOT'] > zl)]
        tempData = dataCat[(dataCat['zPHOT'] < zh) & (dataCat['zPHOT'] > zl)]


        allSimCCZbin = allSimCC[(allSimCC['zPHOT'] < zh) & (allSimCC['zPHOT'] > zl)]
        allSimIaZbin = allSimIa[(allSimIa['zPHOT'] < zh) & (allSimIa['zPHOT'] > zl)]

        print "all Sim CC Zbin/IaZbin"
        print allSimCCZbin.shape[0]
        print allSimIaZbin.shape[0]

        allDataZbin = allData[(allData['zPHOT'] < zh) & (allData['zPHOT'] > zl)]


        #binchoice = np.abs(binCent - np.mean(tempData['zPHOT']))

        #fracCCCent = fracCCData[np.argmin(binchoice)]
        
        if type(muresBins == int):
            histD, muresBins = np.histogram(tempData['MURES'], bins = muresBins)
            binsD = muresBins
            #histDAll, muresBins = np.histogram(allDataZbin['MURES'], bins = binsD)
        else: 
            histD, binsD = np.histogram(tempData['MURES'], bins = muresBins)
            #histDAll = np.histogram(allDataZbin['MURES'], bins = muresBins)




        histS, binsS  = np.histogram(tempSim['MURES'], bins = muresBins)
        tempSimCC = tempSim[tempSim['SIM_TYPE_INDEX'] != 1]
        tempSimIa = tempSim[tempSim['SIM_TYPE_INDEX'] == 1]
        histSCC, binsSCC  = np.histogram(tempSimCC['MURES'], bins = muresBins)
        if Rate_Model == 'discrete':
            histSIa, binsSIa  = np.histogram(tempSimIa['MURES'], bins = muresBins)
            histSIa = histSIa*f_Js
        else:
            histSIa, binsSIa  = np.histogram(tempSimIa['MURES'], bins = muresBins, weights = (1+tempSimIa['zPHOT'])**Beta)
        #histSAllCC, binsSAllCC  = np.histogram(allSimCCZbin['MURES'], bins = muresBins)
        #histSAllIa, binsSAllIa  = np.histogram(allSimIaZbin['MURES'], bins = muresBins)

        


        #histSIa = histSIa*(1+np.mean(tempSim['zPHOT'])**Beta)
        
        print " MURES tail"
        print histS
        print histSCC
        print histSIa

        print "MURES bins"
        print binsS
        print binsSCC
        print binsSIa

        print "MURES tail scaled by fnorm2"
        print histS*fnorm2
        print histSCC*fnorm2
        print histSIa*fnorm2

        print "Data tail"
        print histD

        print "fnorm2"
        print fnorm2

        

        print "data events outliers only and total"
        print histD
        print allDataZbin.shape[0]

        #R = np.array(histD).astype(float)/np.array(histDAll).astype(float)
        R = np.array(histD).astype(float)/float(allDataZbin.shape[0])

        print "R"

        print R

        print "Hist CC, outlier and total"

        print histSCC
        #print histSAllCC
        print allSimCCZbin.shape[0]

        print "pre Beta Correction allSimIa"
        print allSimIaZbin.shape[0]

        #print "Beta Correction Factor"
        #print (1+np.mean(allSimIaZbin['zPHOT']))**Beta
        #print "BetaCorrected allSimIa"
        #print allSimIaZbin.shape[0]*(1+np.mean(allSimIaZbin['zPHOT']))**Beta
        #betaCorrAllSimIaZbin = allSimIaZbin.shape[0]*(1+np.mean(allSimIaZbin['zPHOT']))**Beta
        if Rate_Model == 'discrete':
            hist, bins = np.histogram(allSimIaZbin['zPHOT'], bins = 10)
            print 'fJ shape'
            print f_Js.shape
            print f_Js
            print hist
            print bins
            betaCorrAllSimIaZbin =np.sum(hist*f_Js)
        else:
            betaCorrAllSimIaZbin =np.sum((1+ allSimIaZbin['zPHOT'])**Beta)
        #S = float(np.array(R*histSAllIa) - np.array(histSIa))/float(np.array(histSCC) - np.array(R*histSAllCC))

        try:
            print "Test S"
            print R
            print betaCorrAllSimIaZbin
            print histSIa
            print histSCC
            print allSimCCZbin.shape
            print 'EEE'
            print np.array(R*betaCorrAllSimIaZbin)
            print 'DDD'
            print np.array(histSIa)
            print 'CCC'
            print (np.array(histSCC) - np.array(R*allSimCCZbin.shape[0]))
            print "AAA"
            print (np.array(R*betaCorrAllSimIaZbin) - np.array(histSIa))/(np.array(histSCC) - np.array(R*allSimCCZbin.shape[0]))
            print "BBB"
            S = (np.array(R*betaCorrAllSimIaZbin) - np.array(histSIa))/(np.array(histSCC) - np.array(R*allSimCCZbin.shape[0]))
        except: 
            S = np.nan

        print "S WTF"
        print S


        print "Uncertainty Related Bullshit"
        '''
        print "Delta R"

        dR = np.sqrt(histD + histDAll)

        print dR

        num1 = np.sqrt(np.sqrt((dR/R)**2 + histSAllIa) + histSIa)

        num2 = np.sqrt(np.sqrt((dR/R)**2 + histSAllCC) + histSCC)

        den1 = (R*histSAllIa - histSIa)

        den2 = (histSCC - R*histSAllCC)


        dS = np.sqrt((num1/den1)**2 + (num2/den2)**2)
        '''
        #ddnCC = np.sqrt(histSCC)*(histSIa - histSAllIa*R)/(histSCC - R*histSAllCC)**2

        #ddNCC = np.sqrt(histSAllCC)*R*(histSAllIa*R - histSIa)/(histSCC - R*histSAllCC)**2

        #ddnIa = np.sqrt(histSIa)/(histSCC - R*histSAllCC)
        #ddNIa = np.sqrt(histSAllIa)*R/(histSCC - R*histSAllCC)

        ddnCC = np.sqrt(histSCC)*(histSIa - allSimIaZbin.shape[0]*R)/(histSCC - R*allSimCCZbin.shape[0])**2

        ddNCC = np.sqrt(allSimCCZbin.shape[0])*R*(allSimIaZbin.shape[0]*R - histSIa)/(histSCC - R*allSimCCZbin.shape[0])**2

        ddnIa = np.sqrt(histSIa)/(histSCC - R*allSimCCZbin.shape[0])
        ddNIa = np.sqrt(allSimIaZbin.shape[0])*R/(histSCC - R*allSimCCZbin.shape[0])

        #ddR = (histSAllIa*histSCC - histSAllCC * histSIa)/(histSCC - R*histSAllCC)**2

        dS = np.sqrt(ddnCC**2 + ddNCC**2 + ddnIa**2 + ddNIa**2)# + ddR**2)

        print "ddnCC"

        print ddnCC

        print "ddNCC"

        print ddNCC

        print "ddnIa"

        print ddnIa

        print "ddNIa"

        print ddNIa

        #print "ddR"

        #print ddR

        print "Delta S"

        print dS

        #assert(S > 0)
        if S < 0: 
            S = np.nan

        CCScales.append(S)
        CCScaleErrs.append(dS[0])



   

        plt.step((simBinsCC[1:] + simBinsCC[:-1])/2.0, simHistCC, c = 'b', where = 'mid',  label = 'prescaled Sim CC')
        plt.step((simBinsCC[1:] + simBinsCC[:-1])/2.0, CCScales[0]*simHistCC*fnorm2, c = 'g', where = 'post',  label = 'postscaledSimCC')
        plt.step((muresBins[1:] + muresBins[:-1])/2.0, histD, c = 'r', where = 'mid',  label = 'data')
        plt.savefig(outfilePrefix + 'ScaledHist.png')
        plt.clf()
    return CCScales, CCScaleErrs


def applyCCScale(NCC, CCScales, datazbins = None,  zbins = None):

    NCCScaled = CCScales*NCC


    return NCCScaled

if __name__ == '__main__':
    from sys import argv
    print "argv"
    print argv
    datadir = argv[1]
    simdir = argv[2]
    dataname = argv[3]
    print "dataname"
    simname = argv[4]
    print simname
    simgenfile = argv[5]
    print simgenfile
    NNCut = False
    cheatType = bool(int(argv[6]))
    cheatZ = bool(int(argv[7]))
    trueBeta = float(argv[8])
    paramFile = argv[9]
    cutFiles = argv[10:]
    
    if( ('Combine' in simdir) or ('SALT2' in simdir)) &  (('Combine' in datadir) or ('SALT2' in simdir)):
        NNCut = True
        NNProbCut = 0.95
    
    #if len(argv) > 6:
    #    NNCut = True
    #    NNProbCut = 0.9
    #    NNData = argv[6]
    #    NNSim = argv[7]

    
    #default params

    zmin = 0.1
    zmax = 1.2
    MJDMin = 0.0
    MJDMax = np.inf
    bins = "equalSize" 
    runFit = True
    fracContamCuts = [-1]
    fixBeta = True
    fixK = False
    nbins = 10
    ScaleMuResCutLow = -1
    ScaleMuResCutHigh = 1
    muresBins = 1
    muresBinsLow  = 3
    muresBinsHigh = 3
    scaleZBins = [0.0, 1.2]
    cheatCCSub = False
    cheatCCScale = False
    ZSysFlag = False

    #override file

    params = open(paramFile, 'r').readlines()

    for p in params:
        exec(p)

    kmean = []
    ksigma = []
    kErr = []
    BetaMean = []
    BetaSigma= []
    BetaErr = []
    Chi2Mean = []
    Chi2Sigma = []

    CCScaleStorageGlobal = []
    CCScaleErrStorageGlobal = []


    #MURES_Cuts = [2.0]
    #MURES_Cuts = [1.0, 1.5, 2.0, 3.0, 4.0, 99.0, 2.0]
    #for MURES_Cut in MURES_Cuts:
    fcc = -1
    for cf in cutFiles:
        cuts = [] #    cuts = [('FITPROB', 0.01, np.inf), ('NN_PROB_IA', NNProbCut, np.inf)]

        cutlist = open(cf, 'r').readlines()
        for l in cutlist:
            spl = l.split()
            cuts.append(('{0}'.format(spl[0]), float('{0}'.format(spl[1])), float('{0}'.format(spl[2]))))

        ks = []
        kErrs = []
        Betas = []
        BetaErrs = []
        Chi2s = []

        CCScaleStorage = []
        CCScaleErrStorage = []


        nFail = 0
        simLoaded = False
        #print "FUCK MPI"
        #if Rate_Model == 'discrete':
        #    subprocess.call(['python', 'constructChi2Func.py', str(nbins)], shell = False)
        #print "MPI Fucked"
        if '{' in datadir:
            #nfile = 100
            nfile = 49
        else:
            nfile = 2
        for simInd in xrange(1,nfile):
            

            #print "Sim {0}".format(simInd)
            #SimBeta = 2.1 # simdir.split('_')[-3]
            #SimR0 = 1.7*10**-5 #simdir.split('_')[-5]
            #print "Sim R0 = {1}; Sim Beta = {0}".format(SimBeta, SimR0)

            
            print datadir.format(simInd)
            if simLoaded:
                try:
                    if NNCut:
                        RateTest.newData(datadir.format(simInd), dataname.format(simInd), simIndex = simInd)
                        if ZSysFlag:
                            RateTest.zSystematic()
                    else:
                        RateTest.newData(datadir.format(simInd), dataname.format(simInd), simIndex = simInd)
                        if ZSysFlag:
                            RateTest.zSystematic()


                   

                    #RateTest.effCalc(nbins = nbins, fracContamCut = fcc, simIndex = simInd)
                    #RateTest.effCalc(nbins = 20)
                    BetaIter = []
                    BetaErrIter = []
                    CCIter = []
                    CCErrIter = []
                    RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, simIndex = simInd, trueBeta = trueBeta - 2.11, CCScale = 1.0, TrueCCScale = TrueCCScale)
                    BetaIter.append(RateTest.Beta)
                    BetaErrIter.append(RateTest.BetaErr)
                    for iteration in xrange(nIter):
                        CCScale, CCScaleErr =  getCCScale(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBins = muresBins, Beta = RateTest.Beta, binList = RateTest.binList, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname,Rate_Model = Rate_Model, f_Js =RateTest.fJList)
                        CCIter.append(CCScale[0])
                        CCErrIter.append(CCScaleErr[0])
                        RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - 2.11, CCScale = CCScale, TrueCCScale = TrueCCScale, BetaInit = RateTest.Beta, kInit = RateTest.k, BetaErr = RateTest.BetaErr, kErr = RateTest.kErr, f_Js =RateTest.fJList)
                        BetaIter.append(RateTest.Beta)
                        BetaErrIter.append(RateTest.BetaErr)

                    CCScale, CCScaleErr =  getCCScale(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBins = muresBins, Beta = RateTest.Beta, binList = RateTest.binList, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname,Rate_Model = Rate_Model, f_Js =RateTest.fJList)
                    CCIter.append(CCScale[0])
                    CCErrIter.append(CCScaleErr[0])
                    print "Beta Progression"
                    print BetaIter
                    print "Beta Err Progressions"
                    print BetaErrIter
                    print "CCScale Progression"
                    print CCIter
                    print "CCScale Err Progression"
                    print CCErrIter
                    print "Mean Betas"
                    print np.mean(BetaIter)
                    print "Mean CCScales"
                    print np.mean(CCIter)
                    

                    print "AAA CC Scales"
                    
                    CCScale, CCScaleErr =  getCCScale(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBins = muresBins, Beta = RateTest.Beta, binList = RateTest.binList, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname, Rate_Model = Rate_Model, f_Js =RateTest.fJList)
                    print CCScale
                    CCScaleStorage.append(CCScale[0])
                    CCScaleErrStorage.append(CCScaleErr[0])
                    


                    ks.append(RateTest.k)
                    kErrs.append(RateTest.kErr)
                    Betas.append(RateTest.Beta)
                    BetaErrs.append(RateTest.BetaErr)
                    Chi2s.append(RateTest.chi2)
                    print "CCScale Storage Iter {0}".format(simInd)
                    print CCScaleStorage
                    print CCScale
                    print CCScale[0]

                    print "CCScaleErr Storage Iter {0}".format(simInd)
                    print CCScaleErrStorage
                    print CCScaleErr
                    print CCScaleErr[0]

                    dnamestr = datadir.format(simInd)

                    cutdnamestr = dnamestr.split('.')[0] + '+CUTS.FITRES.gz'
                    if saveCuts:
                        np.savetxt(cutdnamestr, RateTest.realcat.Catalog, delimiter = ' ', fmt='%s')

                    #with open(cutdnamestr, 'rb') as f_in:
                    #    with gzip.open(cutdnamestr + '.gz', 'wb') as f_out:
                    #        shutil.copyfileobj(f_in, f_out)
                except Exception, e:
                    print "FAILURE"
                    print e
                    traceback.print_exc()
                    nFail +=1
            else:
                try:

                    RateTest = Rate_Fitter(datadir.format(simInd), dataname.format(simInd), simdir, simname,simgenfile, 2.1, zmin = 0.1, zmax = 1.2, cheatZ = cheatZ, cheatType = cheatType, cuts = cuts, cheatCCSub = cheatCCSub, cheatCCScale = cheatCCScale, Rate_Model = Rate_Model)# , MJDMin = 0, MJDMax = np.inf)
                    
                    if ZSysFlag:
                            RateTest.zSystematic()
                    simLoaded = True
                    RateTest.effCalc(nbins = nbins, fracContamCut = fcc)
                    #RateTest.effCalc(nbins = 20)
                    BetaIter = []
                    CCIter = []
                    RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, simIndex = simInd, trueBeta = trueBeta - 2.11, CCScale = 1.0, TrueCCScale = TrueCCScale)
                    BetaIter.append(RateTest.Beta)
                    BetaErrIter.append(RateTest.BetaErr)
                    for iteration in xrange(nIter):
                        CCScale, CCScaleErr =  getCCScale(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBins = muresBins, Beta = RateTest.Beta, binList = RateTest.binList, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname, Rate_Model = Rate_Model, f_Js =RateTest.fJList)
                        CCIter.append(CCScale[0])
                        CCErrIter.append(CCScaleErr[0])

                        RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - 2.11, CCScale = CCScale, TrueCCScale = TrueCCScale, BetaInit = RateTest.Beta, kInit = RateTest.k, BetaErr = RateTest.BetaErr, kErr = RateTest.kErr)
                        BetaIter.append(RateTest.Beta)
                        BetaErrIter.append(RateTest.BetaErr)

                    CCScale, CCScaleErr =  getCCScale(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBins = muresBins, Beta = RateTest.Beta, binList = RateTest.binList, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname, Rate_Model = Rate_Model, f_Js =RateTest.fJList)
                    CCIter.append(CCScale[0])
                    CCErrIter.append(CCScaleErr[0])
                    print "Beta Progression"
                    print BetaIter
                    print "Beta Err Progressions"
                    print BetaErrIter
                    print "CCScale Progression"
                    print CCIter
                    print "CCScale Err Progression"
                    print CCErrIter
                    print "Mean Betas"
                    print np.mean(BetaIter)
                    print "Mean CCScales"
                    print np.mean(CCIter)
                    
                    
                    print "AAA CC Scales"
                    
                    CCScale, CCScaleErr =  getCCScale(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBins = muresBins, Beta = RateTest.Beta, binList = RateTest.binList, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname, f_Js =RateTest.fJList)
                    print 'CC Scale'
                    print CCScale
                    CCScaleStorage.append(CCScale[0])
                    CCScaleErrStorage.append(CCScaleErr[0])

                    dnamestr = datadir.format(simInd)

                    cutdnamestr = dnamestr.split('.')[0] + '+CUTS.FITRES.gz'

                    np.savetxt(cutdnamestr, RateTest.realcat.Catalog, delimiter = ' ', fmt='%s')

                    #with open(cutdnamestr, 'rb') as f_in:
                    #    with gzip.open(cutdnamestr + '.gz', 'wb') as f_out:
                    #        shutil.copyfileobj(f_in, f_out)



                    cutsnamestr = simname.split('.')[0] + '+CUTS.FITRES.gz'

                    np.savetxt(cutsnamestr, RateTest.realcat.Catalog, delimiter = ' ', fmt = '%s')

                    #with open(cutsnamestr, 'rb') as f_in:
                    #    with gzip.open(cutsnamestr + '.gz', 'wb') as f_out:
                    #        shutil.copyfileobj(f_in, f_out)


                    ks.append(RateTest.k)
                    kErrs.append(RateTest.kErr)
                    Betas.append(RateTest.Beta)
                    BetaErrs.append(RateTest.BetaErr)
                    Chi2s.append(RateTest.chi2)
                    print "CCScale Storage Iter {0}".format(simInd)
                    print CCScaleStorage
                    print CCScale
                    print CCScale[0]



                except Exception, e:
                    print "FAILURE"
                    print e
                    traceback.print_exc()
                    nFail +=1

        print "Number of Failures"
        print nFail
        print "mean k"
        print np.mean(ks)
        print "mean kerrs"
        print np.mean(kErrs)
        print "std. k"
        print np.std(ks)
        print "Mean beta"
        print np.mean(Betas)
        print "Mean betaerrs"
        print np.mean(BetaErrs)
        print "std. beta"
        print np.std(Betas)
        kmean.append(np.mean(ks))
        ksigma.append(np.std(ks))
        kErr.append(np.mean(kErrs))
        BetaMean.append(np.mean(Betas))
        BetaSigma.append(np.std(Betas))
        BetaErr.append(np.mean(BetaErrs))
        Chi2Mean.append(np.mean(Chi2s))
        Chi2Sigma.append(np.std(Chi2s))

        bins0 = [0.0, 1.0, 3.0, 6.0, 9.0, 11.0, 14.0, 18.0, 25.0]

        hist, bins = np.histogram(Chi2s, bins = bins0)
        xs = (bins[1:] + bins[:-1])/2.0

        plt.bar(xs, hist, width = bins[1:] - bins[:-1])

        chi2s = scipy.stats.chi2.pdf(xs, 11)

        norm = np.max(hist)*1.0/np.max(chi2s)

        plt.plot(xs, scipy.stats.chi2.pdf(xs, 11)*norm, color = 'g')
        if cheatType and not cheatZ:
            plt.savefig(dataname +'Chi2Plot_CheatType.png')
        elif cheatZ and not cheatType:
            plt.savefig(dataname  +'Chi2Plot_CheatZ.png')
        elif cheatZ and cheatType:
            plt.savefig(dataname +'Chi2Plot_CheatTypeZ.png')
        else:
            plt.savefig(dataname +'Chi2Plot.png')

        
        print "AAA CC Scale means (weighted, unweighted)"
        print np.average(ma.masked_invalid(np.array(CCScaleStorage)),weights = 1.0/ma.masked_invalid(CCScaleErrStorage)**2, axis = 0)
        print np.mean(ma.masked_invalid(np.array(CCScaleStorage)), axis = 0)
        print "AAA CC Scale stds"
        print np.nanstd(np.array(CCScaleStorage), axis = 0)
        CCScaleStorageGlobal.append(CCScaleStorage)

        



if cheatType:
    print "THESE RESULTS ONLY INCLUDE TRUE Ias BECAUSE WE CHEATED AND USED THE SIM INFORMATION"
if cheatZ:
    print "THESE RESULTS Use Simulated Redshift info"
print "lengths of lists"

print len(RateTest.globalNDataStorage)
print len(RateTest.globalChi2Storage)
print len(RateTest.globalZPhotBinStorage)
print len(RateTest.globalNDataIaPhotBinStorage)
plt.clf()
plt.scatter(RateTest.globalNDataStorage, RateTest.globalChi2Storage)
plt.xlabel('nData')
plt.ylabel('chi2 in bin')
string = ''
if cheatType: string += 'CheatType'
if cheatZ: string += 'CheatZ'
print 'string here'
print string
plt.savefig(RateTest.realName + 'Chi2VsnData' + string +'.png')
plt.clf()


plt.scatter(RateTest.globalZPhotBinStorage, RateTest.globalChi2Storage)
plt.xlabel('zPhot bin center')
plt.ylabel('chi2 in bin')
plt.savefig(RateTest.realName + 'Chi2VsZPhot' + string +'.png')
plt.clf()

plt.clf()
plt.scatter(RateTest.globalZPhotBinStorage, RateTest.globalNDataIaPhotBinStorage, s = 1, c = 'r', label = 'Type Ia Data, zPhot')
plt.scatter(RateTest.globalZPhotBinStorage, RateTest.globalNDataCCPhotBinStorage, s = 1, c = 'b', label = 'CC Data, zPhot')
plt.scatter(RateTest.globalZTrueBinStorage, RateTest.globalNDataIaTrueBinStorage, s = 1, c = 'Pink', label = 'Type Ia Data, zTrue')
plt.scatter(RateTest.globalZTrueBinStorage, RateTest.globalNDataCCTrueBinStorage, s = 1, c = 'Cyan', label = 'CC  Data, zTrue')
plt.yscale('log')
plt.xlabel('redshift either true or phot')
plt.legend()
plt.savefig(RateTest.realName + 'AggregateZDistro' + string +'.png')


#print "MURES CUTS"
#print MURES_Cuts
print "Frac Contam Cuts"
print fracContamCuts

print "Kmeans"
print kmean
print "Ksigmas"
print ksigma
print "BetaMeans"
print BetaMean
print "BetaSigmas"
print BetaSigma
print "Chi2Means"
print Chi2Mean
print "Chi2Sigma"
print Chi2Sigma

assert(fracContamCuts[0] == -1)
outfile = dataname

print "outfile Pre Prefix"
print outfile

if cheatType:
    outfile = outfile + '_CheatType'
    if cheatZ:
        outfile = outfile + 'Z'
elif cheatZ:
    outfile = outfile + '_CheatZ'

outfile = outfile + '.txt'
print "Outfile Name"
if not(os.path.isfile(outfile)):
    output = open(outfile, 'w')
    output.write('#Date Date/time at which job finished\n')
    output.write('#DataBeta Input beta for the simulated data sample. Will be 0.0 for real data.\n')
    output.write('#N_sims Number of datalike sims that go into the subsequent means\n')
    output.write('#delta_Beta mean difference between large MC sim beta (2.11 for the time being) and the measured beta for the data (not the beta in column 2.\n')
    output.write('#sigma_Beta std. error in the mean of delta_Beta over N_sims sims\n')
    output.write('#Beta_err mean statistical error on beta\n')
    output.write('#meanZ mean photoZ of the large MC sim\n')
    output.write('#sigmaZ std. deviation of the photoZs for the large Sim\n')
    output.write('#sigmaDZ std. deviation of (zSim - zPHOT)\n')
    output.write('#NCC/NTot overall CC Contamination\n')
    output.write('#TypeChoice Internal Diagnostic, check code comments\n')
    output.write('#NNProbCut Threshold for NN probability of Ia\n')
    output.write('#Date \t\tDataBeta N_sims delta_Beta sigma_Beta BetaStatErr meanZ sigmaZ sigmaDz NCC/NTot TypeChoice NNProbCut\n')
else:
    output = open(outfile, 'a')
print 'outfile'
print outfile

cat = RateTest.simcat.Catalog
t = time.strftime('%b-%d-%H:%M')
N_Sims = len(ks)
BetaStdErr = float(BetaSigma[0])/np.sqrt(N_Sims)
meanZ =  np.mean(cat['zPHOT'])
sigZ = np.std(cat['zPHOT'])
sigDZ = np.std(cat['zPHOT'] - cat['SIM_ZCMB'])
contam = np.sum(cat['SIM_TYPE_INDEX'] !=1).astype(float)/ float(cat.shape[0])
output.write('{0}\t\t{1:.2f}\t{2}\t{3:.3f}\t{4:.3f}\t{5:.3f}\t{6:.3f}\t{7:.3f}\t{8:.3f}\t{9:.3f}\t{10}\t{11:.3f}\n'.format(t, trueBeta, N_Sims, BetaMean[0], BetaStdErr, BetaErrs[0],meanZ, sigZ, sigDZ, contam,  RateTest.typeString, NNProbCut))
print "BetaMean[0]"
print BetaMean[0]
print BetaMean

print "Individual Scales"
print CCScaleStorage
print "Individual ScaleErrs"
print CCScaleErrStorage
print "average ScaleErrs"
print np.mean(CCScaleErrStorage)
print "AAA CC Scale means (weighted, unweighted)2"
print np.average(ma.masked_invalid(np.array(CCScaleStorage)), weights = 1.0/ma.masked_invalid(CCScaleErrStorage)**2)
print np.mean(ma.masked_invalid(np.array(CCScaleStorage)))

print "AAA CC Scale stds"
print np.nanstd(np.array(CCScaleStorage))
plt.clf()
hist, bins = np.histogram(CCScaleStorage, bins = np.linspace(0.0, 5.0, 10))
plt.step((bins[1:]+bins[:-1])/2.0, hist, where = 'mid', c = 'g')
plt.savefig(dataname + 'ScaleDistro.png')
plt.clf()


print "nIter"
print nIter
'''
        argList = ''
        minObjList = ''
        chi2Initargs = ''
        for i in xrange(zCenters.shape[0]):
            argList += 'f{0},'.format(i)
            minObjList += 'f{0} = 1.0, error_f{0} = 0.1, limit_f{0} = (0.0, None),'.format(i)
            chi2Initargs += '1.0,'
        argList = argList[:-1]
        minObjList = minObjList[:-1]
        chi2Initargs = chi2Initargs[:-1]
        #print argList
        #print minObjList
        #print chi2Initargs

        exec('''
'''
def chi2func(nData, nSim, effmat, fnorm, zCenters, {0}, dump = False, complexdump = False):

    Chi2Temp = 0.0
    f_Js = [{0}]
    chi2Mat = np.zeros((self.nbins))
    adjNMC = np.zeros((self.nbins))
    #print f_Js
    #Check if I am scaling errors down with increasing MC size. Make MC twice as large as "Data" to test.
    for row, nDataI, i in zip(effmat, nData, xrange(self.nbins)):
        #if dump:
        #    print "nDataI"
        #    print nDataI
        JSumTemp = 0.0
        for eff, nSimJ, f_J, j in zip(row, nSim, f_Js, xrange(self.nbins)):
            JSumTemp += nSimJ*f_J*eff*fnorm
            if dump and i == j:
                print "nDataI"
                print nDataI
                print "Bin Contribution to scaled nSim"
                print nSimJ*f_J*eff*fnorm
                #print "Product of nSimJ, f_J, eff, fnorm"
                #print nSimJ
                #print f_J
                #print eff
                #print fnorm
        if nDataI > 1E-11 or JSumTemp > 1E-11:
            if dump and i == j:
                print "nDataI"
                print nDataI
                print "scaled nSim"
                print JSumTemp
                print "fnorm"
                print fnorm
                print "error"
                print nDataI + JSumTemp*fnorm
                if (nDataI + JSumTemp*fnorm) <= 0:
                    print (nDataI + JSumTemp*fnorm)
                    assert(0)
            Chi2Temp += ((nDataI - JSumTemp)**2/(nDataI + JSumTemp*fnorm))#*fnorm**2

    return Chi2Temp
                ''''''.format(argList), locals())
        fnorm = float(np.sum(nData))/float(self.simcat.Catalog['zPHOT'].shape[0])

        #print type(chi2func)
        #print 'lamChi2 = lambda {0}: chi2func(nData, nSim, self.effmat, fnorm, zCenters, {0})'.format(argList)
        exec('lamChi2 = lambda {0}: chi2func(nData, nSim, self.effmat, fnorm, zCenters, {0})'.format(argList),locals())
        exec('lamChi2Dump = lambda {0}: chi2func(nData, nSim, self.effmat, fnorm, zCenters, {0}, dump = True)'.format(argList),locals())
        #print type(lamChi2)
        #print type(lamChi2Dump)
        #print 'MinObj = M(lamChi2, {0})'.format(minObjList)
        exec('MinObj = M(lamChi2, {0})'.format(minObjList),locals())
        exec('chi2Init = lamChi2Dump({0})'.format(chi2Initargs),locals())
        #print "Chi2 init = {0}".format(round(chi2Init, 4))



        MinObj.set_strategy(2)
        MinObj.migrad()
        #MinObj.minos()
        zCenters = (simBins[1:] + simBins[:-1])/2.0
        print MinObj.values
        fJs = []
        fJErrs = []
        for v in MinObj.values.keys():
            fJs.append(MinObj.values[v])
            fJErrs.append(MinObj.errors[v])

        
        exec('lamChi22 = lambda k, Beta: self.chi2V2(fJs, fJErrs, zCenters, k, Beta)',locals())
        exec('MinObj2 = M(lamChi22, k = 1.0, error_k = 0.1, limit_k = (0.0, None), Beta = 0.0, error_Beta = 0.1)',locals())


#print "Large Perfect Sim {0}".format(simInd)
    #print "Sim R0 = 1.7E-5; Sim Beta = 4.2"
    ##print "Sim Beta = 1.5; Data Beta = 1.5"
    ##RateTest = Rate_Fitter('DES_FULLSURVEY_TEST/JLDESFULLSURVEYIaOnly+zPHOT+smearC11/FITOPT000+SALT2mu.FITRES', 'JLDESFULLSURVEYIaOnly+zPHOT+smearC11','JLDES_R0_7E-5_Beta_1-5_Shallow/JLDES_R0_7E-5_Beta_1-5_Shallow/FITOPT000+SALT2mu.FITRES', 'JLDES_R0_7E-5_Beta_1-5_Shallow','/project/rkessler/SN/SNDATA_ROOT/SIM/JLDES_R0_7E-5_Beta_1-5_Shallow/JLDES_R0_7E-5_Beta_1-5_Shallow.DUMP')
    #print '/project/rkessler/jlasker/Rate_Analysis/TestSameK2Beta/outFit_datasize/JLDES_R0_1-7E-5_Beta_4-2_Datasize_Perfect-00{0:02d}/FITOPT000.FITRES'.format(simInd)

    #RateTest = Rate_Fitter('/project/rkessler/jlasker/Rate_Analysis/TestSameK2Beta/outFit_datasize/JLDES_R0_1-7E-5_Beta_4-2_Datasize_Perfect-00{0:02d}/FITOPT000.FITRES'.format(simInd), 'TestSameK2Beta/JLDES_R0_1-7E-5_Beta_4-2-00{0:02d}'.format(simInd),'/project/rkessler/jlasker/Rate_Analysis/outFit_datalike/JLDES_R0_1-7E-5_Beta_2-1_Datalike_PERFECT/FITOPT000.FITRES', 'JLDES_R0_1-7E-5_Beta_2-1_DataLikePhotZ','/scratch/midway2/rkessler/SNDATA_ROOT/SIM/JLDES_R0_1-7E-5_Beta_2-1_Datalike_PERFECT/JLDES_R0_1-7E-5_Beta_2-1_Datalike_PERFECT.DUMP', 2.1, zmin = 0.1, zmax = 1.3)# , MJDMin = 0, MJDMax = np.inf)


    #RateTest.effCalc(nbins = 12)
    ##RateTest.effCalc(nbins = 20)
    #RateTest.fit_rate()


    #ksPerf.append(RateTest.k)
    #kErrsPerf.append(RateTest.kErr)
    #BetasPerf.append(RateTest.Beta)
    #BetaErrsPerf.append(RateTest.BetaErr)
    #print "Sim Beta = 1.5; Data Beta = 1.5"
    #RateTest = Rate_Fitter('DES_FULLSURVEY_TEST/JLDESFULLSURVEYIaOnly+zPHOT+smearC11/FITOPT000+SALT2mu.FITRES', 'JLDESFULLSURVEYIaOnly+zPHOT+smearC11','JLDES_R0_7E-5_Beta_1-5_Shallow/JLDES_R0_7E-5_Beta_1-5_Shallow/FITOPT000+SALT2mu.FITRES', 'JLDES_R0_7E-5_Beta_1-5_Shallow','/project/rkessler/SN/SNDATA_ROOT/SIM/JLDES_R0_7E-5_Beta_1-5_Shallow/JLDES_R0_7E-5_Beta_1-5_Shallow.DUMP')


    try:
        optfname = argv[1]
        opts = open(optfname, 'r')
        optlist = opts.readlines()

        zmin = None; zmax = None; MJDMin = None; MJDMax = None; bins = None; runFit = None

        for opt in optlist:
            try: 
                optName, optVal = opt.split()
            except:
                print "{0} not formatted correctly".format(opt)
                continue

            if (optName.lower() == 'zmin') & (not zmin): zmin = optVal
            if (optName.lower() == 'zmax') & (not zmax): zmax = optVal
            if (optName.lower() == 'mjdmin') & (not MJDMin): MJDMin = optVal
            if (optName.lower() == 'mjdmax') & (not MJDMax): MJDMax = optVal
            if (optName.lower() == 'bins') & (not bins): zmin = optVal
            if (optName.lower() == 'runfit') & (not runFit == None): zmin = optVal

        if zmin == None: zmin = 0.1
        if zmax == None: zmax = 1.2
        if MJDMin == None: MJDMin = 0.0
        if MJDMax == None: MJDMax = np.inf
        if bins == None: bins = "equalSize"
        if runFit == None: runFit = True

    except:
        print "Option File not working/Nonexistent. Using default values"
'''