#!/software/python-2.7-2014q3-el6-x86_64/bin/python
import SNANA_Reader as simread
import REAL_Reader as dataread
#import astropy.cosmology as cosmo
import traceback
import scipy
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import Cosmology
import scipy.stats.mstats as mstats
import scipy.stats as stats
from scipy.interpolate import UnivariateSpline
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
import pandas as pd


class Rate_Fitter:
    def __init__(self, realfilename, realName, simfilename, simName, simgenfilename, MCBeta, MCK, zminSamp=0.1, zmaxSamp=1.20 , zminFit = 0.1, zmaxFit = 1.20, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95,  Rate_Model = 'powerlaw',  cheatType = False, cheatZ = False, cheatCCSub = False, cheatCCScale = False, cuts = None, nprint = 5, MURESCuts = None, noCCMC = False, priorRate = None, priorZEff = None, ratePriorErrUp = None, ratePriorErrDown =None, ratePriorErrAll = None, fixCCScale = False):
        print "Rate_Fitter"
        print "np version {0}".format(np.__version__)
        
        self.zminSamp = zminSamp
        self.zmaxSamp = zmaxSamp
        self.zminFit = zminFit
        self.zmaxFit = zmaxFit
        self.MCBeta = MCBeta
        self.MCK = MCK
        self.Rate_Model = Rate_Model
        self.cheatType = cheatType
        self.cheatZ = cheatZ
        self.cheatCCSub = cheatCCSub
        self.cheatCCScale = cheatCCScale
        self.cuts = cuts
        self.nprint = nprint
        self.MURESCuts = MURESCuts
        self.priorRate = priorRate
        self.priorZEff = priorZEff
        self.ratePriorErrUp = ratePriorErrUp
        self.ratePriorErrDown = ratePriorErrDown
        self.ratePriorErrAll = ratePriorErrAll
        self.fixCCScale = fixCCScale

        #print "PRIORS"
        #print priorRate
        #print priorZEff
        #print ratePriorErrUp
        #print ratePriorErrDown

        if self.cheatZ:
            self.ztype = 'SIM_ZCMB'
        else:
            #self.ztype = 'zHD'
            self.ztype = 'zPHOT'

        self.shiftFlagData = False
        self.shiftFlagSim = False


        self.globalChi2Storage = []
        self.globalNDataStorage = []
        '''
        
        self.globalZPhotBinStorage = []
        self.globalNDataIaPhotBinStorage = []
        self.globalNDataCCPhotBinStorage = []
        self.globalZTrueBinStorage = []
        self.globalNDataIaTrueBinStorage = []
        self.globalNDataCCTrueBinStorage = []
        '''
        print 'a'
        try: 
            self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
        except:
            try:
                self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95, skip_header = 5)

            except: 
                self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95, skip_header = 6)
        print 'b'   
        self.simName = simName
        self.simgencat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
        print 'c'   
        try:
            #with np.load(simgenfilename+'.npz', allow_pickle = True) as data0:
            #    SIMGEN = data0['a']
        
            SIMGEN = np.load(simgenfilename + '.npy', allow_pickle = True)
        except:
        
            SIMGEN = np.genfromtxt(simgenfilename, dtype=None, names = True, skip_footer=3, invalid_raise=False)
            print "Compress save A"
            SIMGEN.dtype.names = map(str, SIMGEN.dtype.names)
            #np.savez_compressed(simgenfilename+'.npz', a = SIMGEN)
            np.save(simgenfilename+'.npy', SIMGEN)
        
            print "WHY DO YOU HATE ME WHEN I SHOW YOU NOTHING BUT LOVE"
            print simgenfilename
        #SIMGEN = pd.read_csv(simgenfilename, delim_whitespace=True, comment="#").to_records(index = False)
        print 'd'
        SIMGEN = SIMGEN[SIMGEN['GENZ'] != 'GENZ']

        self.simgencat.params = {'flat':True, 'H0': simH0, 'Om0':simOmegaM, 'Ob0': simOb0, 'sigma8': simSigma8, 'ns': simNs}
        self.simgencat.cosmo = Cosmology.setCosmology('simCosmo', self.simcat.params)
        self.simgencat.OrigCatalog = np.copy(SIMGEN)
        self.simgencat.Catalog = np.copy(SIMGEN)
        self.simgencat.Catalog = self.simgencat.Catalog[self.simgencat.Catalog['GENZ'] != 'GENZ']
        self.simgencat.simname = simName
        self.simgencat.NSN = self.simgencat.Catalog['GENZ'].shape[2]

        print "SIMGEN NUMBER"
        print self.simgencat.NSN
        print "TEST2"
        print self.simgencat.Catalog['GENZ'].shape[0]
        print self.simgencat.Catalog['GENZ'].shape[1]
        print self.simgencat.Catalog['GENZ'].shape[2]
        print "SIMGENCAT FILE"
        print simfilename

        self.realName = realName
        try:
            print 'q'
            self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95, skip_header = 6)
        except:
            #self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
            try:
                print 'r'
                self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
            except:
                print 's'
                self.realcat = dataread.REAL_Cat(realfilename, realName, skip_header =11 )

        if self.cheatType:
            print "WARNING, THE FITTER IS CHEATING AND ELIMINATED NON-IAs USING SIM INFO"
            self.realcat.Catalog = self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1]
            self.simcat.Catalog = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1]

        print "Pre cut Catalog"
        print self.realcat.Catalog.shape
        for cut in cuts:
            print 'a'
            print cut
            print self.realcat.Catalog.shape
            self.realcat.Catalog = self.realcat.Catalog[(self.realcat.Catalog[cut[0]].astype(type(cut[1])) > cut[1]) & (self.realcat.Catalog[cut[0]].astype(type(cut[2])) < cut[2])]
            self.simcat.Catalog = self.simcat.Catalog[(self.simcat.Catalog[cut[0]].astype(type(cut[1])) > cut[1]) & (self.simcat.Catalog[cut[0]].astype(type(cut[2])) < cut[2])]
            print 'b'
            print cut
            print self.realcat.Catalog.shape

        self.postCutRealCat = np.copy(self.realcat.Catalog)
        self.postCutSimCat = np.copy(self.simcat.Catalog)

        self.realcat.Catalog = self.realcat.Catalog[(self.realcat.Catalog[self.ztype].astype(float) > self.zminSamp) & (self.realcat.Catalog[self.ztype].astype(float) < self.zmaxSamp)]
        self.simcat.Catalog = self.simcat.Catalog[(self.simcat.Catalog[self.ztype].astype(float) > self.zminSamp) & (self.simcat.Catalog[self.ztype].astype(float) < self.zmaxSamp)]
        print 'zCut Pre MURESCut'
        print np.sum((self.realcat.Catalog[self.ztype].astype(float) > self.zminFit) & (self.realcat.Catalog[self.ztype].astype(float) < self.zmaxFit))
        print 'MURESCUT'
        print self.MURESCuts
        print self.realcat.Catalog.shape

        if not (self.MURESCuts is None):
            '''
            #MURES Cut format: (zmin, zmax, neg Cut, pos Cut)

            for mc in self.MURESCuts:

                realCond = (self.realcat.Catalog[self.ztype] < mc[0]) | (self.realcat.Catalog[self.ztype] > mc[1])| ((self.realcat.Catalog['MURES'] > mc[2])& (self.realcat.Catalog['MURES'] < mc[3]))

                simCond = (self.simcat.Catalog[self.ztype] < mc[0]) | (self.simcat.Catalog[self.ztype] > mc[1])| ((self.simcat.Catalog['MURES'] > mc[2])& (self.simcat.Catalog['MURES'] < mc[3]))

                self.realcat.Catalog = self.realcat.Catalog[realCond]
                self.simcat.Catalog = self.simcat.Catalog[simCond]
                '''

            self.realcat.Catalog = self.realcat.Catalog[ np.abs( self.realcat.Catalog['MURES'] * 1.0 / self.realcat.Catalog['MUERR'] ) < MURESCuts]
            self.simcat.Catalog = self.simcat.Catalog[ np.abs( self.simcat.Catalog['MURES'] * 1.0 / self.simcat.Catalog['MUERR'] ) < MURESCuts]
        print "PostMURESCut Shape"
        print self.realcat.Catalog.shape
        print 'zCut Post MURESCut'
        print np.sum((self.realcat.Catalog[self.ztype].astype(float) > self.zminFit) & (self.realcat.Catalog[self.ztype].astype(float) < self.zmaxFit))

        print "Post cut Catalog"

        print self.realcat.Catalog.shape

        if noCCMC:
            self.simgencat.Catalog = self.simgencat.Catalog[self.simgencat.Catalog['GENTYPE'] == 1]
            self.simcat.Catalog = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1]

        
        
    def newData(self, realfilename, realName, simInd =100):
        self.realName = realName
        self.shiftFlagData = False
        try:
            self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
        except:
            self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95, skip_header = 6 )
        if self.cheatType:
            print "WARNING, THE FITTER IS CHEATING AND ELIMINATED NON-IAs USING SIM INFO"
            self.realcat.Catalog = self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1]
        if simInd < self.nprint:
            print 'N precuts'
            print self.realcat.Catalog['FITPROB'].shape
        print "Pre cut Catalog"
        print self.realcat.Catalog.shape

        for cut in cuts:
            self.realcat.Catalog = self.realcat.Catalog[(self.realcat.Catalog[cut[0]].astype(type(cut[1])) > cut[1]) & (self.realcat.Catalog[cut[0]].astype(type(cut[2])) < cut[2])]

        self.realcat.Catalog = self.realcat.Catalog[(self.realcat.Catalog[self.ztype].astype(float) > self.zminSamp) & (self.realcat.Catalog[self.ztype].astype(float) < self.zmaxSamp)]
        print "Post cut Catalog"
        print self.realcat.Catalog.shape   
        

        self.postCutRealCat = np.copy(self.realcat.Catalog)
        print 'MURESCUT'
        print self.MURESCuts
        print self.realcat.Catalog.shape
        if not (self.MURESCuts is None):
        
            #MURES Cut format: (zmin, zmax, neg Cut, pos Cut)
            '''
            for mc in self.MURESCuts:
            
                realCond = (self.realcat.Catalog[self.ztype] < mc[0]) | (self.realcat.Catalog[self.ztype] > mc[1])| ((self.realcat.Catalog['MURES'] > mc[2])& (self.realcat.Catalog['MURES'] < mc[3]))

                self.realcat.Catalog = self.realcat.Catalog[realCond]
            '''
            self.realcat.Catalog = self.realcat.Catalog[np.abs(self.realcat.Catalog['MURES']*1.0/self.realcat.Catalog['MUERR']) < MURESCuts]
        print "PostMURESCut Shape"
        print self.realcat.Catalog.shape

        
        if simInd < self.nprint:
            print "Minimum Fitprob"
            print np.min(self.realcat.Catalog['FITPROB'])
            print 'N postcuts'
            print self.realcat.Catalog['FITPROB'].shape

    def zSystematic(self, binList = None, nbins = None):
        assert(0)
        if nbins is None:
            try: 
                self.nbins = len(binList) - 1
                self.binList = binList
            except:
                self.nbins = binList.shape[0] - 1
                self.binList = binList
        else:
            binList = np.linspace(self.zmin, self.zmax, nbins+1)
            self.nbins = nbins
            self.binList = binList
        if self.shiftFlagData:
            print "DONT DOUBLE SHIFT"
            return 0
        if not self.shiftFlagSim:
        
            oldsimz = self.simcat.Catalog['zPHOT']
            oldsimtruez = self.simcat.Catalog['SIM_ZCMB']
            stat, bins, binnum = stats.binned_statistic(oldsimz, oldsimz - oldsimtruez, bins = self.binList, statistic = 'mean')
            self.zBiasShifts = stat
            newsimz = oldsimz - stat[binnum]
            assert(np.sum(np.abs(newsimz - oldsimz)) > 0)
            assert((oldzshape - np.arange(0, oldz.shape[0]).shape[0])< 1)
            self.shiftFlagSim = True
        oldz = self.realcat.Catalog['zPHOT']
        _,_, binnum = stats.binned_statistic(oldz, oldz , bins = self.binList, statistic = 'mean')
        newz = oldz - self.zBiasShifts[binnum]
        oldzshape = oldz.shape[0]
        self.realcat.Catalog['zPHOT'].put(np.arange(0, oldz.shape[0]), newz)
        assert(np.sum(np.abs(newz - oldz)) > 0)
        assert((oldzshape - np.arange(0, oldz.shape[0]).shape[0])< 1)
        self.simFlagData = True
        
    def effCalc(self, fracContamCut = 0.0, nbinsSamp = None, nbinsFit = None, binListSamp = None, binListFit = None, simInd =100):
        #### Do we want SNIas or all SN for efficiency?
        import matplotlib as mpl
        if nbinsSamp is None:
            try: 
                self.nbinsSamp = len(binListSamp) - 1
                self.binListSamp = binListSamp
            except:
                self.nbinsSamp = binListSamp.shape[0] - 1
                self.binListSamp = binListSamp
        else:
            binListSamp = np.linspace(self.zminSamp, self.zmaxSamp, nbinsSamp+1)
            self.nbinsSamp = nbinsSamp
            self.binListSamp = binListSamp

        if nbinsFit is None:
            try: 
                self.nbinsFit = len(binListFit) - 1
                self.binListFit = binListFit
            except:
                self.nbinsFit = binListFit.shape[0] - 1
                self.binListFit = binListFit
        else:
            binListFit = np.linspace(self.zminFit, self.zmaxFit, nbinsFit+1)
            self.nbinsFit = nbinsFit
            self.binListFit = binListFit

        
        self.typeString = ''

        #if self.cheatZ:
        #    self.ztype = 'SIM_ZCMB'
        #else:
        #    self.ztype = 'zPHOT'

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
        zPHOTs = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1][self.ztype].astype(float)

        zTRUEs = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) == 1]['SIM_ZCMB'].astype(float)

        self.typeString = self.typeString + 'A1'
        
        
        print "Type Location A"
        print "Choice A1"
        print zPHOTs.shape
        print zTRUEs.shape
        print binList
        
        counts, zPhotEdges, zTrueEdges, binnumber = scipy.stats.binned_statistic_2d(zPHOTs, zTRUEs, zTRUEs, statistic = 'count', bins =  (self.binListFit, self.binListSamp))
        assert(zPhotEdges.shape[0] == (self.nbinsFit + 1))
        print "Type Location B"
        print "Choice B1"
        
        self.typeString = self.typeString + 'B1'
        zGenHist, zGenBins = np.histogram(self.simgencat.Catalog[self.simgencat.Catalog['GENTYPE'].astype(int) == 1]['GENZ'].astype(float), bins = self.binListSamp)

        #zSim1Hist, zSim1Bins = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) ==1]['SIM_ZCMB'].astype(float), bins = self.binListSamp)
        
       
        
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
        self.effmat = np.zeros((self.nbinsFit, self.nbinsSamp))
        xMax = zPhotEdges.shape[0] - 2
        yMax = zTrueEdges.shape[0] - 2
        print zGenHist
        print counts.astype(int)
        '''
        for zPhotLedge, zPhotRedge, row, i in zip(zPhotEdges[:-1], zPhotEdges[1:], counts, xrange(xMax + 1)):
            zPhotCenter = (zPhotLedge + zPhotRedge)/2.0
            for zTrueLedge, zTrueRedge, count, j in zip(zTrueEdges[:-1], zTrueEdges[1:], row, xrange(yMax + 1)):
                zTrueCenter = (zTrueLedge + zTrueRedge)/2.0
                inCell = (zPHOTs > zPhotLedge) & (zPHOTs < zPhotRedge) & (zTRUEs > zTrueLedge)& (zTRUEs < zTrueRedge)
                zPhotCell = zPHOTs[inCell];zTrueCell = zTRUEs[inCell]
                self.effmat[i][j] = count # np.sum(inCell)
                #print "inCell"
                #print np.sum(inCell)
                #print "count"
                #print count
                #try:
                #    assert(np.abs(np.sum(inCell) - count) < 2)
                #except:
                #    print "CHECK ABOVE"
                
        for row, i in zip(self.effmat, xrange(self.effmat.shape[0])):
            for j in xrange(row.shape[0]):
                self.effmat[i][j] /= zGenHist[j]
        '''
        self.effmat = counts/zGenHist

        #if simInd < self.nprint:
        print 'effmat'
        print self.effmat




        extent = [zPhotEdges[0], zPhotEdges[-1], zTrueEdges[0], zTrueEdges[-1]]
        if simInd == 1:
            plt.figure()
            plt.imshow(np.flipud(counts.T), extent = extent, cmap = 'Blues')
            plt.colorbar()
            plt.savefig(self.realName + 'redshiftDistro.png')
            plt.clf()
            plt.close()
            plt.figure()
            plt.imshow(np.flipud(self.effmat.T), extent = extent, cmap = 'Blues', norm=mpl.colors.LogNorm())
            plt.colorbar()
            plt.savefig(self.realName + 'efficiencyMatrixLog.png')
            plt.clf()
            plt.close()
            plt.figure()
            plt.imshow(np.flipud(self.effmat.T), extent = extent, cmap = 'Blues')
            plt.colorbar()
            plt.savefig(self.realName + 'efficiencyMatrix.png')
            plt.clf()
            plt.close()
            
    def fit_rate(self, fixK = False, fixBeta = False, simInd =100, trueBeta = 0, CCScale = 1.0, CCScaleErr = None, TrueCCScale = 1.0, BetaInit = 0.0, kInit = 1.0, BetaErr = 1, kErr = 1, f_Js = None, CCZbins = None, scaleZBins = None, Blind = False):
        #import iminuit as iM
        #from iminuit import Minuit as M
        #import numpy as np
        #import matplotlib as mpl
        #import matplotlib.pyplot as plt
        #if self.cheatZ:
        #    self.ztype = 'SIM_ZCMB'
        #else:
        #    self.ztype = 'zPHOT'
        plt.switch_backend('Agg')

        if simInd < self.nprint:
            print "Type Location C"
            print "Choice C1"

        if len(self.typeString) <= 4:
            self.typeString = self.typeString + 'C1'


        nSim, simBins = np.histogram(self.simgencat.Catalog[self.simgencat.Catalog['GENTYPE'].astype(int) == 1]['GENZ'].astype(float), bins=self.binListSamp)
        if simInd < self.nprint:
            print "nSim1"
            print nSim
            print self.simgencat.Catalog.shape
        
        print "FIGURE OUT WHY YOU MADE THIS ASSERT STATEMENT LATER"
        #assert(0)
        nSim2, simBins2 = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) ==1][self.ztype].astype(float), bins=self.binListFit)
        
        
       
        nSim3, simBins3 = np.histogram(self.simcat.Catalog[self.ztype].astype(float), bins=self.binListFit)
        

        NCC , _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1][self.ztype].astype(float), bins=self.binListFit)
        if self.fixCCScale:
            print "Fix CC Scale at 1"
        else:
            if simInd < self.nprint:
                print "nSim2"
                print nSim2
                print "nSim3"
                print nSim3
                print "nCC"
                print NCC
            OrigNCC = np.copy(NCC)
            if self.cheatCCSub:
                if self.cheatCCScale:
                    print "WARNING: Only cheating on CC Subtraction not scale"
                print "Setting NCC to infinity to make sure that cheating correctly"
                print "Diagnostics after this point may be nonsense"
                print self.cheatCCSub
                print "NCC BeforeFck"
                print NCC
                NCC = NCC*1E100
                print "NCC AfterFck"
                print NCC    
            elif self.cheatCCScale:
                print "NCC Before1"
                print NCC
                print TrueCCScale
                NCC = applyCCScale(NCC, TrueCCScale, CCScaleErr, zbins = CCZbins, datazbins = self.binListFit)
                print "NCC After1"
                print NCC
            else: 
                print "NCC Before2"
                print NCC
                print CCScale
                NCC = applyCCScale(NCC, CCScale, CCScaleErr, zbins = CCZbins, datazbins = self.binListFit)
                print "NCC After2"
                print NCC
            #assert(0)

        
        NIa , _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1][self.ztype].astype(float), bins=self.binListFit)
        '''
        DebugNIaPhot, _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1]['zPHOT'].astype(float), bins=self.binListFit)
        DebugNCCPhot, _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1]['zPHOT'].astype(float), bins=self.binListFit)
        DebugNCCPhot = applyCCScale(DebugNCCPhot, CCScale, CCScaleErr, zbins = scaleZBins, datazbins = self.binListFit)
        DebugNIaTrue, _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1]['SIM_ZCMB'].astype(float), bins=self.binListSamp)
        DebugNCCTrue, _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1]['SIM_ZCMB'].astype(float), bins=self.binListSamp)
        DebugNCCTrue = applyCCScale(DebugNCCTrue, CCScale, CCScaleErr, zbins = scaleZBins, datazbins = self.binListSamp)

        uselessCtr = 0
        for niap, nccp, niat, ncct, zp, zt in zip(DebugNIaPhot, DebugNCCPhot, DebugNIaTrue, DebugNCCTrue,(self.binListFit[1:] + self.binListFit[:-1])/2.0, (self.binListSamp[1:] + self.binListSamp[:-1])/2.0 ):
            uselessCtr +=1
            self.globalZTrueBinStorage.append(zt)
            self.globalZPhotBinStorage.append(zp)
            self.globalNDataIaPhotBinStorage.append(niap)
            self.globalNDataCCPhotBinStorage.append(nccp)
            self.globalNDataIaTrueBinStorage.append(niat)
            self.globalNDataCCTrueBinStorage.append(ncct)
        print "UselessCtr"
        print uselessCtr
        
        '''

        try:
            TrueNCC, _ = np.histogram(self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'] !=1][self.ztype].astype(float), bins=self.binListFit)
            if simInd < self.nprint:

                print "True NCC Data"
                print TrueNCC
        except:
            print "Using real data"

            TrueNCC = 0.0

        nData, dataBins = np.histogram(self.realcat.Catalog[self.ztype].astype(float), bins=self.binListFit)
        print "nData"
        print nData
        if not(self.cheatCCSub):
            FracBad = NCC*1.0/(1.0*(NCC+NIa))
            nCCData = nData*FracBad
        else: 
            nCCData = TrueNCC*1.0
            FracBad = TrueNCC*1.0/nData
        if simInd < self.nprint:
            print "PreScale NCC/nSim"
            print OrigNCC*1.0/(OrigNCC+NIa)
            
            print "PreScale Pred NCC Data"
            print OrigNCC*1.0/(OrigNCC+NIa)*nData

            print "PreScale Pred NCC Data if 2NCC"
            print OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData

            print "TrueNCC"
            print TrueNCC
        if type(TrueNCC) != int:
            if simInd < self.nprint:
                print "PreScale PredNCCData - TrueNCCData"
                print OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData - TrueNCC

                print "PreScale PredNCCData - TrueNCCData/ PredNCCData"
                print (OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData - TrueNCC)/(OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData)
        else:
            print "Using real data"
        
        print "Mean of PreScale PredNCCData - TrueNCCData/ PredNCCData"
        print np.nanmean((OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData - TrueNCC)/(OrigNCC*2.0/(2.0*OrigNCC+NIa)*nData))

        print "PostScale NCC/nData"
        print NCC*1.0/(NCC+NIa)

        if simInd < self.nprint:
            print "Fraction of CCs in each bin"
            print FracBad

            print 'NCC'
            print NCC

            print 'nSim2'
            print nSim2
            print "nData, dataBins, realcat shape pre contam correction"
            print nData
            print dataBins
            print np.sum(self.realcat.Catalog[self.ztype].astype(float) > self.zmaxFit)
            print np.sum(self.realcat.Catalog[self.ztype].astype(float) < self.zminFit)
            print self.realcat.Catalog[self.ztype].shape
            
            print "Ratio nData/nSim"
            print 1.0*nData/(1.0*nSim3)
    

            print "Ratio nSim2/nData"
            print 1.0*nSim3/(1.0*nData)

            print "FracBad"
            print FracBad
            print 'NCCData'
            print nCCData

        if simInd < self.nprint:

            print "overall Contam"
            print np.sum(NCC)*1.0/(np.sum(nSim3)*1.0)
        
        def chi2func(nData, nSim, effmat, fnorm, zCentersSamp, zCentersFit, k = 1.0, Beta = 0.0, zBreak = 1.0, dump = False, complexdump = False, modelError = False, nIA = None, nCC = None, Rate_Model = 'powerlaw', zbins = None, simInd = 100, BetaPrior = (-3, 3), KPrior = (0.0, 50.0), priorRate = None, priorZEff = None, ratePriorErrUp = None, ratePriorErrDown =None, ratePriorErrAll = None,  TrueNCCData = None, f_1 = 1.0, f_2 = 1.0, f_3 = 1.0, f_4 = 1.0, f_5 = 1.0, f_6 = 1.0, f_7 = 1.0, f_8 = 1.0, f_9 = 1.0, f_10 = 1.0, f_11 = 1.0):
            print "PRIORS2"
            print priorRate
            print priorZEff
            print ratePriorErrUp
            print ratePriorErrDown
            Chi2Temp = 0.0
            if Rate_Model == 'powerlaw':
                f_Js = k*(1+zCentersSamp)**Beta
            elif Rate_Model == 'discrete':
                f_Js = np.array([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11])
            elif (Rate_Model == 'brokenpowerlaw') | (Rate_Model == 'brokenpowerlawVar'):
                f_Js = []
                #zCenters = (zbins[1:]+zbins[:-1])/2.0
                temp = None
                for zC in zCentersSamp:
                    if zC < zBreak:
                        f_Js.append(k*(1+zC)**Beta)
                    elif not(temp is None):
                        f_Js.append(temp)
                    else:
                        temp = f_Js[-1]
                        f_Js.append(temp)
                f_Js = np.array(f_Js)
            else: 
                assert(0)
            if simInd >  self.nprint:
                if Rate_Model == 'discrete':
                    print "f_Js init"
                    print f_Js
                else:
                    print "Beta init"
                    print Beta
                    print "k init"
                    print k
            #chi2Mat = np.zeros((self.nbinsFit))
            #adjNMC = np.zeros((self.nbinsFit))
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
                    if dump:
                        print 'f_J adjusted CC Cut'
                        print Rate_Model
                        print nCC
                        print nIA
                        print np.array(f_Js)
                    fracCCData = (nCC*1.0)/((1.0*nCC + nIA*np.array(f_Js)))
                    print fracCCData
                else:
                    if dump:
                        print "Beta Adjusted CC Cut"
                        print Rate_Model
                    #BetaRatio = k*(1+zCenters)**(Beta)#/(1+zCenters)**MCBeta
                    BetaRatio = (1+zCentersFit)**(Beta)#/(1+zCenters)**MCBeta
                    if dump:
                        print "Beta Ratio"
                        print BetaRatio
                        print "BadFracCCData"
                        print (nCC*1.0)/((1.0*nCC + nIA*BetaRatio))
                        print "bad NCCData"
                        print (nCC*1.0)/((1.0*nCC + nIA*BetaRatio))*nData
                    fracCCData = (nCC*1.0)/((1.0*nCC + nIA*BetaRatio))
            
            

            if dump and (self.nprint > simInd):
                print 'abc'
                print "fracCCData2"
                print fracCCData
                print "unscaled fracCCData"
                print (1.0*nCC)/(1.0*(nCC+nIA))
            if self.cheatCCSub:
                nCCData = TrueNCCData
                if dump and (self.nprint < simInd):

                    print "Cheating CC Sub"
                    assert(not(TrueNCCData is None))

            elif dump and (self.nprint > simInd):
                print 'def'
                print "Normal CC Sub"
            if not self.cheatCCSub:
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
            if dump: JSumTempNumStor = []
            if dump: JSumTempDenStor = []

            if dump:
                print "actually used NCC"
                #print nCC
                print nCCData
            if simInd < self.nprint:
                print "effmat"
                print effmat
                print "nData"
                print nData
                print "nCCData"
                print nCCData
                print "nSim"
                print nSim

            print nCCData
            for row, nDataI, nCCDataI, i, zc in zip(effmat, nData, nCCData, range(self.nbinsFit), zCentersFit):
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
                if simInd < self.nprint:
                    print "nBinsSamp"
                    print self.nbinsSamp
                assert(row.shape[0] == self.nbinsSamp)
                assert(nSim.shape[0] == self.nbinsSamp)
                assert(len(f_Js) == self.nbinsSamp)
                for eff, nSimJ, f_J, j in zip(row, nSim, f_Js, range(self.nbinsSamp)):
                    if dump  and (self.nprint < simInd):
                        print 'NGen J'
                        print nSimJ
                        print 'JSumTempNum contr'
                        print nSimJ*f_J*eff*fnorm
                        print 'JSumTempDen contr'
                        print nSimJ*f_J*eff*fnorm*f_J*fnorm
                    #if dump and (i != j) and self.cheatZ and (self.nprint < simInd):
                    #    if nSimJ*f_J*eff*fnorm > 0:
                    #        print " This should be zero but isnt "
                    #        print nSimJ*f_J*eff*fnorm
                    #        assert(0)
                    JSumTempNum += nSimJ*f_J*eff*fnorm
                    JSumTempDen += nSimJ*f_J*eff*fnorm*f_J*fnorm
                dataFunc = np.maximum(nDataI ,1)
                #CCFunc = np.ceil(np.maximum(nCCDataI, 1))
                CCFunc = np.maximum(nCCDataI, 1)
                c2t = (nDataI - nCCDataI - JSumTempNum)**2/( dataFunc + CCFunc + JSumTempDen) 
                if dump:
                    JSumTempNumStor.append(JSumTempNum)
                    JSumTempDenStor.append(JSumTempDen)

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
            if dump and (self.nprint < simInd):
                print "JSumTempNum/Den"
                print JSumTempNumStor
                print JSumTempDenStor

            if dump:
                if (self.nprint < simInd):
                    print Chi2Temp
                    print kprior
                    print betaprior
                    print chi2Storage

                    
                    print "nData"
                    print nData
                    print "nCCData"
                    print nCCData
                if priorRate is None:

                    return Chi2Temp+kprior+betaprior , chi2Storage 
                else:
                    print "PRIORS3"
                    print priorRate
                    print "fit k"
                    print k
                    print 'MCK'
                    print self.MCK
                    print "fit beta"
                    print Beta
                    print 'MCBeta'
                    print self.MCBeta
                    print ratePrior(k*self.MCK, Beta + self.MCBeta, priorRate, priorZEff, ratePriorErrUp, ratePriorErrDown, ratePriorErrAll)

                    return Chi2Temp+kprior+betaprior + ratePrior(k*self.MCK, Beta+self.MCBeta, priorRate, priorZEff, ratePriorErrUp, ratePriorErrDown, ratePriorErrAll), chi2Storage 
            else:
                if dump and (self.nprint < simInd):
                    print 'C2T'
                    print Chi2Temp
                    print kprior
                    print betaprior

                if priorRate is None:

                    return Chi2Temp+kprior+betaprior 
                else:
                    print "PRIORS3"
                    print priorRate
                    print "fit k"
                    print k
                    print 'MCK'
                    print self.MCK
                    print "fit beta"
                    print Beta
                    print 'MCBeta'
                    print self.MCBeta
                    print ratePrior(k*self.MCK, Beta+self.MCBeta, priorRate, priorZEff, ratePriorErrUp, ratePriorErrDown, ratePriorErrAll)

                    return Chi2Temp+kprior+betaprior + ratePrior(k*self.MCK, Beta+self.MCBeta, priorRate, priorZEff, ratePriorErrUp, ratePriorErrDown, ratePriorErrAll)
        
        zCentersSamp = (self.binListSamp[1:] + self.binListSamp[:-1])/2.0
        zCentersFit = (self.binListFit[1:] + self.binListFit[:-1])/2.0
 
        #Is this right? Everything else in the other side of the chi2 function should be Ia only
        if self.cheatCCSub:
            self.fracCCData = TrueNCC*1.0/nData
        else:
            self.fracCCData = (NCC*1.0)/(1.0*(NCC + NIa))
        if (self.nprint < simInd):
            print "nSim"
            print nSim
            print 'fracCCData'
            print self.fracCCData
            print "nData"
            print nData
        #fnorm = float(np.sum(nData*(1-self.fracCCData)))/float(np.sum(nSim))
        fnorm = 1.0/240.0
        #print "PRIORS"
        #print self.priorZEff
        #print self.priorRate
        #print self.ratePriorErrUp
        #print self.ratePriorErrDown
        if self.Rate_Model == 'powerlaw':
            lamChi2 = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, k, Beta, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC, priorRate = self.priorRate, priorZEff = self.priorZEff, ratePriorErrUp = self.ratePriorErrUp, ratePriorErrDown =self.ratePriorErrDown, ratePriorErrAll = self.ratePriorErrAll)#, zbins = self.binListFit)
            lamChi2Dump = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, k, Beta, dump = True, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC, priorRate = self.priorRate, priorZEff = self.priorZEff, ratePriorErrUp = self.ratePriorErrUp, ratePriorErrDown =self.ratePriorErrDown, ratePriorErrAll = self.ratePriorErrAll)#, zbins = self.binListFit)
            MinObj = M(lamChi2, k = kInit, error_k = kErr , Beta = BetaInit, error_Beta = BetaErr, limit_k = (0.0, None), limit_Beta = (-100, 100), fix_k = fixK, fix_Beta = fixBeta)
            c2i, _ = lamChi2Dump(1.0, 0.0)

            print "Chi2 init = {0}".format(round(c2i, 4))
        elif self.Rate_Model == 'brokenpowerlaw':
            lamChi2 = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, k, Beta, 1.0, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC, Rate_Model = 'brokenpowerlaw', priorRate = self.priorRate, priorZEff = self.priorZEff, ratePriorErrUp = self.ratePriorErrUp, ratePriorErrDown =self.ratePriorErrDown, ratePriorErrAll = self.ratePriorErrAll)#, zbins = self.binListFit)
            lamChi2Dump = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, k, Beta, 1.0, dump = True, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC, Rate_Model = 'brokenpowerlaw', priorRate = self.priorRate, priorZEff = self.priorZEff, ratePriorErrUp = self.ratePriorErrUp, ratePriorErrDown =self.ratePriorErrDown, ratePriorErrAll = self.ratePriorErrAll)#, zbins = self.binListFit)

            MinObj = M(lamChi2, k = kInit, error_k = kErr , Beta = BetaInit, error_Beta = BetaErr, limit_k = (0.0, None), limit_Beta = (-100, 100), fix_k = fixK, fix_Beta = fixBeta)
            c2i, _ = lamChi2Dump(1.0, 0.0)

            print "Chi2 init = {0}".format(round(c2i, 4))
        elif self.Rate_Model == 'brokenpowerlawVar':
            lamChi2 = lambda k, Beta, zBreak: chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, k, Beta, zBreak, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC, Rate_Model = 'brokenpowerlawVar', priorRate = self.priorRate, priorZEff = self.priorZEff, ratePriorErrUp = self.ratePriorErrUp, ratePriorErrDown =self.ratePriorErrDown, ratePriorErrAll = self.ratePriorErrAll)#, zbins = self.binListFit)
            lamChi2Dump = lambda k, Beta, zBreak: chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, k, Beta, zBreak, dump = True, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC, Rate_Model = 'brokenpowerlawVar', priorRate = self.priorRate, priorZEff = self.priorZEff, ratePriorErrUp = self.ratePriorErrUp, ratePriorErrDown =self.ratePriorErrDown, ratePriorErrAll = self.ratePriorErrAll)#, zbins = self.binListFit)

            MinObj = M(lamChi2, k = kInit, error_k = kErr , Beta = BetaInit, error_Beta = BetaErr, limit_k = (0.0, None), limit_Beta = (-100, 100), fix_k = fixK, fix_Beta = fixBeta, zBreak = 1.0, error_zBreak = 0.1, limit_zBreak = (self.zminFit, self.zmaxFit))
            c2i, _ = lamChi2Dump(1.0, 0.0)

            print "Chi2 init = {0}".format(round(c2i, 4))

            
        elif self.Rate_Model == 'discrete':
            
            lamChi2 = lambda f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11: chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, 1.0, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC, f_1 = f_1, f_2 = f_2,f_3 = f_3, f_4 = f_4,f_5 = f_5, f_6 = f_6,f_7 = f_7, f_8 = f_8,f_9 = f_9, f_10 = f_10, f_11 = f_11, Rate_Model = 'discrete', priorRate = self.priorRate, priorZEff = self.priorZEff, ratePriorErrUp = self.ratePriorErrUp, ratePriorErrDown =self.ratePriorErrDown, ratePriorErrAll = self.ratePriorErrAll)#, zbins = self.binListFit )
            lamChi2Dump = lambda f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11: chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, 1.0, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC, f_1 = f_1, f_2 = f_2,f_3 = f_3, f_4 = f_4,f_5 = f_5, f_6 = f_6,f_7 = f_7, f_8 = f_8,f_9 = f_9, f_10 = f_10, f_11 = f_11, dump = True, Rate_Model = 'discrete', priorRate = self.priorRate, priorZEff = self.priorZEff, ratePriorErrUp = self.ratePriorErrUp, ratePriorErrDown =self.ratePriorErrDown, ratePriorErrAll = self.ratePriorErrAll)#, zbins = self.binListFit)

            c2i, _ = lamChi2Dump(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

            print "Chi2 init = {0}".format(round(c2i, 4))

            MinObj = M(lamChi2, f_1 = 1.0, error_f_1 = 1.0, limit_f_1 = (0.0, None), f_2 = 1.0, error_f_2 = 1.0, limit_f_2 = (0.0, None), f_3 = 1.0, error_f_3 = 1.0, limit_f_3 = (0.0, None), f_4 = 1.0, error_f_4 = 1.0, limit_f_4 = (0.0, None), f_5 = 1.0, error_f_5 = 1.0, limit_f_5 = (0.0, None), f_6 = 1.0, error_f_6 = 1.0, limit_f_6 = (0.0, None), f_7 = 1.0, error_f_7 = 1.0, limit_f_7 = (0.0, None), f_8 = 1.0, error_f_8 = 1.0, limit_f_8 = (0.0, None), f_9 = 1.0, error_f_9 = 1.0, limit_f_9 = (0.0, None), f_10 = 1.0, error_f_10 = 1.0, limit_f_10 = (0.0, None), f_11 = 1.0,error_f_11 = 1.0, limit_f_11 = (0.0, None))

        if self.Rate_Model == 'discrete':
            c2f, c2stor = lamChi2Dump(MinObj.values['f_1'],MinObj.values['f_2'],MinObj.values['f_3'],MinObj.values['f_4'],MinObj.values['f_5'],MinObj.values['f_6'],MinObj.values['f_7'],MinObj.values['f_8'],MinObj.values['f_9'],MinObj.values['f_10'],MinObj.values['f_11'])
        else: 
            print "TEST DUMP HERE"
            c2f, c2stor = lamChi2Dump(MinObj.values['k'], MinObj.values['Beta'])

        

        
        #MinObj = M(lamChi2, k = 1.0, fix_k = True, Beta = 0.0, error_Beta = 0.1)
        

        MinObj.set_strategy(2)

        fmin, param = MinObj.migrad(nsplit= 10)
        #fmin, param = MinObj.migrad()
        #ErrDict = MinObj.minos()

        
        self.covar =  MinObj.np_covariance()

        ErrDict = MinObj.minos(maxcall = 1000)
        

        #plt.scatter(nData, c2stor)
        #plt.xlabel('nData')
        #plt.ylabel('chi2 in bin')
        #plt.savefig(self.realName + 'Chi2VsnData.png')
        #plt.clf()
        if self.nprint < simInd:
            print "Shapes of things"
            print len(c2stor)
            print nData.shape

            print dataBins.shape

            print self.binListFit.shape
            print self.binListSamp.shape
            #print DebugNIaPhot.shape
            #print DebugNCCPhot.shape
            #print DebugNIaTrue.shape
            #print DebugNCCTrue.shape

        for c2, nd in zip(c2stor, nData):
            self.globalChi2Storage.append(c2)
            self.globalNDataStorage.append(nd)

        if self.Rate_Model == 'discrete':
            fJList = [MinObj.values['f_1'],MinObj.values['f_2'],MinObj.values['f_3'],MinObj.values['f_4'],MinObj.values['f_5'],MinObj.values['f_6'],MinObj.values['f_7'],MinObj.values['f_8'],MinObj.values['f_9'],MinObj.values['f_10'],MinObj.values['f_11']]
            fJErrList = [MinObj.errors['f_1'],MinObj.errors['f_2'],MinObj.errors['f_3'],MinObj.errors['f_4'],MinObj.errors['f_5'],MinObj.errors['f_6'],MinObj.errors['f_7'],MinObj.errors['f_8'],MinObj.errors['f_9'],MinObj.errors['f_10'],MinObj.errors['f_11']]

            
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
            #kErr = MinObj.errors['k']
            kErr = (np.abs(ErrDict['k']['lower']) + np.abs(ErrDict['k']['upper']))/2.0
            Beta = MinObj.values['Beta']
            #BetaErr = MinObj.errors['Beta']
            BetaErr = (np.abs(ErrDict['Beta']['lower']) + np.abs(ErrDict['Beta']['upper']))/2.0
            if self.Rate_Model == 'brokenpowerlawVar':
                zBreak = MinObj.values['zBreak']
                zBreakErr = MinObj.values['zBreakErr']
            self.k = k
            self.Beta = Beta
            self.kErr = kErr
            self.BetaErr = BetaErr
            #/(self.nbins - 2)
            self.BetaRatio = (1+zCentersFit)**(Beta)
            self.fJList = None
            self.fracCCData = (NCC*1.0)/(1.0*(1.0*NCC + NIa*self.BetaRatio))
            self.fracCCDataTot = (np.sum(NCC)*1.0)/(1.0*(1.0*np.sum(NCC) + np.sum(NIa*self.BetaRatio)))


            #print self.fracCCDataTot
            #print type(self.fracCCDataTot)
            #assert(type(self.fracCCDataTot) == float)
            print "Chi2 final = {0}".format(round(lamChi2Dump(self.k, self.Beta)[0], 4))
        self.chi2 = fmin.fval
        print "Chi2final? = {0}".format(round(fmin.fval, 4))


        if not(self.priorRate is None):
            ratePriorFinalVal =  ratePrior(self.k*self.MCK, self.Beta+self.MCBeta, self.priorRate, self.priorZEff, self.ratePriorErrUp, self.ratePriorErrDown, self.ratePriorErrAll  )
            c2NoPrior = chi2func(nData, nSim, self.effmat, fnorm, zCentersSamp, zCentersFit, self.k, self.Beta, dump = False, nIA = NIa, nCC = NCC, simInd =simInd, TrueNCCData = TrueNCC)
            print "RATE PRIOR FINAL"
            print ratePriorFinalVal
            print "Chi2final? = {0}".format(round(fmin.fval, 4))
            print "Chi2FinalNoPrior"
            print c2NoPrior

        #fJs = np.ones(zCenters.shape)
        
        try:
            if (Rate_Model != 'discrete'):
                plt.clf()
                MinObj.draw_contour('k','Beta', nsigma=3)
                plt.savefig('{0}_{1}_k_beta_contour.png'.format(self.realName, self.simName))
                if Blind:
                    locs, labels = plt.xticks()
                    labels = locs + np.cos(cosVal)
                    plt.xticks(labels)
                    locs, labels = plt.yticks()
                    labels = locs + np.cos(cosVal)
                    plt.yticks(labels)
                plt.clf()
                
                #xgrid,ygrid, sigma, rawdata = MinObj.mncontour_grid('k', 'Beta', numpoints=30, sigma_res = 1, nsigma = 2.0)
                #fig, ax = plt.subplots(1)
                #plt.clf()
                #CS = ax.contour(xgrid, ygrid + self.MCBeta, sigma, levels = [ 1.0, 2.0])
                #ax.clabel(CS, fontsize=7, inline=1)
                #ax.set_xlabel('k')
                #ax.set_ylabel('Beta')
                #if Blind:
                #    ax.set_xticklabels([])
                #    ax.set_yticklabels([])
                #plt.savefig('{0}_{1}_k_beta_contour.png'.format(self.realName, self.simName))
                #plt.close()
        except: 
            print "Plot Fail A"

        try:
            if (Rate_Model != 'discrete'):
                plt.clf()
                MinObj.draw_profile('Beta', text = False)
                if Blind:

                    locs, labels = plt.xticks()
                    labels = locs + np.cos(cosVal)
                    plt.xticks(labels)
                plt.savefig('{0}_{1}_beta_contour.png'.format(self.realName, self.simName))
                plt.clf()
        except:
            print "Plot Fail C"
        try:
            if Rate_Model != 'discrete':
                Betas = np.linspace(self.Beta - 0.5, self.Beta + 0.5, 51)
                FCNs = []
                for bTemp in Betas:
                    FCN = lamChi2( self.k, bTemp)
                    FCNs.append(FCN)

                plt.plot(Betas, FCNs, c = 'k', label = 'Non Minuit Contour')
                plt.legend()
                plt.xlabel('Beta')
                plt.ylabel('Chi2')
                if Blind:

                    locs, labels = plt.xticks()
                    labels = locs + np.cos(cosVal)
                    plt.xticks(labels)
                plt.savefig('{0}_{1}_beta_mycontour.png'.format(self.realName, self.simName))
                plt.clf()


        except:
            print "Plot Fail D"

        if Rate_Model != 'discrete':
            plt.clf()
            ax = plt.axes()
            Betas = np.linspace(self.Beta - 0.1, self.Beta + 0.1, 501)
            FCNs = []
            for bTemp in Betas:
                FCN = lamChi2( self.k, bTemp)
                FCNs.append(FCN)

            plt.plot(Betas, FCNs, c = 'k', label = 'Non Minuit Contour')
            plt.legend()
            plt.xlabel('Beta')
            plt.ylabel('Chi2')
            if Blind:

                locs, labels = plt.xticks()
                labels = locs + np.cos(cosVal)
                ax.set_xticklabels(labels)
            print "FCNs"
            print FCNs
            plt.savefig('{0}_{1}_beta_myzoomcontour.png'.format(self.realName, self.simName))
            plt.clf()


            plt.clf()
            ax = plt.axes()
            ks = np.linspace(self.k - 0.1, self.k + 0.1, 501)
            FCNs = []
            for kTemp in ks:
                FCN = lamChi2(  kTemp,self.Beta)
                FCNs.append(FCN)

            plt.plot(ks, FCNs, c = 'k', label = 'Non Minuit Contour')
            plt.legend()
            plt.xlabel('k')
            plt.ylabel('Chi2')
            
            print "FCNs"
            print FCNs
            plt.savefig('{0}_{1}_k_myzoomcontour.png'.format(self.realName, self.simName))
            plt.clf()



            df = np.array(FCNs[1:]) - np.array(FCNs[:-1])
            inds = np.where(df > 0)[0]
            print 'inds'
            print inds
            print inds < 250
            print np.where(inds < 250)
            inds = inds[np.where(inds < 250)]
            print 'inds'
            print inds
            print "INDSSHAPE"
            print inds.shape
            if inds.shape[0]:
                print "MINUIT IS PROBABLY MAD. HERES WHY"
                print inds
                print Betas[inds]
                if inds.shape[0] > 1:
                    inds = inds[-1]
                print inds
                print Betas[inds]

                lamChi2Dump(self.k, Betas[inds -3])
                print "MINUIT MAD 2"
                lamChi2Dump(self.k, Betas[inds -2])
                print "MINUIT MAD 3"
                lamChi2Dump(self.k, Betas[inds -1])

                print "MINUIT MAD 4"
                lamChi2Dump(self.k, Betas[inds])
                print "MINUIT MAD 5"
                lamChi2Dump(self.k, Betas[inds + 1])
                print "MINUIT MAD 6"
                lamChi2Dump(self.k, Betas[inds + 2])
                print "MINUIT MAD 7"
                lamChi2Dump(self.k, Betas[inds + 3])
                print "END MINUIT MAD"
        



        try:
            if (Rate_Model != 'discrete'):
                plt.clf()
                MinObj.draw_mncontour('k','Beta', nsigma=3)
                plt.savefig('{0}_{1}_k_beta_mncontour.png'.format(self.realName, self.simName))
                if Blind:
                    locs, labels = plt.xticks()
                    labels = locs + np.cos(cosVal)
                    plt.xticks(labels)
                    locs, labels = plt.yticks()
                    labels = locs + np.cos(cosVal)
                    plt.yticks(labels)
                plt.clf()
                MinObj.draw_mnprofile('Beta', text = False, subtract_min = True)
                if Blind:
                    

                    locs, labels = plt.xticks()
                    labels = locs + np.cos(cosVal)
                    plt.xticks(labels)
                plt.savefig('{0}_{1}_beta_mncontour.png'.format(self.realName, self.simName))
                plt.clf()
                #xgrid,ygrid, sigma, rawdata = MinObj.mncontour_grid('k', 'Beta', numpoints=30, sigma_res = 1, nsigma = 2.0)
                #fig, ax = plt.subplots(1)
                #plt.clf()
                #CS = ax.contour(xgrid, ygrid + self.MCBeta, sigma, levels = [ 1.0, 2.0])
                #ax.clabel(CS, fontsize=7, inline=1)
                #ax.set_xlabel('k')
                #ax.set_ylabel('Beta')
                #if Blind:
                #    ax.set_xticklabels([])
                #    ax.set_yticklabels([])
                #plt.savefig('{0}_{1}_k_beta_contour.png'.format(self.realName, self.simName))
                #plt.close()
        except: 
            print "Plot Fail B"
            pass
        
        

    
        #plt.axhline(y = self.MCBeta, c = 'k', label = 'True Beta')
        #plt.axhline(y = Beta + self.MCBeta, c = 'g', label= 'Best Fit Beta')
        #plt.axvline(x = k, label = 'Best Fit k')
        
    '''
    def chi2V2(self, fJs, fJErrs, zCenters, k, Beta):
            fitfJs = k*(1+zCenters)**Beta
            Chi2Temp = 0
            for fJ, fitfJ, fJErr in zip(fJs, fitfJs, fJErrs):
                Chi2Temp += (fJ - fitfJ)**2/(fJ + fJErr)
            return Chi2Temp
    '''

def weakPrior(value, priorTuple):
    if value < priorTuple[1]:
        if value > priorTuple[0]:
            return 1
        else: 
            return (value - priorTuple[0])**4
    else:
        return (value - priorTuple[1])**4

def ratePrior(fitK, fitBeta, priorRate, zEffPrior, priorRateErrUp = None, priorRateErrDown = None, priorRateErrAll = None):

    print "PRIOR"
    print priorRate
    print zEffPrior
    print priorRateErrUp
    print priorRateErrDown
    print "Fit Beta/k"
    print fitBeta
    print fitK
    fitRate = fitK*(1+zEffPrior)**fitBeta
    print 'Fit Rate'
    print fitRate
    print "PriorChi2"

    if fitRate > priorRate:

        if not (priorRateErrUp is None):
            print (fitRate - priorRate)**2/priorRateErrUp**2
            return (fitRate - priorRate)**2/priorRateErrUp**2
        else:
            print (fitRate - priorRate)**2/priorRateErrAll**2
            return (fitRate - priorRate)**2/priorRateErrAll**2
    else:
        if not (priorRateErrDown is None):
            print (fitRate - priorRate)**2/priorRateErrDown**2
            return (fitRate - priorRate)**2/priorRateErrDown**2
        else:
            print (fitRate - priorRate)**2/priorRateErrAll**2
            return (fitRate - priorRate)**2/priorRateErrAll**2







def getCCScale(simCat, dataCat, MURESWindow = (-1, 1), zbins = [0.0, 0.3, 0.6, 0.9, 1.2], Beta = None, binList = None, fracCCData = None, outfilePrefix = 'Test', Rate_Model = 'powerlaw', f_Js = None, returnHist = False, debug = False, simInd = 100, ztype = 'zPHOT'):
    #import iminuit as iM
    #from iminuit import Minuit as M
    if debug:
        print "Check this"
        print Rate_Model
        print f_Js
        print Beta
        print fracCCData
        print "Done Checking"
    CCScales = []
    CCScaleErrs = []
    simIaHists = []
    simCCHists = []
    dataHists = []
    if not(f_Js is None):
        f_Js = np.array(f_Js)

    allSimCC = simCat[simCat['SIM_TYPE_INDEX'].astype(int) != 1]
    allSimIa = simCat[simCat['SIM_TYPE_INDEX'].astype(int) == 1]
    allData = np.copy(dataCat)


    #fnorm2 = float(dataCat.shape[0])/float(np.sum(simHist))
    
    simCat = simCat[(simCat['MURES'] < MURESWindow[0]) | (simCat['MURES'] > MURESWindow[1]) ]
    dataCat = dataCat[(dataCat['MURES'] < MURESWindow[0]) | (dataCat['MURES'] > MURESWindow[1]) ]
    

    for zl, zh in zip(zbins[:-1], zbins[1:]):

        tempSim = simCat[(simCat[ztype] < zh) & (simCat[ztype] > zl)]
        tempData = dataCat[(dataCat[ztype] < zh) & (dataCat[ztype] > zl)]


        allSimCCZbin = allSimCC[(allSimCC[ztype] < zh) & (allSimCC[ztype] > zl)]
        allSimIaZbin = allSimIa[(allSimIa[ztype] < zh) & (allSimIa[ztype] > zl)]
        if debug:
            print "all Sim CC Zbin/IaZbin"
            print allSimCCZbin.shape[0]
            print allSimIaZbin.shape[0]

        allDataZbin = allData[(allData[ztype] < zh) & (allData[ztype] > zl)]



        tempSimCC = tempSim[tempSim['SIM_TYPE_INDEX'] != 1]
        tempSimIa = tempSim[tempSim['SIM_TYPE_INDEX'] == 1]

        R = float(tempData.shape[0])/float(allDataZbin.shape[0])
        if debug:
            print "R"

            print R

            print "Hist CC, outlier and total"
            print tempSim.shape[0]
            print allSimCCZbin.shape[0]


            print "pre Beta Correction allSimIa"
            print tempData.shape[0]
            print allSimIaZbin.shape[0]

        if Rate_Model == 'discrete':
            hist, bins = np.histogram(allSimIaZbin[ztype], bins = 11)
            if debug:
                print 'fJ shape'
                print f_Js.shape
                print f_Js
                print hist
                print bins
            betaCorrAllSimIaZbin =np.sum(hist*f_Js)
        else:
            betaCorrAllSimIaZbin =np.sum((1+ allSimIaZbin[ztype])**Beta)
        #S = float(np.array(R*histSAllIa) - np.array(tempSimIa.shape[0]))/float(np.array(tempSimCC.shape[0]) - np.array(R*histSAllCC))

        try:
            if debug:
                print "Test S"
                print R
                print betaCorrAllSimIaZbin
                print tempSimIa.shape[0]
                print tempSimCC.shape[0]
                print allSimCCZbin.shape
                print 'EEE'
                print np.array(R*betaCorrAllSimIaZbin)
                print 'DDD'
                print np.array(tempSimIa.shape[0])
                print 'CCC'
                print (np.array(tempSimCC.shape[0]) - np.array(R*allSimCCZbin.shape[0]))
                print "AAA"
                print (np.array(R*betaCorrAllSimIaZbin) - np.array(tempSimIa.shape[0]))/(np.array(tempSimCC.shape[0]) - np.array(R*allSimCCZbin.shape[0]))
                print "BBB"
                #S = (np.array(R*betaCorrAllSimIaZbin) - np.array(tempSimIa.shape[0]))/(np.array(tempSimCC.shape[0]) - np.array(R*allSimCCZbin.shape[0]))
            S = float(np.array(R*betaCorrAllSimIaZbin) - np.array(tempSimIa.shape[0]))/float(np.array(tempSimCC.shape[0]) - np.array(R*allSimCCZbin.shape[0]))
        except: 
            S = np.nan
        if debug:
            print "S WTF"
            print S


            print "Uncertainty Related Bullshit"
            '''
            print "Delta R"

            dR = np.sqrt(histD + histDAll)

            print dR

            num1 = np.sqrt(np.sqrt((dR/R)**2 + histSAllIa) + tempSimIa.shape[0])

            num2 = np.sqrt(np.sqrt((dR/R)**2 + histSAllCC) + tempSimCC.shape[0])

            den1 = (R*histSAllIa - tempSimIa.shape[0])

            den2 = (tempSimCC.shape[0] - R*histSAllCC)


            dS = np.sqrt((num1/den1)**2 + (num2/den2)**2)
            '''
        #ddnCC = np.sqrt(tempSimCC.shape[0])*(tempSimIa.shape[0] - histSAllIa*R)/(tempSimCC.shape[0] - R*histSAllCC)**2

        #ddNCC = np.sqrt(histSAllCC)*R*(histSAllIa*R - tempSimIa.shape[0])/(tempSimCC.shape[0] - R*histSAllCC)**2

        #ddnIa = np.sqrt(tempSimIa.shape[0])/(tempSimCC.shape[0] - R*histSAllCC)
        #ddNIa = np.sqrt(histSAllIa)*R/(tempSimCC.shape[0] - R*histSAllCC)

        ddnCC = np.sqrt(tempSimCC.shape[0])*(tempSimIa.shape[0] - allSimIaZbin.shape[0]*R)/(tempSimCC.shape[0] - R*allSimCCZbin.shape[0])**2

        ddNCC = np.sqrt(allSimCCZbin.shape[0])*R*(allSimIaZbin.shape[0]*R - tempSimIa.shape[0])/(tempSimCC.shape[0] - R*allSimCCZbin.shape[0])**2

        ddnIa = np.sqrt(tempSimIa.shape[0])/(tempSimCC.shape[0] - R*allSimCCZbin.shape[0])
        ddNIa = np.sqrt(allSimIaZbin.shape[0])*R/(tempSimCC.shape[0] - R*allSimCCZbin.shape[0])

        #ddR = (histSAllIa*tempSimCC.shape[0] - histSAllCC * tempSimIa.shape[0])/(tempSimCC.shape[0] - R*histSAllCC)**2

        dS = np.sqrt(ddnCC**2 + ddNCC**2 + ddnIa**2 + ddNIa**2)# + ddR**2)

        if debug:

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
        if np.isnan(S):
            print 'SCALE IS NAN'
            if len(CCScales) > 0:
                #CCScales.append(CCScales[-1])
                CCScales.append(1.0)
            else: 
                CCScales.append(1.0)
        else:
            CCScales.append(S)
        if type(dS) == np.ndarray:
            if np.isnan(dS[0]):
                CCScaleErrs.append(1.0)
            else:
                CCScaleErrs.append(dS[0])
        else:
            if np.isnan(dS):
                CCScaleErrs.append(1.0)
            else:
                CCScaleErrs.append(dS)

        #if debug:
        #    print "CC PlotDebug"
        #    print (simBinsCC[1:] + simBinsCC[:-1])/2.0
        #    print simHistCC
        #   print CCScales[0]
        #    print dS
        #    print fnorm2
        #    print histD
        #    print (muresBins[1:] + muresBins[:-1])/2.0
   
        #if simInd ==1:
        #    plt.step((simBinsCC[1:] + simBinsCC[:-1])/2.0, simHistCC*fnorm2, c = 'b', where = 'mid',  label = 'prescaled Sim CC')
        #    plt.step((simBinsCC[1:] + simBinsCC[:-1])/2.0, CCScales[0]*simHistCC*fnorm2, c = 'g', where = 'post',  label = 'postscaledSimCC')
        #    plt.step((muresBins[1:] + muresBins[:-1])/2.0, histD, c = 'r', where = 'mid',  label = 'data')
        #    plt.legend()
        #    plt.savefig(outfilePrefix + 'ScaledHist.png')
        #    plt.clf()
    if debug:
        print "CCScaleErrs"
        print CCScaleErrs
    if returnHist:
        return CCScales, CCScaleErrs, simIaHists, simCCHists, dataHists
    return CCScales, CCScaleErrs

def applyCCScale(NCC, CCScales, CCScaleErrs, datazbins = None,  zbins = None):
    if not(zbins is None):
        zbins = np.array(zbins)
    if not (datazbins is None):
        datazbins = np.array(datazbins)
    if type(CCScaleErrs) == list:
        CCScaleErrs = np.array(CCScaleErrs)
    if type(CCScales) == list:
        CCScales = np.array(CCScales)
    print 'CCScaleErrs'
    print CCScaleErrs
    print datazbins
    print zbins


    
    if type(CCScales) == np.ndarray:
        if CCScales.shape[0] == 1:
            NCCScaled = CCScales[0]*NCC
        else:
            if (datazbins is None) | (zbins is None):
                assert(0)
            if CCScales.shape[0] < 4:
                k = CCScales.shape[0] -1
            else:
                k = 3
            
            nancond = np.isnan(CCScales)
            if np.sum(nancond) > 0:
                CCScales[nancond] = 1.
                CCScaleErrs[nancond] = 1.

            zCenters = (zbins[1:]+ zbins[:-1])/2.0
            print zCenters
            print CCScales
            
            #spline = UnivariateSpline(zbins, CCScales, w = 1.0/CCScaleErrs, k = k)
            spline = UnivariateSpline(zCenters, CCScales, w = 1.0/CCScaleErrs, k = k)

            print datazbins.shape
            print datazbins
            print NCC.shape

            datazcents = (datazbins[1:]+ datazbins[:-1])/2.0

            NCCScaled = spline(datazcents)*NCC

    elif (type(CCScales) == int) | (type(CCScales) == float):
        NCCScaled = CCScales*NCC
    else:
        assert(0)

    NCCScaled = NCCScaled.clip(0)
    print NCCScaled

    assert(not bool(np.sum(NCCScaled < 0)))


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
    cutFiles = [argv[10]]
    try:
        debug = bool(int(argv[11]))
    except:
        debug = False

    
    #if( ('Combine' in simdir) or ('SALT2' in simdir)) &  (('Combine' in datadir) or ('SALT2' in simdir)):
        #NNCut = True
        #NNProbCut = 0.95
    
    #if len(argv) > 6:
    #    NNCut = True
    #    NNProbCut = 0.9
    #    NNData = argv[6]
    #    NNSim = argv[7]

    
    #default params

    zminFit = 0.1
    zmaxFit = 1.2
    zminSamp = 0.1
    zmaxSamp = 1.2
    MJDMin = 0.0
    MJDMax = np.inf
    bins = "equalSize" 
    runFit = True
    fracContamCuts = [-1]
    fixBeta = True
    fixK = False
    nbins = None
    binList = None
    ScaleMuResCutLow = -1
    ScaleMuResCutHigh = 1
    #muresBins = 1
    muresBinsLow  = 3
    muresBinsHigh = 3
    scaleZBins = [0.0, 1.2]
    nScaleZBins = None
    cheatCCSub = False
    cheatCCScale = False
    ZSysFlag = False
    Blind = False
    Rate_Model = 'powerlaw'
    MURESCuts = 2.0 #[(0.0, 0.8, -0.5, 0.5), (0.8, 1.5, -1, 1)]
    noCCMC = False
    fixCCScale = False
    trueMCBeta = 1.65
    trueMCK = 1.97E-5

    priorRate = None
    priorZEff = None
    ratePriorErrUp = None
    ratePriorErrDown =None
    ratePriorErrAll = None
    priors = None

    #override file

    params = open(paramFile, 'r').readlines()

    for p in params:

        print p
        exec(p)

    if nScaleZBins is None :
        redoScaleZBinFlag = False

    else:
        redoScaleZBinFlag = True

    if not(priors is None):
        if len(priors) == 3:
            priorRate, priorZEff, ratePriorErrAll = priors
            ratePriorErrUp = None
            ratePriorErrDown = None
        elif len(priors) == 4:
            priorRate, priorZEff, ratePriorErrUp, ratePriorErrDown = priors
            ratePriorErrAll =None




    cosVal = 47392945716038.134971247
    kmean = []
    ksigma = []
    kErr = []
    BetaMean = []
    #BetaWeightMean = []
    #KWeightMean = []
    BetaSigma= []
    BetaErr = []
    zBreakMeans = []
    zBreakSigmas =[]
    zBreakErrs = []
    Chi2Mean = []
    Chi2Sigma = []
    f_JStorage = []
    f_JErrStorage = []
    SampleSizes = []

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
        zBreaks =[]
        zBreakErrs = []
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
            nfile = 48
        else:
            nfile = 2
        for simInd in range(1,nfile):
            

            #print "Sim {0}".format(simInd)
            #SimBeta = 2.1 # simdir.split('_')[-3]
            #SimR0 = 1.7*10**-5 #simdir.split('_')[-5]
            #print "Sim R0 = {1}; Sim Beta = {0}".format(SimBeta, SimR0)

            
            print datadir.format(simInd)
            if simLoaded:
                try:
                
                    RateTest.newData(datadir.format(simInd), dataname.format(simInd), simInd =simInd)
                    if ZSysFlag:
                        assert(0)
                        RateTest.zSystematic(nbins = nbins, binList = binList)


                    if redoScaleZBinFlag:

                        RealCat = RateTest.postCutRealCat 
                        RealOutlierCat = RealCat[(RealCat['MURES'] > muresBinsHigh)| (RealCat['MURES'] < muresBinsLow)]

                        zArray =RealOutlierCat[RateTest.ztype]
                        zArray.sort()

                        splitZs = np.array_split(zArray, nScaleZBins)

                        #[(0[0], (0[-1] + 1[0]), (1[-1] + 2[0]), 2[1]]

                        scaleZBins = [splitZs[0][0]]

                        
                        for i in range(1,nScaleZBins):

                            scaleZBins.append((splitZs[i-1][-1] + splitZs[i][0] )/2.0)
                        scaleZBins.append(splitZs[i][-1])


                    #RateTest.effCalc(nbins = nbins, fracContamCut = fcc, simInd =simInd)
                    #RateTest.effCalc(nbins = 20)
                    BetaIter = []
                    BetaErrIter = []
                    CCIter = []
                    CCErrIter = []
                    RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, simInd =simInd, trueBeta = trueBeta - trueMCBeta, CCScale = 1.0, TrueCCScale = TrueCCScale, scaleZBins = scaleZBins, Blind = Blind)
                    if Rate_Model != 'discrete':
                        if Blind:
                            print "Blinding A"
                            BetaIter.append(RateTest.Beta+ np.cos(cosVal))
                        else:
                            BetaIter.append(RateTest.Beta)
                        BetaErrIter.append(RateTest.BetaErr)

                    for iteration in range(nIter):
                        if not fixCCScale:
                            if not noCCMC:
                                CCScale, CCScaleErr =  getCCScale(RateTest.postCutSimCat, RateTest.postCutRealCat, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, Beta = RateTest.Beta, binList = RateTest.binListFit, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname,Rate_Model = Rate_Model, f_Js =RateTest.fJList, simInd = simInd, debug = debug, ztype = RateTest.ztype)
                                CCIter.append(CCScale)
                                CCErrIter.append(CCScaleErr)
                                RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - trueMCBeta, CCScale = CCScale, CCScaleErr = CCScaleErr, TrueCCScale = TrueCCScale, BetaInit = RateTest.Beta, kInit = RateTest.k, BetaErr = RateTest.BetaErr, kErr = RateTest.kErr, f_Js =RateTest.fJList, CCZbins = scaleZBins , scaleZBins = scaleZBins, Blind = Blind)
                            else:
                                CCIter.append(0.0)
                                CCErrIter.append(0.0)
                                RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - trueMCBeta, CCScale = 0.0, CCScaleErr = 1.0, TrueCCScale = 0.0, BetaInit = RateTest.Beta, kInit = RateTest.k, BetaErr = RateTest.BetaErr, kErr = RateTest.kErr, f_Js =RateTest.fJList, CCZbins = scaleZBins , scaleZBins = scaleZBins, Blind = Blind)
                        else:
                            CCIter.append(1.0)
                            CCErrIter.append(0.0)
                            RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - trueMCBeta, CCScale = 1.0, CCScaleErr = 1.0, TrueCCScale = 0.0, BetaInit = RateTest.Beta, kInit = RateTest.k, BetaErr = RateTest.BetaErr, kErr = RateTest.kErr, f_Js =RateTest.fJList, CCZbins = scaleZBins , scaleZBins = scaleZBins, Blind = Blind)

                        if Blind:
                            print "Blinding b"
                            BetaIter.append(RateTest.Beta+ np.cos(cosVal))
                        else:
                            BetaIter.append(RateTest.Beta)
                        BetaErrIter.append(RateTest.BetaErr)
                    if not fixCCScale:
                        if not noCCMC:
                            CCScale, CCScaleErr =  getCCScale(RateTest.postCutSimCat, RateTest.postCutRealCat, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, Beta = RateTest.Beta, binList = RateTest.binListFit, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname,Rate_Model = Rate_Model, f_Js =RateTest.fJList, simInd = simInd, debug = debug, ztype = RateTest.ztype)
                            CCIter.append(CCScale)
                            CCErrIter.append(CCScaleErr)
                    else:
                        CCIter.append(1.0)
                        CCErrIter.append(0.0)
                    
                    print "CCScale Progression"
                    print CCIter
                    print "CCScale Err Progression"
                    print CCErrIter
                    if Rate_Model != 'discrete':
                        print "Beta Progression"
                        print BetaIter
                        print "Beta Err Progressions"
                        print BetaErrIter
                        print "Mean Betas"
                        print np.nanmean(BetaIter)

                        print "Mean CCScales"
                        print np.nanmean(CCIter)
                    else:
                        f_JStorage.append(RateTest.fJList)
                        f_JErrStorage.append(RateTest.fJErrList)

                    #print "AAA CC Scales"
                    if not fixCCScale:

                        if not noCCMC:
                            CCScale, CCScaleErr =  getCCScale(RateTest.postCutSimCat, RateTest.postCutRealCat, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, Beta = RateTest.Beta, binList = RateTest.binListFit, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname, Rate_Model = Rate_Model, f_Js =RateTest.fJList, simInd = simInd, debug = debug, ztype = RateTest.ztype)
                            print CCScale
                            CCScaleStorage.append(CCScale)
                            CCScaleErrStorage.append(CCScaleErr)
                        else:
                            CCScaleStorage.append(0.0)
                            CCScaleErrStorage.append(1.0)
                    else:
                        CCScaleStorage.append(1.0)
                        CCScaleErrStorage.append(1.0)
                    


                    ks.append(RateTest.k)
                    kErrs.append(RateTest.kErr)
                    if Blind:
                        print "Blinding c"
                        Betas.append(RateTest.Beta+ np.cos(cosVal))

                    else:
                        Betas.append(RateTest.Beta)
                    BetaErrs.append(RateTest.BetaErr)
                    if Rate_Model == 'brokenpowerlawVar':
                        zBreaks.append(Rate_Fitter.zBreak)
                        zBreakErrs.append(Rate_Fitter.zBreakErr)

                    Chi2s.append(RateTest.chi2)
                    print "CCScale Storage Iter {0}".format(simInd)
                    print CCScaleStorage
                    if not noCCMC:
                        print CCScale
                        print CCScale[0]

                    
                    dnamestr = datadir.format(simInd)

                    cutdnamestr = dnamestr.split('.')[0] + '+CUTS.FITRES.gz'
                    #if saveCuts:
                    #    np.savetxt(cutdnamestr, RateTest.realcat.Catalog, delimiter = ' ', fmt='%s')

                    lowzCut = zminFit
                    highzCut = zmaxFit
                    SampleSizes.append(  RateTest.realcat.Catalog[(RateTest.realcat.Catalog[RateTest.ztype] < zmaxFit) & (RateTest.realcat.Catalog[RateTest.ztype] > zminFit)].shape[0])
                    if saveCuts:
                        np.savetxt(cutdnamestr, RateTest.realcat.Catalog[(RateTest.realcat.Catalog[RateTest.ztype] < zmaxFit) & (RateTest.realcat.Catalog[RateTest.ztype] > zminFit)], delimiter = ' ', fmt='%s')
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

                    RateTest = Rate_Fitter(datadir.format(simInd), dataname.format(simInd), simdir, simname,simgenfile, trueMCBeta, trueMCK, zminSamp =zminSamp, zmaxSamp =zmaxSamp, zminFit =zminFit, zmaxFit =zmaxFit, cheatZ = cheatZ, cheatType = cheatType, cuts = cuts, cheatCCSub = cheatCCSub, cheatCCScale = cheatCCScale, Rate_Model = Rate_Model, MURESCuts = MURESCuts, noCCMC = noCCMC, priorRate = priorRate, priorZEff = priorZEff, ratePriorErrUp = ratePriorErrUp, ratePriorErrDown =ratePriorErrDown, ratePriorErrAll = ratePriorErrAll)# , MJDMin = 0, MJDMax = np.inf)
                    
                    if ZSysFlag:
                            RateTest.zSystematic(nbins = nbins, binList = binList)
                    simLoaded = True

                    RateTest.effCalc(nbinsSamp = nbinsSamp,nbinsFit = nbinsFit,  fracContamCut = fcc)
                    #RateTest.effCalc(nbins = 20)
                    BetaIter = []
                    BetaErrIter = []
                    CCIter = []
                    CCErrIter = []

                    if redoScaleZBinFlag:

                        RealCat = RateTest.postCutRealCat 
                        RealOutlierCat = RealCat[(RealCat['MURES'] > muresBinsHigh)| (RealCat['MURES'] < muresBinsLow)]

                        zArray =RealOutlierCat[RateTest.ztype]
                        zArray.sort()

                        print 'zArray'
                        print zArray
                        print 'nScaleZBins'
                        print nScaleZBins

                        splitZs = np.array_split(zArray, nScaleZBins)

                        #[(0[0], (0[-1] + 1[0]), (1[-1] + 2[0]), 2[1]]

                        scaleZBins = [splitZs[0][0]]

                        
                        for i in range(1,nScaleZBins):

                            scaleZBins.append((splitZs[i-1][-1] + splitZs[i][0] )/2.0)
                        scaleZBins.append(splitZs[i][-1])


                    RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, simInd =simInd, trueBeta = trueBeta - trueMCBeta, CCScale = 1.0, TrueCCScale = TrueCCScale, scaleZBins = scaleZBins, Blind = Blind)
                    if Rate_Model != 'discrete':
                        if Blind:
                            print "Blinding d"
                            BetaIter.append(RateTest.Beta+ np.cos(cosVal))
                        else:
                            BetaIter.append(RateTest.Beta)
                        BetaErrIter.append(RateTest.BetaErr)
                    for iteration in range(nIter):
                        print "interation Number"
                        print iteration
                        if not fixCCScale:
                            if not noCCMC:
                                CCScale, CCScaleErr =  getCCScale(RateTest.postCutSimCat, RateTest.postCutRealCat, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, Beta = RateTest.Beta, binList = RateTest.binListFit, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname, Rate_Model = Rate_Model, f_Js =RateTest.fJList, simInd = simInd, debug = debug, ztype = RateTest.ztype)
                                CCIter.append(CCScale)
                                CCErrIter.append(CCScaleErr)
                                RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - trueMCBeta, CCScale = CCScale, CCScaleErr = CCScaleErr, TrueCCScale = TrueCCScale, BetaInit = RateTest.Beta, kInit = RateTest.k, BetaErr = RateTest.BetaErr, kErr = RateTest.kErr, CCZbins = scaleZBins, scaleZBins = scaleZBins, Blind = Blind)
                            else:
                                CCIter.append(0.0)
                                CCErrIter.append(1.0)

                                RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - trueMCBeta, CCScale = 0.0, CCScaleErr = 1.0, TrueCCScale = 0.0, BetaInit = RateTest.Beta, kInit = RateTest.k, BetaErr = RateTest.BetaErr, kErr = RateTest.kErr, CCZbins = scaleZBins, scaleZBins = scaleZBins, Blind = Blind)
                        else:
                            CCIter.append(1.0)
                            CCErrIter.append(1.0)
                            RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - trueMCBeta, CCScale = 1.0, CCScaleErr = 1.0, TrueCCScale = 0.0, BetaInit = RateTest.Beta, kInit = RateTest.k, BetaErr = RateTest.BetaErr, kErr = RateTest.kErr, CCZbins = scaleZBins, scaleZBins = scaleZBins, Blind = Blind)

                        
                        if Rate_Model != 'discrete':
                            if Blind:
                                print "Blinding e"
                                BetaIter.append(RateTest.Beta+ np.cos(cosVal))
                            else:
                                BetaIter.append(RateTest.Beta)
                            BetaErrIter.append(RateTest.BetaErr)
                    if not fixCCScale:
                        if not noCCMC:
                            CCScale, CCScaleErr =  getCCScale(RateTest.postCutSimCat, RateTest.postCutRealCat, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, Beta = RateTest.Beta, binList = RateTest.binListFit, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname, Rate_Model = Rate_Model, f_Js =RateTest.fJList, simInd = simInd, debug = debug, ztype = RateTest.ztype)
                            CCIter.append(CCScale)
                            CCErrIter.append(CCScaleErr)
                    if Rate_Model != 'discrete':
                        print "Beta Progression"
                        print BetaIter
                        print "Beta Err Progressions"
                        print BetaErrIter
                       
                        print "Mean Betas"
                        print np.nanmean(BetaIter)

                    else:
                        f_JStorage.append(RateTest.fJList)
                        f_JErrStorage.append(RateTest.fJErrList)
                        
                    print "CCScale Progression"
                    print CCIter
                    print "CCScale Err Progression"
                    print CCErrIter
                    print "Mean CCScales"
                    print np.nanmean(CCIter)
                    if not fixCCScale:
                        if not noCCMC:
                            print "AAA CC Scales"
                            
                            CCScale, CCScaleErr =  getCCScale(RateTest.postCutSimCat, RateTest.postCutRealCat, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, Beta = RateTest.Beta, binList = RateTest.binListFit, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname, f_Js =RateTest.fJList, Rate_Model = Rate_Model, simInd = simInd, debug = debug, ztype = RateTest.ztype)
                            print 'CC Scale'
                            print CCScale
                            CCScaleStorage.append(CCScale)
                            CCScaleErrStorage.append(CCScaleErr)
                        else: 
                            CCScaleStorage.append(0.0)
                            CCScaleErrStorage.append(1.0)
                    else:
                        CCScaleStorage.append(1.0)
                        CCScaleErrStorage.append(1.0)

                    dnamestr = datadir.format(simInd)

                    cutdnamestr = dnamestr.split('.')[0] + '+CUTS.FITRES.gz'

                    np.savetxt(cutdnamestr, RateTest.realcat.Catalog, delimiter = ' ', fmt='%s')

                    #with open(cutdnamestr, 'rb') as f_in:
                    #    with gzip.open(cutdnamestr + '.gz', 'wb') as f_out:
                    #        shutil.copyfileobj(f_in, f_out)



                    cutsnamestr = simname.split('.')[0] + '+CUTS.FITRES.gz'

                    np.savetxt(cutsnamestr, RateTest.realcat.Catalog[(RateTest.realcat.Catalog[RateTest.ztype] < zmaxFit) & (RateTest.realcat.Catalog[RateTest.ztype] > zminFit)], delimiter = ' ', fmt = '%s')

                    lowzCut = zminFit
                    highzCut = zmaxFit
                    SampleSizes.append(  RateTest.realcat.Catalog[(RateTest.realcat.Catalog[RateTest.ztype] < zmaxFit) & (RateTest.realcat.Catalog[RateTest.ztype] > zminFit)].shape[0])

                    #with open(cutsnamestr, 'rb') as f_in:
                    #    with gzip.open(cutsnamestr + '.gz', 'wb') as f_out:
                    #        shutil.copyfileobj(f_in, f_out)


                    ks.append(RateTest.k)
                    kErrs.append(RateTest.kErr)
                    if Rate_Model != 'discrete':
                        if Blind:
                            print "Blinding f"
                            Betas.append(RateTest.Beta+ np.cos(cosVal))
                        else:
                            Betas.append(RateTest.Beta)
                        BetaErrs.append(RateTest.BetaErr)

                    if Rate_Model == 'brokenpowerlawVar':
                        zBreaks.append(Rate_Fitter.zBreak)
                        zBreakErrs.append(Rate_Fitter.zBreakErr)

                    Chi2s.append(RateTest.chi2)
                    print "CCScale Storage Iter {0}".format(simInd)
                    print CCScaleStorage
                    if not noCCMC:
                        print CCScale
                        print CCScale[0]
                    if Rate_Model != 'discrete':
                        if np.isnan(RateTest.Beta):
                            nFail +=1

                except Exception, e:
                    print "FAILURE"
                    print e
                    traceback.print_exc()
                    nFail +=1
        #if Blind:
        #    Betas = np.array(Betas) + np.cos(47392945716038.134971247)
        print "Number of Failures"
        print nFail
        if Rate_Model != 'discrete':

            badSims = np.invert(np.isfinite(Betas) & (BetaErrs > 0) & np.isfinite(ks) & (kErrs > 0))
            mBetas = ma.masked_array(Betas, mask=badSims)
            mBetaErrs = ma.masked_array(BetaErrs, mask=badSims)
            mks = ma.masked_array(ks, mask=badSims)
            mkErrs = ma.masked_array(kErrs, mask=badSims)
            print "mean k"
            print np.nanmean(ks)
            print "mean kerrs"
            print np.nanmean(kErrs)
            print "std. k"
            print np.nanstd(ks)
            print "Mean beta"
            print np.nanmean(Betas)
            print "Mean betaerrs"
            print np.nanmean(BetaErrs)
            print "std. beta"
            print np.nanstd(Betas)
            if len(Betas) == 1:
                kmean.append(ks[0])
                ksigma.append(0.0)
                kErr.append(kErrs[0])
                BetaMean.append(Betas[0])
                BetaSigma.append(0.0)
                BetaErr.append(BetaErrs[0])
            else:
                print "test here"
                print ks
                print mks
                print Betas
                print mBetas
                print 'end test here'
                kmean.append(np.average(mks, weights = 1.0/mkErrs**2))
                ksigma.append(np.std(mks))
                kErr.append(np.mean(mkErrs))
                BetaMean.append(np.average(mBetas, weights = 1.0/mBetaErrs**2))
                #BetaWeightMean.append(np.average(Betas, weights = 1.0/ma.masked_invalid(BetaErrs)**2))
                #KWeightMean.append(np.average(ks, weights = 1.0/ma.masked_invalid(kErrs)**2))
                BetaSigma.append(np.std(mBetas))
                BetaErr.append(np.mean(mBetaErrs))
        else:
            print "mean f_Js"
            print np.nanmean(f_JStorage, axis =0)
            print "mean f_JErrs"
            print np.nanmean(f_JErrStorage, axis =0)
        if Rate_Model == 'brokenpowerlawVar':
            zBreakMeans.append(np.nanmean(zBreaks))
            zBreakSigmas.append(np.nanstd(zBreaks))

        Chi2Mean.append(np.nanmean(Chi2s))
        Chi2Sigma.append(np.nanstd(Chi2s))

        

        
        if simInd == 1:
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

        if not noCCMC:
            print "AAA CC Scale means (weighted, unweighted)"
            #print np.average(ma.masked_invalid(np.array(CCScaleStorage)),weights = 1.0/ma.masked_invalid(CCScaleErrStorage)**2, axis = 0)
            #print np.nanmean(ma.masked_invalid(np.array(CCScaleStorage)), axis = 0)
            #print CCScaleStorage
            #print CCScaleErrStorage
            print np.average(np.array(CCScaleStorage),weights = 1.0/ma.masked_invalid(CCScaleErrStorage)**2, axis = 0)
            print np.nanmean(np.array(CCScaleStorage), axis = 0)
            print "AAA CC Scale stds"
            print np.nanstd(np.array(CCScaleStorage), axis = 0)
            CCScaleStorageGlobal.append(CCScaleStorage)

        

    print "All Betas"
    print Betas

    if cheatType:
        print "THESE RESULTS ONLY INCLUDE TRUE Ias BECAUSE WE CHEATED AND USED THE SIM INFORMATION"
    if cheatZ:
        print "THESE RESULTS Use Simulated Redshift info"
    '''
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

    '''
    #print "MURES CUTS"
    #print MURES_Cuts
    print "Frac Contam Cuts"
    print fracContamCuts
    if Rate_Model != 'discrete':
        print "Kmeans"
        print kmean
        print "Ksigmas"
        print ksigma
        print "BetaMeans"
        print BetaMean
        print "BetaSigmas"
        print BetaSigma
        print "BetaErrs"
        print BetaErr
    else: 
        print "f_J mean unweighted"
        print np.mean(f_JStorage, axis = 0)
        print "f_J mean weighted"
        print np.average(f_JStorage, weights = 1.0/(np.array(f_JErrStorage))**2, axis = 0)

        print "f_J Errors"
        print np.mean(f_JErrStorage, axis = 0)

    if Rate_Model == 'brokenpowerlawVar':
        print "mean powerlaw break z"
        print zBreakMeans
        print "st. dev powerlaw break z"
        print zBreakSigmas
    print "Chi2Means"
    print Chi2Mean
    print "Chi2Sigma"
    print Chi2Sigma

    assert(fracContamCuts[0] == -1)
    outfile = dataname
    if Rate_Model != 'discrete':
        print "outfile Pre Prefix"
        print outfile

        if cheatType:
            outfile = outfile + '_CheatType'
            if cheatZ:
                outfile = outfile + 'Z'
        elif cheatZ:
            outfile = outfile + '_CheatZ'

        outfile1 = outfile + '.txt'
        outfile2 = outfile + '-IndivBetaK.txt'
        output2 = open(outfile2, 'w')
        output2.write('i Beta_i k_i BetaErr_i kErr_i\n')
        for i, b, k, berr, kerr in zip(range(len(Betas)),Betas, ks, BetaErrs, kErrs):
            output2.write('{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f}\n'.format(i, b, k, berr, kerr))
        output2.close()
        print "Outfile Name"
        if not(os.path.isfile(outfile1)):
            output = open(outfile1, 'w')
            output.write('#Date Date/time at which job finished\n')
            output.write('#DataBeta Input beta for the simulated data sample. Will be 0.0 for real data.\n')
            output.write('#N_sims Number of datalike sims that go into the subsequent means\n')
            output.write('#SampleSize Mean Number of Events in data post cut\n')
            output.write('#delta_Beta mean difference between large MC sim beta (2.11 for the time being) and the measured beta for the data (not the beta in column 2.\n')
            output.write('#sigma_Beta stdev of delta_Beta over N_sims sims\n')
            output.write('#BetaStdErr std. error in the mean of delta_Beta over N_sims sims\n')
            output.write('#Beta_err mean statistical error on beta\n')
            output.write('#K mean ratio between large MC sim K (1.7E-5 for the time being) and the measured K for the data \n')
            output.write('#sigma_K stdev of K over N_sims sims\n')
            output.write('#KStdErr std. error in the mean of K over N_sims sims\n')
            output.write('#KStaterr mean statistical error on K\n')
            output.write('#meanZ mean photoZ of the large MC sim\n')
            output.write('#sigmaZ std. deviation of the photoZs for the large Sim\n')
            output.write('#sigmaDZ std. deviation of (zSim - zPHOT)\n')
            output.write('#NCC/NTotScaled overall CC Contamination after adjusting CC Frac to data\n')
            output.write('#NCC/NTot overall CC Contamination in sim only\n')
            output.write('#CCScales relative sim vs. CC rate in z-bins \n')
            output.write('#TypeChoice Internal Diagnostic, check code comments\n')
            output.write('#NNProbCut Threshold for NN probability of Ia\n')
            output.write('#NBins Number of Analysis Bins\n')
            output.write('#MRSLow Threshold for Neg Mures Outliers\n')
            output.write('#MRSHigh Threshold for Pos Mures Outliers\n')
            output.write('#FitprobCut Lowest Fitprob in sim\n')
            output.write('#MRSCut NSigma Hubble residual cut\n')
            output.write('#Chi2 minimum value of Chi2 function\n')
            output.write('#Correlation cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])\n')
            output.write('#Date \t\tDataBeta N_sims SampleSize delta_Beta sigma_Beta BetaStdErr BetaStatErr K sigma_K KStdErr KStatErr meanZ sigmaZ sigmaDz NCC/NTotScaled NCC/NTot CCScales TypeChoice NNProbCut NBins MRSLow MRSHigh FitprobCut MRSCut Chi2 Correlation\n')
        else:
            output = open(outfile1, 'a')
        print 'outfile'
        print outfile



        cat = RateTest.simcat.Catalog
        t = time.strftime('%b-%d-%H:%M')
        N_Sims = np.sum(np.invert(np.isnan(ks)))
        SigBeta = float(BetaSigma[0])
        SigK = float(ksigma[0])
        kStdErr = float(ksigma[0])/np.sqrt(N_Sims)
        BetaStdErr = float(BetaSigma[0])/np.sqrt(N_Sims)
        meanZ =  np.nanmean(cat[RateTest.ztype])
        sigZ = np.nanstd(cat[RateTest.ztype])
        sigDZ = np.nanstd(cat[RateTest.ztype] - cat['SIM_ZCMB'])
        lowzCut = zminFit
        highzCut = zmaxFit
        contam2 = np.sum(cat[(cat[RateTest.ztype] > lowzCut) & (cat[RateTest.ztype] < highzCut)]['SIM_TYPE_INDEX'] !=1).astype(float)/ float(cat[(cat[RateTest.ztype] > lowzCut) & (cat[RateTest.ztype] < highzCut)].shape[0])
        contam = RateTest.fracCCDataTot
        ccscales = np.average(np.array(CCScaleStorage),weights = 1.0/ma.masked_invalid(CCScaleErrStorage)**2, axis = 0)
        cov = RateTest.covar
        correlation = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
        print "Outfile debug"
        print t
        print trueBeta
        print N_Sims
        print BetaMean[0]
        print BetaStdErr
        print BetaErrs[0]
        print meanZ
        print sigZ
        print sigDZ
        print contam
        print RateTest.typeString
        print RateTest.postCutSimCat['NN_PROB_IA'].min()
        print SigBeta
        print kmean[0]
        print kErrs[0]
        print kStdErr
        print SigK
        print np.nanmean(SampleSizes)
        print int(nbinsFit)
        print ScaleMuResCutLow
        print ScaleMuResCutHigh
        print RateTest.postCutSimCat['FITPROB'].min()
        print MURESCuts
        print np.mean(Chi2Mean)
        print contam2
        print ccscales
        print correlation
        ccscales = ','.join(str(ccscales).split())
        output.write('{0}\t\t{1:.2f}\t{2}\t{17:.3f}\t{3:.3f}\t{12:.3f}\t{4:.3f}\t{5:.3f}\t{13:.3f}\t{14:.3f}\t{15:.3f}\t{16:.3f}\t{6:.3f}\t{7:.3f}\t{8:.3f}\t{9:.3f}\t{24:.3f}\t{25}\t{10}\t{11:.3f}\t{18:d}\t{19:.3f}\t{20:.3f}\t{21:.3f}\t{22:.2f}\t{23:.3f}\t{26:.3f}\n'.format(t, trueBeta, N_Sims, BetaMean[0], BetaStdErr, BetaErrs[0],meanZ, sigZ, sigDZ, contam,  RateTest.typeString, RateTest.postCutSimCat['NN_PROB_IA'].min(), SigBeta, kmean[0], kErrs[0], kStdErr,  SigK, np.nanmean(SampleSizes), int(nbinsFit), ScaleMuResCutLow, ScaleMuResCutHigh, RateTest.postCutSimCat['FITPROB'].min(), MURESCuts, np.mean(Chi2Mean), contam2, ccscales, correlation) )
        print "BetaMean[0]"
        print BetaMean[0]
        print BetaMean
        print "KMean[0]"
        print kmean[0]
        print kmean
        print "Correlation"

        print correlation
        #print "BetaWeightMean[0]"
        #print BetaWeightMean[0]
        #print BetaWeightMean
        #print "KWeightMean[0]"
        #print KWeightMean[0]
        #print KWeightMean
    if not noCCMC:
        print "Individual Scales"
        print CCScaleStorage
        print "Individual ScaleErrs"
        print CCScaleErrStorage
        print "average ScaleErrs"
        print np.nanmean(CCScaleErrStorage)
        print "AAA CC Scale means (weighted, unweighted)2"
        print np.average(ma.masked_invalid(np.array(CCScaleStorage)), weights = 1.0/ma.masked_invalid(CCScaleErrStorage)**2)
        print np.nanmean(ma.masked_invalid(np.array(CCScaleStorage)))

        print "AAA CC Scale stds"
        print np.nanstd(np.array(CCScaleStorage))
    if simInd == 1:
        plt.clf()
        hist, bins = np.histogram(CCScaleStorage, bins = np.linspace(0.0, 5.0, 10))
        plt.step((bins[1:]+bins[:-1])/2.0, hist, where = 'mid', c = 'g')
        plt.savefig(dataname + 'ScaleDistro.png')
        plt.clf()


    print "nIter"
    print nIter
    if not (priorRate is None):
        kPriorPlots = np.linspace(0.8, 1.5, 300)
        kPriors = []
        for ktemp in kPriorPlots:
            kPriors.append(ratePrior(ktemp*trueMCK, BetaMean[0]*trueMCBeta, priorRate, priorZEff, priorRateErrUp = ratePriorErrUp, priorRateErrDown = ratePriorErrDown, priorRateErrAll = ratePriorErrAll))


        betaPriorPlots = np.linspace(-0.5, 0.5, 300)
        betaPriors = []
        for btemp in betaPriorPlots:
            betaPriors.append(ratePrior(kmean[0]*trueMCK, b*trueMCBeta, priorRate, priorZEff, priorRateErrUp = ratePriorErrUp, priorRateErrDown = ratePriorErrDown, priorRateErrAll = ratePriorErrAll))

        actualPrior =  ratePrior(kmean[0]*trueMCK, BetaMean[0]*trueMCBeta, priorRate, priorZEff, priorRateErrUp = ratePriorErrUp, priorRateErrDown = ratePriorErrDown, priorRateErrAll = ratePriorErrAll)


        kPriors = np.array(kPriors)
        betaPriors = np.array(betaPriors)

        plt.clf()
        plt.figure()
        
        plt.plot(kPriorPlots, np.log10(kPriors) )
        plt.hlines(np.log10(actualPrior), kPriorPlots[0], kPriorPlots[-1], label = 'Best Fit Prior = {0:.03f}'.format(actualPrior))
        plt.vlines(kmean[0], np.log10(kPriors).min(), np.log10(kPriors).max(), label = 'Best Fit K = {0:.03f}'.format(kmean[0]))
        plt.xlabel('k')
        plt.ylabel('ratePrior')
        plt.legend()
        plt.savefig(dataname + '_LogKPriorPlot.png')

        

        plt.clf()
        plt.figure()
        plt.plot(kPriorPlots, kPriors)
        plt.hlines(actualPrior, kPriorPlots[0], kPriorPlots[-1], label = 'Best Fit Prior = {0:.03f}'.format(actualPrior))
        plt.vlines(kmean[0], kPriors.min(), kPriors.max(), label = 'Best Fit K = {0:.03f}'.format(kmean[0]))
        plt.xlabel('k')
        plt.ylabel('ratePrior')
        plt.legend()
        plt.savefig(dataname + '_KPriorPlot.png')

        plt.clf()
        plt.figure()
        plt.plot(betaPriorPlots, betaPriors)
        plt.hlines(actualPrior, betaPriorPlots[0], betaPriorPlots[-1], label = 'Best Fit Prior = {0:.03f}'.format(actualPrior))
        plt.vlines(BetaMean[0], betaPriors.min(), betaPriors.max(), label = 'Best Fit Beta = {0:.03f}'.format(BetaMean[0]))
        plt.xlabel('beta')
        plt.ylabel('ratePrior')
        plt.legend()
        plt.savefig(dataname + '_BetaPriorPlot.png')

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