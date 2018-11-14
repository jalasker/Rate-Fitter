#!/software/python-2.7-2014q3-el6-x86_64/bin/python
import SNANA_Reader as simread
import REAL_Reader as dataread
import astropy.cosmology as cosmo
import traceback
import scipy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import Cosmology
import emcee as MC
import corner
import scipy.stats.mstats as mstats
from sys import argv
#import inspect



class Rate_Fitter:
    def __init__(self, realfilename, realName, simfilename, simName, simgenfilename, MCBeta, zmin=0.1, zmax=1.20 , simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95,  Rate_Model = 'powerlaw',  cheatType = False, cheatZ = False, cuts = None):#, M, JDMin = 56535.0, MJDMax = 56716.0):
        self.zmin = zmin
        self.zmax = zmax
        self.MCBeta = MCBeta
        self.Rate_Model = Rate_Model
        self.cheatType = cheatType
        self.cheatZ = cheatZ
        self.cuts = cuts
        #print "Sim File Name"
        #print simfilename

        self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)#, MJDMin = MJDMin, MJDMax = MJDMax )
        #print self.simcat.Catalog.shape
        self.simName = simName
        self.simgencat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)#, MJDMin = MJDMin, MJDMax = MJDMax )
        SIMGEN = np.genfromtxt(simgenfilename, dtype=None, names = True, skip_footer=3, skip_header=1, invalid_raise=False)
        SIMGEN = SIMGEN[SIMGEN['GENZ'] != 'GENZ']
        #SIMGEN = SIMGEN[(SIMGEN['PEAKMJD'].astype(float) > MJDMin) & (SIMGEN['PEAKMJD'].astype(float) < MJDMax)]
        self.simgencat.params = {'flat':True, 'H0': simH0, 'Om0':simOmegaM, 'Ob0': simOb0, 'sigma8': simSigma8, 'ns': simNs}
        self.simgencat.cosmo = Cosmology.setCosmology('simCosmo', self.simcat.params)
        self.simgencat.OrigCatalog = np.copy(SIMGEN)
        self.simgencat.Catalog = np.copy(SIMGEN)
        self.simgencat.Catalog = self.simgencat.Catalog[self.simgencat.Catalog['GENZ'] != 'GENZ']
        self.simgencat.simname = simName
        self.simgencat.NSN = int(len(self.simgencat.Catalog['GENZ']))
        ##print self.simgencat.NSN
        ##print self.simcat.Catalog.shape
        print "SIMGEN NUMBER"
        print self.simgencat.NSN
        print "SIMGENCAT FILE"
        print simfilename
        #plt.figure()
        #MJDHist, MJDBins = np.histogram(self.simgencat.Catalog['PEAKMJD'].astype(float), bins = 100) #np.linspace(56535.0,  56716.0, 100))
        #plt.bar((MJDBins[1:] + MJDBins[:-1])/2.0, MJDHist, width = (MJDBins[1:] - MJDBins[:-1]))
        #plt.savefig('MJD_Sim_Dist.png')
        #plt.close()
        

        #obsCat = self.simgencat.Catalog[self.simgencat.Catalog['NOBS'].astype(int) > 0.01]

        #plt.figure()
        #MJDHist, MJDBins = np.histogram(obsCat['PKMJD'].astype(float), bins = 100) #np.linspace(56535.0,  56716.0, 100))
        #plt.bar((MJDBins[1:] + MJDBins[:-1])/2.0, MJDHist, width = (MJDBins[1:] - MJDBins[:-1]))
        #plt.savefig('MJD_SimObs_Dist.png')
        #plt.close()
        #print 'realfilename'
        #print realfilename
        self.realName = realName
        self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95 )

        #if (not (NNSimFile is None)) & (not (NNDataFile is None)) & (not (NNProbCut is None)): 
        #    self.NNSimFile = np.genfromtxt(NNSimFile, names = True, dtype = None, skip_header = 4)
        #    self.NNDataFile = np.genfromtxt(NNDataFile, names = True, dtype = None, skip_header = 4)
        #    self.realcat.Catalog = self.realcat.Catalog[self.NNDataFile['NN_PROB_IA'] > NNProbCut]
        #    self.simcat.Catalog = self.simcat.Catalog[self.NNSimFile['NN_PROB_IA'] > NNProbCut]
        
        #if not NNProbCut is None:
        #    self.realcat.Catalog = self.realcat.Catalog[self.realcat.Catalog['NN_PROB_IA'] > NNProbCut]
        #    self.simcat.Catalog = self.simcat.Catalog[self.simcat.Catalog['NN_PROB_IA'] > NNProbCut]



        if self.cheatType:
            print "WARNING, THE FITTER IS CHEATING AND ELIMINATED NON-IAs USING SIM INFO"
            self.realcat.Catalog = self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'] == 1]
            self.simcat.Catalog = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1]
        print 'N precuts'
        print self.realcat.Catalog['FITPROB'].shape
        for cut in cuts:
            self.realcat.Catalog = self.realcat.Catalog[(self.realcat.Catalog[cut[0]] > cut[1]) & (self.realcat.Catalog[cut[0]] < cut[2])]
            self.simcat.Catalog = self.simcat.Catalog[(self.simcat.Catalog[cut[0]] > cut[1]) & (self.simcat.Catalog[cut[0]] < cut[2])]
        print "Minimum Fitprob"
        print np.min(self.realcat.Catalog['FITPROB'])
        print 'N postcuts'
        print self.realcat.Catalog['FITPROB'].shape
        #if not(MURES_Cut is None):
        #    print "NO MORE STRAIGHT MURES CUTS"
        #    assert(0)
        '''
        try:

            #LCDM = cosmo.FlatLambdaCDM(70.0, 0.3)
            #MuRealModel = LCDM.distmod(self.realcat.Catalog['zPHOT'])
            #MuSimModel = LCDM.distmod(self.simcat.Catalog['zPHOT'])
            #alf = 0.146; bet = 3.03; M = -19.3
            #MuRealObs = self.realcat.Catalog['mB'] + 19.3 + alf*self.realcat.Catalog['x1'] - bet*self.realcat.Catalog['c']
            #MuSimObs = self.simcat.Catalog['mB'] + 19.3 + alf*self.simcat.Catalog['x1'] - bet*self.simcat.Catalog['c']
            #self.realcat.Catalog = self.realcat.Catalog[np.abs(MuRealModel - MuRealObs) < MURES_Cut]
            #self.simcat.Catalog = self.simcat.Catalog[np.abs(MuSimModel - MuSimObs) < MURES_Cut]

            self.realcat.Catalog = self.realcat.Catalog[np.abs(self.realcat.Catalog['MURES']) < MURES_Cut]
            self.simcat.Catalog = self.simcat.Catalog[np.abs(self.simcat.Catalog['MURES']) < MURES_Cut]
        except Exception, e: 
            print self.realcat.Catalog.dtype
            print "NO HUBBLE RESIDUALS. USE SALT2mu FIRST"
            print e
            traceback.print_exc()
            assert(0)
        '''
        #print "REALCAT SHAPE"
        #print self.realcat.Catalog.shape
        
    def newData(self, realfilename, realName):
        self.realName = realName
        self.realcat = simread.SNANA_Cat(realfilename, realName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95 )

        if self.cheatType:
            print "WARNING, THE FITTER IS CHEATING AND ELIMINATED NON-IAs USING SIM INFO"
            self.realcat.Catalog = self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'] == 1]
        
        print 'N precuts'
        print self.realcat.Catalog['FITPROB'].shape
        for cut in cuts:
            self.realcat.Catalog = self.realcat.Catalog[(self.realcat.Catalog[cut[0]] > cut[1]) & (self.realcat.Catalog[cut[0]] < cut[2])]
            #self.simcat.Catalog = self.simcat.Catalog[(self.simcat.Catalog[cut[0]] > cut[1]) & (self.simcat.Catalog[cut[0]] < cut[2])]
        print "Minimum Fitprob"
        print np.min(self.realcat.Catalog['FITPROB'])
        print 'N postcuts'
        print self.realcat.Catalog['FITPROB'].shape
        #if not(self.MURES_Cut is None):
        #    print "NO MORE STRAIGHT MURES CUTS"
        #    assert(0)
        '''
        try:
            self.realcat.Catalog = self.realcat.Catalog[np.abs(self.realcat.Catalog['MURES']) < MURES_Cut]
        except Exception, e: 
            print self.realcat.Catalog.dtype
            print "NO HUBBLE RESIDUALS. USE SALT2mu FIRST"
            print e
            traceback.print_exc()
            assert(0)
        '''
    def effCalc(self, fracContamCut = 0.0, nbins = 10):
        #### Do we want SNIas or all SN for efficiency?
        self.nbins = nbins
        #zPHOTs = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] > 0.5]['zHD'].astype(float)
        if self.cheatZ:
            ztype = 'SIM_ZCMB'
        else:
            ztype = 'zPHOT'
        #zTRUEs = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] > 0.5]['SIM_ZCMB'].astype(float)
        if (fracContamCut > 0.000000001) & (fracContamCut < 1.0):
            print "BPQWDSF Cutting based on Frac Contam"
            histTot, binsX, binsY = np.histogram2d(self.simcat.Catalog[ztype], self.simcat.Catalog['MURES'], bins = nbins)
            #histIa, binsX, binsY = np.histogram2D(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1]['zPHOT'], self.simcat.Catalog['MURES'], bins = self.bins)
            histCC, binsX, binsY = np.histogram2d(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1][ztype], self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1]['MURES'], bins = (binsX, binsY))

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
        
        zPHOTs = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1][ztype].astype(float)

        zTRUEs = self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1]['SIM_ZCMB'].astype(float)
        '''
        zPHOTs = self.simcat.Catalog[ztype].astype(float)

        zTRUEs = self.simcat.Catalog['SIM_ZCMB'].astype(float)
        '''
        print zPHOTs.shape
        print zTRUEs.shape
        
        binList = np.linspace(self.zmin, self.zmax, nbins+1)
        print binList
        self.binList = binList
        counts, zPhotEdges, zTrueEdges, binnumber = scipy.stats.binned_statistic_2d(zPHOTs, zTRUEs, zTRUEs, statistic = 'count', bins =  self.binList)
        assert(zPhotEdges.shape[0] == (self.nbins + 1))
        #zGenHist, zGenBins = np.histogram(self.simgencat.Catalog['GENZ'].astype(float), bins = self.binList)
        #zSim1Hist, zSim1Bins = np.histogram(self.simcat.Catalog['SIM_ZCMB'].astype(float), bins = self.binList)
        
        zGenHist, zGenBins = np.histogram(self.simgencat.Catalog[self.simgencat.Catalog['GENTYPE'] == 1]['GENZ'].astype(float), bins = self.binList)
        zSim1Hist, zSim1Bins = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] ==1]['SIM_ZCMB'].astype(float), bins = self.binList)
        

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


        extent = [zPhotEdges[0], zPhotEdges[-1], zTrueEdges[0], zTrueEdges[-1]]
        plt.figure()
        plt.imshow(np.flipud(counts), extent = extent, cmap = 'Blues')
        plt.colorbar()
        plt.savefig('redshiftDistro.png')
        plt.clf()
        plt.close()
        plt.figure()
        plt.imshow(np.flipud(self.effmat), extent = extent, cmap = 'Blues', norm=mpl.colors.LogNorm())
        plt.colorbar()
        plt.savefig('efficiencyMatrixLog.png')
        plt.clf()
        plt.close()
        plt.figure()
        plt.imshow(np.flipud(self.effmat), extent = extent, cmap = 'Blues')
        plt.colorbar()
        plt.savefig('efficiencyMatrix.png')
        plt.clf()
        plt.close()
        print 'effmat'
        print self.effmat
    def fit_rate(self):
        import iminuit as iM
        from iminuit import Minuit as M
        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        if self.cheatZ:
            ztype = 'SIM_ZCMB'
        else:
            ztype = 'zPHOT'
        plt.switch_backend('Agg')


        nSim, simBins = np.histogram(self.simgencat.Catalog['GENZ'].astype(float), bins=self.binList)
        nSim2, simBins2 = np.histogram(self.simcat.Catalog[ztype].astype(float), bins=self.binList)
        

        #nSim, simBins = np.histogram(self.simgencat.Catalog[self.simgencat.Catalog['GENTYPE'] == 1]['GENZ'].astype(float), bins=self.binList)
        #nSim2, simBins2 = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] ==1][ztype].astype(float), bins=self.binList)

        NCC , _ = np.histogram(self.simcat.Catalog[ztype][self.simcat.Catalog['SIM_TYPE_INDEX'] != 1].astype(float), bins=self.binList)

        FracBad = NCC*1.0/(1.0*nSim2)

        print "Fraction of CCs in each bin"
        print FracBad

        print 'NCC'
        print NCC

        print 'nSim2'
        print nSim2

        print "Fraction of CCs in each bin"
        print FracBad

        #### CHANGE WHEN USED WITH REAL DATA
        #bins, nSim = np.histogram(self.realcat., bins=np.linspace(self.zmin, self.zmax, self.nbins)) 
        #nData, dataBins = np.histogram(self.realcat.Catalog['zPHOT'].astype(float), bins=np.linspace(self.zmin, self.zmax, self.nbins + 1))
        '''
        if not(NNProbCut is None):
            NNProbs = self.realcat.Catalog['NN_PROB_IA']
            nData0, dataBins = np.histogram(self.realcat.Catalog['zHD'].astype(float), bins=self.binList)
            nData, databins, binnums = scipy.stats.binned_statistic(self.realcat.Catalog['zHD'].astype(float), NNProbs, bins = self.binList)
            nData = nData*(np.sum(nData0)*1.0/np.sum(nData))
        else:
            nData, dataBins = np.histogram(self.realcat.Catalog['zHD'].astype(float), bins=self.binList)
        '''
        nData, dataBins = np.histogram(self.realcat.Catalog[ztype].astype(float), bins=self.binList)

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
        

        nCCData = nData*FracBad
        print 'NCCData'
        print nCCData

        print "overall Contam"
        print np.sum(NCC)*1.0/(np.sum(nSim2)*1.0)
        '''
        nData = nData*FracGood

        print "nData, dataBins, realcat shape post contam correction"
        print nData
        print dataBins
        print np.sum(self.realcat.Catalog['zHD'].astype(float) > zmax)
        print np.sum(self.realcat.Catalog['zHD'].astype(float) < zmin)
        print self.realcat.Catalog['zHD'].shape
        '''
        if self.Rate_Model == 'powerlaw':
            def chi2func(nData, nSim, effmat, fnorm, zCenters, k, Beta, dump = False, complexdump = False, nCCData = None):
                #zBreak = 1.0
                Chi2Temp = 0.0
                #f_Js = []
                #temp = None
                #for zC in zCenters:
                #    if zC < zBreak:
                #        f_Js.append(k*(1+zC)**Beta)
                #    elif not(temp is None):
                #        f_Js.append(temp)
                #    else:
                #        temp = f_Js[-1]
                #        f_Js.append(temp)
                f_Js = k*(1+zCenters)**Beta
                chi2Mat = np.zeros((self.nbins))
                adjNMC = np.zeros((self.nbins))
                if nCCData is None:
                    nCCData = np.zeros(nData.shape)
                #print f_Js
                #Check if I am scaling errors down with increasing MC size. Make MC twice as large as "Data" to test.
                if dump: chi2Storage = []
                if dump: scaledNSimStor = []
                for row, nDataI, nCCDataI, i in zip(effmat, nData, nCCData, xrange(self.nbins)):
                    if dump:
                        print 'effmat row'
                        print row
                        print 'nDataI'
                        print nDataI
                        print 'nCCDataI'
                        print nCCDataI
                        scaledNSimTemp = 0.0
                    #if dump:
                    #    print "nDataI"
                    #    print nDataI
                    JSumTempNum = 0.0
                    JSumTempDen = 0.0
                    for eff, nSimJ, f_J, j in zip(row, nSim, f_Js, xrange(self.nbins)):
                        if dump and (i == j):
                            print 'NGen J'
                            print nSimJ
                            print 'JSumTempNum contr'
                            print nSimJ*f_J*eff*fnorm
                            print 'JSumTempDen contr'
                            print nSimJ*f_J*eff*fnorm*f_J*fnorm
                        if dump and (i != j) and self.cheatZ:
                            if nSimJ*f_J*eff*fnorm > 0:
                                print " This should be zero but isnt "
                                print nSimJ*f_J*eff*fnorm
                                assert(0)
                        JSumTempNum += nSimJ*f_J*eff*fnorm
                        JSumTempDen += nSimJ*f_J*eff*fnorm*f_J*fnorm

                    if dump:
                        print i
                        print 'fnorm'
                        print fnorm
                        print "JSumTempNum"
                        print JSumTempNum
                        print "JSumTempDen"
                        print JSumTempDen
                        print "Chi2Bin"
                        dataFunc = np.maximum(nDataI ,1)
                        #print dataFunc
                        c2t = (nDataI - nCCDataI - JSumTempNum)**2/( dataFunc + nCCDataI + JSumTempDen)
                        print c2t
                    if nDataI > 1E-11 or JSumTempDen > 1E-11:
                        dataFunc = np.maximum(nDataI ,1)
                        Chi2Temp += ((nDataI - nCCDataI - JSumTempNum)**2/(dataFunc + nCCDataI + JSumTempDen))#*fnorm**2
                return Chi2Temp
        elif self.Rate_Model == 'brokenpowerlaw':
            def chi2func(nData, nSim, effmat, fnorm, zCenters, k, Beta, zBreak, dump = False, complexdump = False, nCCData = None):
                #zBreak = 1.0
                Chi2Temp = 0.0
                f_Js = []
                temp = None
                if nCCData is None:
                    nCCData = np.zeros(nData.shape)
                for zC in zCenters:
                    if zC < zBreak:
                        f_Js.append(k*(1+zC)**Beta)
                    elif not(temp is None):
                        f_Js.append(temp)
                    else:
                        temp = f_Js[-1]
                        f_Js.append(temp)
                #f_Js = k*(1+zCenters)**Beta
                chi2Mat = np.zeros((self.nbins))
                adjNMC = np.zeros((self.nbins))
                #print f_Js
                #Check if I am scaling errors down with increasing MC size. Make MC twice as large as "Data" to test.
                if dump: chi2Storage = []
                if dump: scaledNSimStor = []
                for row, nDataI, nCCDataI, i in zip(effmat, nData, nCCData, xrange(self.nbins)):
                    if dump:
                        print 'effmat row'
                        print row
                        scaledNSimTemp = 0.0
                    #if dump:
                    #    print "nDataI"
                    #    print nDataI
                    JSumTempNum = 0.0
                    JSumTempDen = 0.0
                    for eff, nSimJ, f_J, j in zip(row, nSim, f_Js, xrange(self.nbins)):
                        if dump:
                            print 'JSumTempNum contr'
                            print nSimJ*f_J*eff*fnorm
                            print 'JSumTempDen contr'
                            print nSimJ*f_J*eff*fnorm*f_J*fnorm
                        JSumTempNum += nSimJ*f_J*eff*fnorm
                        JSumTempDen += nSimJ*f_J*eff*fnorm*f_J*fnorm
                    if nDataI > 1E-11 or JSumTempDen > 1E-11:
                        if dump:
                            print 'JSumTempNum'
                            print JSumTempNum
                            print 'JSumTempDen'
                            print JSumTempDen
                        Chi2Temp += ((nDataI - nCCDataI - JSumTempNum)**2/(nDataI + nCCDataI + JSumTempDen))#*fnorm**2
                return Chi2Temp
        zCenters = (simBins[1:] + simBins[:-1])/2.0
 
        fnorm = float(np.sum(nData))/float(np.sum(nSim))
        if self.Rate_Model == 'powerlaw':
            lamChi2 = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, nCCData = nCCData)
            lamChi2Dump = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, dump = True, nCCData = nCCData)
        elif self.Rate_Model == 'brokenpowerlaw':
            lamChi2 = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, 1.0, nCCData = nCCData)
            lamChi2Dump = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, 1.0, dump = True, nCCData = nCCData)

        MinObj = M(lamChi2, k = 1.0, error_k = 0.1, Beta = 0.0, error_Beta = 0.1, limit_k = (0.0, None))
        print "Chi2 init = {0}".format(round(lamChi2Dump(1.0, 0.0), 4))
        #MinObj = M(lamChi2, k = 1.0, fix_k = True, Beta = 0.0, error_Beta = 0.1)
        

        MinObj.set_strategy(2)
        fmin, param = MinObj.migrad(nsplit= 10)


        k = MinObj.values['k']
        kErr = MinObj.errors['k']
        Beta = MinObj.values['Beta']
        BetaErr = MinObj.errors['Beta']
        self.k = k
        self.Beta = Beta
        self.kErr = kErr
        self.BetaErr = BetaErr
        self.chi2 = fmin.fval#/(self.nbins - 2)

        print "Chi2 final = {0}".format(round(lamChi2Dump(self.k, self.Beta), 4))

        fJs = np.ones(zCenters.shape)
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


if __name__ == '__main__':
    from sys import argv
    datadir = argv[1]
    simdir = argv[2]
    dataname = argv[3]
    simname = argv[4]
    simgenfile = argv[5]
    NNCut = False
    cheatType = bool(int(argv[6]))
    cheatZ = bool(int(argv[7]))
    if( ('Combine' in simdir) or ('SALT2' in simdir)) &  (('Combine' in datadir) or ('SALT2' in simdir)):
        NNCut = True
        NNProbCut = 0.5
    
    #if len(argv) > 6:
    #    NNCut = True
    #    NNProbCut = 0.9
    #    NNData = argv[6]
    #    NNSim = argv[7]

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
        zmin = 0.1
        zmax = 1.2
        MJDMin = 0.0
        MJDMax = np.inf
        bins = "equalSize" 
        runFit = True
    kmean = []
    ksigma = []
    BetaMean = []
    BetaSigma= []
    Chi2Mean = []
    Chi2Sigma = []
    #MURES_Cuts = [2.0]
    #MURES_Cuts = [1.0, 1.5, 2.0, 3.0, 4.0, 99.0, 2.0]
    #for MURES_Cut in MURES_Cuts:
    fracContamCuts = [-1, 0.5]
    cuts = [('FITPROB', 0.01, np.inf)]
    for fcc in fracContamCuts:
        ks = []
        kErrs = []
        Betas = []
        BetaErrs = []
        Chi2s = []

        #ksPerf = []
        #kErrsPerf = []
        #BetasPerf = []
        #BetaErrsPerf = []
        nFail = 0
        simLoaded = False
        if '{' in datadir:
            nfile = 101
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
                        RateTest.newData(datadir.format(simInd), dataname.format(simInd))
                    else:
                        RateTest.newData(datadir.format(simInd), dataname.format(simInd))

                    RateTest.effCalc(nbins = 13, fracContamCut = fcc)
                    #RateTest.effCalc(nbins = 20)
                    RateTest.fit_rate()


                    ks.append(RateTest.k)
                    kErrs.append(RateTest.kErr)
                    Betas.append(RateTest.Beta)
                    BetaErrs.append(RateTest.BetaErr)
                    Chi2s.append(RateTest.chi2)
                except Exception, e:
                    print "FAILURE"
                    print e
                    traceback.print_exc()
                    nFail +=1
            else:
                try:
                    #if NNCut:
                    #    RateTest = Rate_Fitter(datadir.format(simInd), dataname.format(simInd), simdir, simname,simgenfile, 2.1, zmin = 0.1, zmax = 1.2, NNProbCut = NNProbCut, cheatZ = cheatZ, cheatType = cheatType)
                    #else:
                    RateTest = Rate_Fitter(datadir.format(simInd), dataname.format(simInd), simdir, simname,simgenfile, 2.1, zmin = 0.1, zmax = 1.2, cheatZ = cheatZ, cheatType = cheatType, cuts = cuts)# , MJDMin = 0, MJDMax = np.inf)

                    simLoaded = True
                    RateTest.effCalc(nbins = 13, fracContamCut = fcc)
                    #RateTest.effCalc(nbins = 20)
                    RateTest.fit_rate()


                    ks.append(RateTest.k)
                    kErrs.append(RateTest.kErr)
                    Betas.append(RateTest.Beta)
                    BetaErrs.append(RateTest.BetaErr)
                    Chi2s.append(RateTest.chi2)
                except Exception, e:
                    print "FAILURE"
                    print e
                    traceback.print_exc()
                    nFail +=1

        print "Number of Failures"
        print nFail
        print "mean k"
        print np.mean(ks)
        print "std. k"
        print np.std(ks)
        print "Mean beta"
        print np.mean(Betas)
        print "std. beta"
        print np.std(Betas)
        kmean.append(np.mean(ks))
        ksigma.append(np.std(ks))
        BetaMean.append(np.mean(Betas))
        BetaSigma.append(np.std(Betas))
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
            plt.savefig('Chi2Plot_CheatType.png')
        elif cheatZ and not cheatType:
            plt.savefig('Chi2Plot_CheatZ.png')
        elif cheatZ and cheatType:
            plt.savefig('Chi2Plot_CheatTypeZ.png')
        else:
            plt.savefig('Chi2Plot.png')



if cheatType:
    print "THESE RESULTS ONLY INCLUDE TRUE Ias BECAUSE WE CHEATED AND USED THE SIM INFORMATION"
if cheatZ:
    print "THESE RESULTS Use Simulated Redshift info"

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
'''