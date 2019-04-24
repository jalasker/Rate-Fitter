#!/software/python-2.7-2014q3-el6-x86_64/bin/python
import SNANA_Reader as simread
import REAL_Reader as dataread
import astropy.cosmology as cosmo
import traceback
import scipy
import scipy.stats as stats
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import Cosmology
import emcee as MC
import corner
import scipy.stats.mstats as mstats
from sys import argv
import glob
import time
import os
import gzip
import shutil
import numpy.ma as ma



class Rate_Fitter:
    def __init__(self, realfilename, realName, simfilename, simName, simgenfilename, MCBeta, zmin=0.1, zmax=1.20 , simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95,  Rate_Model = 'powerlaw',  cheatType = False, cheatZ = False, cuts = None, nprint = 5):
        self.zmin = zmin
        self.zmax = zmax
        self.MCBeta = MCBeta
        self.Rate_Model = Rate_Model
        self.cheatType = cheatType
        self.cheatZ = cheatZ
        self.cuts = cuts
        self.nprint = nprint
        try: 
            self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
        except:
            self.simcat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95, skip_header = 5)
        
        self.simName = simName
        self.simgencat = simread.SNANA_Cat(simfilename, simName, simOmegaM=0.3, simOmegaL=0.7, simH0=70.0, simw=-1.0, simOb0=0.049, simSigma8=0.81, simNs=0.95)
       
        try:
            SIMGEN = np.load(simgenfilename + '.npz')['a']
        except:
            SIMGEN = np.genfromtxt(simgenfilename, dtype=None, names = True, skip_footer=3, invalid_raise=False)
            np.savez_compressed(simgenfilename+'.npz', a = SIMGEN)
        
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
        
    def effCalc(self, fracContamCut = 0.0, nbins = 10, simIndex = 100):
        #### Do we want SNIas or all SN for efficiency?
        self.nbins = nbins
        self.typeString = ''

        if self.cheatZ:
            ztype = 'SIM_ZCMB'
        else:
            ztype = 'zPHOT'


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
            
    def fit_rate(self, fixK = False, fixBeta = False, simIndex = 100, trueBeta = 0):
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

        if simIndex < self.nprint:
            print "Type Location C"
            print "Choice C1"

        if len(self.typeString) <= 4:
            self.typeString = self.typeString + 'C1'
        nSim, simBins = np.histogram(self.simgencat.Catalog[self.simgencat.Catalog['GENTYPE'].astype(int) == 1]['GENZ'].astype(float), bins=self.binList)
        nSim2, simBins2 = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'].astype(int) ==1][ztype].astype(float), bins=self.binList)
        
        
       
        nSim3, simBins3 = np.histogram(self.simcat.Catalog[ztype].astype(float), bins=self.binList)
        

        NCC , _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] != 1][ztype].astype(float), bins=self.binList)

        
        NIa , _ = np.histogram(self.simcat.Catalog[self.simcat.Catalog['SIM_TYPE_INDEX'] == 1][ztype].astype(float), bins=self.binList)

        
        FracBad = NCC*1.0/(1.0*(NCC+NIa))
        nData, dataBins = np.histogram(self.realcat.Catalog[ztype].astype(float), bins=self.binList)
        nCCData = nData*FracBad

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


        try:
            TrueNCC, _ = np.histogram(self.realcat.Catalog[self.realcat.Catalog['SIM_TYPE_INDEX'] !=1][ztype].astype(float), bins=self.binList)
            if simIndex < self.nprint:

                print "True NCC Data"
                print TrueNCC
                print "Ratio NCCData/TrueNCC"
                print nCCData*1.0/(TrueNCC*1.0)
                print "Sum of TrueNCC"
                print np.sum(TrueNCC)
                print "Sum of Data"
                print np.sum(nData)
                print "Sum of nCCData"
                print np.sum(nCCData)
                print "True NCC/NData"
                print TrueNCC*1.0/(1.0*nData)
        except:
            print "Using real data"
        if simIndex < self.nprint:

            print "overall Contam"
            print np.sum(NCC)*1.0/(np.sum(nSim3)*1.0)
        
        def chi2func(nData, nSim, effmat, fnorm, zCenters, k, Beta, zBreak = 1.0, dump = False, complexdump = False, modelError = False, nIA = None, nCC = None, Rate_Model = 'powerlaw', zbins = None, simIndex = 100):
            Chi2Temp = 0.0
            if Rate_Model == 'powerlaw':
                f_Js = k*(1+zCenters)**Beta
            else:
                zCenters = (zbins[1:]+zbins[:-1])/2.0
                for zC in zCenters:
                    if zC < zBreak:
                        f_Js.append(k*(1+zC)**Beta)
                    elif not(temp is None):
                        f_Js.append(temp)
                    else:
                        temp = f_Js[-1]
                        f_Js.append(temp)
            chi2Mat = np.zeros((self.nbins))
            adjNMC = np.zeros((self.nbins))

            if (nIA is None) or (nCC is None):
                print "No CC Cut"
                fracCCData = np.zeros(nData.shape)
            else:
                print "Beta Adjusted CC Cut"
                #BetaRatio = k*(1+zCenters)**(Beta)#/(1+zCenters)**MCBeta
                BetaRatio = (1+zCenters)**(Beta)#/(1+zCenters)**MCBeta
                fracCCData = (nCC*1.0)/(1.0*(nCC + nIA*BetaRatio))
                if dump and (self.nprint < simInd):
                    print "fracCCData2"
                    print fracCCData
                    print "unscaled fracCCData"
                    print 1.0*(1.0*nCC)/(1.0*(nCC+nIA))
            
            nCCData = nData*fracCCData
            if dump and (self.nprint < simInd):
                print "nCCData2"
                print nCCData
                
            #print f_Js
            #Check if I am scaling errors down with increasing MC size. Make MC twice as large as "Data" to test.
            if dump: chi2Storage = []
            if dump: scaledNSimStor = []
            for row, nDataI, nCCDataI, i in zip(effmat, nData, nCCData, xrange(self.nbins)):
                #print "FUCKME"
                #print nDataI
                #print nCCDataI
                #assert((nDataI - nCCDataI) > 0)
                if dump and (self.nprint < simInd):
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

                if dump and (self.nprint < simInd):
                    print i
                    print 'nDataI'
                    print nDataI
                    print 'nCCDataI'
                    print nCCDataI
                    print 'fnorm'
                    print fnorm
                    print "JSumTempNum tot"
                    print JSumTempNum
                    print "JSumTempDen tot"
                    print JSumTempDen
                    print "Chi2Bin"
                    dataFunc = np.maximum(nDataI ,1)
                    #print dataFunc
                    c2t = (nDataI - nCCDataI - JSumTempNum)**2/( dataFunc + nCCDataI + JSumTempDen)
                    print c2t
                if nDataI < 5:
                    JSumTempDen += 10*nDataI
                #if nDataI > 1E-11 or JSumTempDen > 1E-11 and modelError:
                dataFunc = np.maximum(nDataI ,1)
                #    Chi2Temp += ((nDataI - nCCDataI - JSumTempNum)**2/(JSumTempNum + JSumTempDen))#*fnorm**2
                if nDataI > 1E-11 or JSumTempDen > 1E-11:
                    Chi2Temp += ((nDataI - nCCDataI - JSumTempNum)**2/(dataFunc + nCCDataI + JSumTempDen))#*fnorm**2

            return Chi2Temp
        
        zCenters = (simBins[1:] + simBins[:-1])/2.0
 
        #Is this right? Everything else in the other side of the chi2 function should be Ia only
        fnorm = float(np.sum(nData))/float(np.sum(nSim))
        if self.Rate_Model == 'powerlaw':
            lamChi2 = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, nIA = NIa, nCC = NCC, simIndex = simIndex)
            lamChi2Dump = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, dump = True, nIA = NIa, nCC = NCC, simIndex = simIndex)
        elif self.Rate_Model == 'brokenpowerlaw':
            lamChi2 = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, 1.0, nCCData = NCCData, simIndex = simIndex)
            lamChi2Dump = lambda k, Beta: chi2func(nData, nSim, self.effmat, fnorm, zCenters, k, Beta, 1.0, dump = True, nCCData = nCCData, simIndex = simIndex)

        MinObj = M(lamChi2, k = 1.0, error_k = 0.1, Beta = trueBeta, error_Beta = 0.1, limit_k = (0.0, None), fix_k = fixK, fix_Beta = fixBeta)
        print "Chi2 init = {0}".format(round(lamChi2Dump(1.0, 0.0), 4))
        #MinObj = M(lamChi2, k = 1.0, fix_k = True, Beta = 0.0, error_Beta = 0.1)
        

        MinObj.set_strategy(2)
        fmin, param = MinObj.migrad(nsplit= 10)

        

        lamChi2Dump(MinObj.values['k'], MinObj.values['Beta'])


        k = MinObj.values['k']
        kErr = MinObj.errors['k']
        Beta = MinObj.values['Beta']
        BetaErr = MinObj.errors['Beta']
        self.k = k
        self.Beta = Beta
        self.kErr = kErr
        self.BetaErr = BetaErr
        self.chi2 = fmin.fval#/(self.nbins - 2)

        self.BetaRatio = (1+zCenters)**(Beta)#/(1+zCenters)**MCBeta
        self.fracCCData = (NCC*1.0)/(1.0*(NCC + NIa*self.BetaRatio))

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

def getCCScaleSplit(simCat, dataCat, MURESWindow = (-1, 1), zbins = [0.0, 0.3, 0.6, 0.9, 1.2], muresBinslow = np.linspace(-2, -1, 4), muresBinshigh = np.linspace(1, 2, 4)):
    import iminuit as iM
    from iminuit import Minuit as M
    CCScaleslow = []
    CCScaleshigh = []
    fnorm2 = float(dataCat.shape[0])/simCat.shape[0]

    simCatlow = simCat[(simCat['MURES'] < MURESWindow[0]) ]#| (simCat['MURES'] > MURESWindow[1]) ]
    dataCatlow = dataCat[(dataCat['MURES'] < MURESWindow[0])]# | (dataCat['MURES'] > MURESWindow[1]) ]
    simCathigh = simCat[(simCat['MURES'] > MURESWindow[1]) ]#| (simCat['MURES'] > MURESWindow[1]) ]
    dataCathigh = dataCat[(dataCat['MURES'] > MURESWindow[1])]# | (dataCat['MURES'] > MURESWindow[1]) ]
    for zl, zh in zip(zbins[:-1], zbins[1:]):
        tempSimlow = simCatlow[(simCatlow['zPHOT'] < zh) & (simCatlow['zPHOT'] > zl)]
        tempDatalow = dataCatlow[(dataCatlow['zPHOT'] < zh) & (dataCatlow['zPHOT'] > zl)]
        tempSimhigh = simCathigh[(simCathigh['zPHOT'] < zh) & (simCathigh['zPHOT'] > zl)]
        tempDatahigh = dataCathigh[(dataCathigh['zPHOT'] < zh) & (dataCathigh['zPHOT'] > zl)]
        if type(muresBinslow == int):
            histDlow, muresBinslow = np.histogram(tempDatalow['MURES'], bins = muresBinslow)
            binsDlow = muresBinslow
        else: 
            histDlow, binsDlow = np.histogram(tempDatalow['MURES'], bins = muresBinslow)


        histSlow, binsSlow  = np.histogram(tempSimlow['MURES'], bins = muresBinslow)
        histSCClow, binsSCClow  = np.histogram(tempSimlow['MURES'][tempSimlow['SIM_TYPE_INDEX'] != 1], bins = muresBinslow)
        histSIalow, binsSIalow  = np.histogram(tempSimlow['MURES'][tempSimlow['SIM_TYPE_INDEX'] == 1], bins = muresBinslow)

        if type(muresBinshigh == int):
            histDhigh, muresBinshigh = np.histogram(tempDatahigh['MURES'], bins = muresBinshigh)
            binsDhigh = muresBinshigh
        else: 
            histDhigh, binsDhigh = np.histogram(tempDatahigh['MURES'], bins = muresBinshigh)

        histShigh, binsShigh  = np.histogram(tempSimhigh['MURES'], bins = muresBinshigh)
        histSCChigh, binsSCChigh  = np.histogram(tempSimhigh['MURES'][tempSimhigh['SIM_TYPE_INDEX'] != 1], bins = muresBinshigh)
        histSIahigh, binsSIahigh  = np.histogram(tempSimhigh['MURES'][tempSimhigh['SIM_TYPE_INDEX'] == 1], bins = muresBinshigh)

        #histDhigh, binsDhigh = np.histogram(tempDatahigh['MURES'], bins = muresBinshigh)
        print "low MURES tail"
        print histSlow
        print histSCClow
        print histSIalow
        print "high MURES tail"
        print histShigh
        print histSCChigh
        print histSIahigh

        print "Data tails"
        print histDlow
        print histDhigh


        def chi2(SimCC, SimIa, Data, scale, fnorm2):
            #ms = ma.masked_array((SimCC*scale + SimIa)*fnorm2, mask = (SimCC+SimIa) == 0)
            ms = (SimCC*scale + SimIa)*fnorm2
            return stats.chisquare(Data, ms)[0]

        lamChi2low = lambda s: chi2(histSCClow, histSIalow, histDlow, s, fnorm2)
        lamChi2high = lambda s: chi2(histSCChigh, histSIahigh, histDhigh, s, fnorm2)
        MinObjlow = M(lamChi2low, s = 1.0, error_s = 1, limit_s = (0.0, None))
        MinObjhigh = M(lamChi2high, s = 1.0, error_s = 1, limit_s = (0.0, None))
        print "Chi2 init low = {0}".format(round(lamChi2low(1.0), 4))
        print "Chi2 init high = {0}".format(round(lamChi2high(1.0), 4))
        #MinObj = M(lamChi2, k = 1.0, fix_k = True, Beta = 0.0, error_Beta = 0.1)
        

        MinObjlow.set_strategy(2)
        MinObjhigh.set_strategy(2)
        fminlow, paramlow = MinObjlow.migrad(nsplit= 10)
        fminhigh, paramhigh = MinObjhigh.migrad(nsplit= 10)

        CCScaleslow.append(MinObjlow.values['s'])
        CCScaleshigh.append(MinObjhigh.values['s'])
    return CCScaleslow, CCScaleshigh

def getCCScale(simCat, dataCat, MURESWindow = (-1, 1), zbins = [0.0, 0.3, 0.6, 0.9, 1.2], muresBins = np.linspace(-2, -1, 4), Beta = 0.0, binList = None, fracCCData = None, outfilePrefix = 'Test'):
    import iminuit as iM
    from iminuit import Minuit as M
    CCScales = []
    CCScaleErrs = []

    tempSimCC = simCat[simCat['SIM_TYPE_INDEX'].astype(int) != 1]
    tempSimIa = simCat[simCat['SIM_TYPE_INDEX'].astype(int) == 1]

    allSimCC = np.copy(tempSimCC)
    allSimIa = np.copy(tempSimIa)
    allData = np.copy(dataCat)
    '''
    print "total sim CC, Ia, all including outliers"
    print tempSimCC.shape
    print tempSimIa.shape
    print tempSimIa.shape[0] + tempSimCC.shape[0]
    NCCTOTSIM = tempSimCC.shape[0]
    NIATOTSIM = tempSimIa.shape[0]
    NTOTSIM = tempSimIa.shape[0] + tempSimCC.shape[0]

    print "total data, all including outliers"
    print dataCat.shape
    NTOTDATA = dataCat.shape[0]
    '''
    simHistCC, simBinsCC  = np.histogram(tempSimCC['MURES'], bins = muresBins)
    simHistIa, simBinsIa  = np.histogram(tempSimIa['MURES'], bins = muresBins)

    simZHistCC, simZBinsCC = np.histogram(tempSimCC['zPHOT'], bins = binList)
    simZHistIa, simZBinsIa = np.histogram(tempSimIa['zPHOT'], bins = binList)

    binCent = (simZBinsIa[1:] + simZBinsIa[:-1])/2.0

    #dataIaHist = dataHist*(1-fracCCData)
    #dataCCHist = dataHist*fracCCData

    print "components of simHist"
    print simZHistIa
    print binCent
    print Beta
    print simHistCC


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
    

    #simCat = simCat[ (simCat['MURES'] > MURESWindow[1]) ]
    #dataCat = dataCat[(dataCat['MURES'] > MURESWindow[1]) ]

    print "postcut"
    print simCat.shape
    print dataCat.shape

    '''
    print "outlier sim CC, Ia, total"
    print np.sum([simCat['SIM_TYPE_INDEX'].astype(int) != 1])
    print np.sum([simCat['SIM_TYPE_INDEX'].astype(int) == 1])
    print simCat.shape
    noutccsim = np.sum([simCat['SIM_TYPE_INDEX'].astype(int) != 1])
    noutiasim = np.sum([simCat['SIM_TYPE_INDEX'].astype(int) == 1])
    noutsim = simCat.shape[0]
    print "outlier data total"
    print dataCat.shape
    noutdata = dataCat.shape[0]
    '''    
    print 'post MURES Cut Ndata and Nsim'
    print dataCat.shape
    print simCat.shape
    '''
    print "---Rick Test Here---"
    R = float(noutdata)/float(NTOTDATA)
    print "R"
    print R

    S = float(R*NIATOTSIM - noutiasim)/float(noutccsim - R*NCCTOTSIM)

    print "S"
    print S

    print "---End Rick Test---"
    '''
    for zl, zh in zip(zbins[:-1], zbins[1:]):

        tempSim = simCat[(simCat['zPHOT'] < zh) & (simCat['zPHOT'] > zl)]
        tempData = dataCat[(dataCat['zPHOT'] < zh) & (dataCat['zPHOT'] > zl)]


        allSimCCZbin = allSimCC[(allSimCC['zPHOT'] < zh) & (allSimCC['zPHOT'] > zl)]
        allSimIaZbin = allSimIa[(allSimIa['zPHOT'] < zh) & (allSimIa['zPHOT'] > zl)]



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
        histSIa, binsSIa  = np.histogram(tempSimIa['MURES'], bins = muresBins)
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

        
        histSIa = histSIa*((1+np.mean(tempSim['zPHOT']))**Beta)

        print "beta scaling"
        print ((1+np.mean(tempSim['zPHOT']))**Beta)

        print "mean zphot"
        print np.mean(tempSim['zPHOT'])

        print "beta"
        print Beta

        print "PostCorr HistSIa "

        print histSIa

        #R = np.array(histD).astype(float)/np.array(histDAll).astype(float)
        R = np.array(histD).astype(float)/float(allDataZbin.shape[0])

        print "R"

        print R

        print "Hist CC, outlier and total"

        print histSCC
        #print histSAllCC
        print allSimCCZbin.shape[0]

        #S = float(np.array(R*histSAllIa) - np.array(histSIa))/float(np.array(histSCC) - np.array(R*histSAllCC))
        S = float(np.array(R*allSimIaZbin.shape[0]) - np.array(histSIa))/float(np.array(histSCC) - np.array(R*allSimCCZbin.shape[0]))

        print "S"
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



        CCScales.append(S)
        CCScaleErrs.append(dS[0])



   

        '''
        def chi2(SimCC, SimIa, Data, scale, fnorm2):
            #ms = ma.masked_array((SimCC*scale + SimIa)*fnorm2, mask = (SimCC+SimIa) == 0)
            #totSimCC = np.sum(simcat['SIM_TYPE_INDEX'] !=1)
            #totSimIa = np.sum(simcat['SIM_TYPE_INDEX'] ==1)
            #fnorm2 = float(datacat.shape[0])/(totSimCC*scale+totSimIa)
            ms = (SimCC*scale + SimIa)*fnorm2

            chi2 = np.sum((Data - ms)**2/(np.maximum(Data*1.0, 1.0)))
            return chi2
            #return stats.chisquare(Data, ms)[0]

        #lamChi2 = lambda s: chi2(histSCC, histSIa, histD, s, dataCat, simCat)
        lamChi2 = lambda s: chi2(histSCC, histSIa, histD, s, fnorm2)
        
        MinObj = M(lamChi2, s = 1.0, error_s = 1, limit_s = (0.0, None))
        
        print "Chi2 init low = {0}".format(round(lamChi2(1.0), 4))
        #MinObj = M(lamChi2, k = 1.0, fix_k = True, Beta = 0.0, error_Beta = 0.1)
        

        MinObj.set_strategy(2)
        fmin, param = MinObj.migrad(nsplit= 10)

        CCScales.append(MinObj.values['s'])
        CCScaleErrs.append(MinObj.errors['s'])
        '''
        #plt.plot((simBinsCC[1:] + simBinsCC[:-1])/2.0, simHistCC, c = 'b', label = 'prescaled Sim CC')
        #plt.plot((simBinsCC[1:] + simBinsCC[:-1])/2.0, CCScales[0]*simHistCC*fnorm2, c = 'g', label = 'postscaledSimCC')
        #plt.plot(muresBins, histD, c = 'r', label = 'data')
        plt.step((simBinsCC[1:] + simBinsCC[:-1])/2.0, simHistCC, c = 'b', where = 'mid',  label = 'prescaled Sim CC')
        plt.step((simBinsCC[1:] + simBinsCC[:-1])/2.0, CCScales[0]*simHistCC*fnorm2, c = 'g', where = 'post',  label = 'postscaledSimCC')
        plt.step((muresBins[1:] + muresBins[:-1])/2.0, histD, c = 'r', where = 'mid',  label = 'data')
        plt.savefig(outfilePrefix + 'ScaledHist.png')
        plt.clf()
    return CCScales, CCScaleErrs


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
    combinedScaleFit = True
    ScaleMuResCutLow = -1
    ScaleMuResCutHigh = 2.5
    muresBins = 1
    muresBinsLow  = 3
    muresBinsHigh = 3
    scaleZBins = [0.0, 1.2]

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
    if combinedScaleFit:

        CCScaleStorageGlobal = []
        CCScaleErrStorageGlobal = []

    else:
        CCScaleStorageLowGlobal = []
        CCScaleStorageHighGlobal = []
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
        if combinedScaleFit:

            CCScaleStorage = []
            CCScaleErrStorage = []

        else:
            CCScaleStorageLow = []
            CCScaleStorageHigh = []
        #ksPerf = []
        #kErrsPerf = []
        #BetasPerf = []
        #BetaErrsPerf = []
        nFail = 0
        simLoaded = False
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
                    else:
                        RateTest.newData(datadir.format(simInd), dataname.format(simInd), simIndex = simInd)


                   

                    RateTest.effCalc(nbins = nbins, fracContamCut = fcc, simIndex = simInd)
                    #RateTest.effCalc(nbins = 20)
                    RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, simIndex = simInd, trueBeta = trueBeta - 2.11)

                    print "AAA CC Scales"
                    if combinedScaleFit:
                        CCScale, CCScaleErr =  getCCScale(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBins = muresBins, Beta = RateTest.Beta, binList = RateTest.binList, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname)
                        print CCScale
                        CCScaleStorage.append(CCScale[0])
                        CCScaleErrStorage.append(CCScaleErr[0])
                    else:

                        CCScaleLow, CCScaleHigh =  getCCScaleSplit(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBinslow = muresBinsLow, muresBinshigh =muresBinsHigh, Beta = RateTest.Beta, binList = RateTest.binList, outfilePrefix =  dataname)
                        print CCScaleLow
                        print CCScaleHigh
                        CCScaleStorageLow.append(CCScaleLow)
                        CCScaleStorageHigh.append(CCScaleHigh)


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
                    #if NNCut:
                    #    RateTest = Rate_Fitter(datadir.format(simInd), dataname.format(simInd), simdir, simname,simgenfile, 2.1, zmin = 0.1, zmax = 1.2, NNProbCut = NNProbCut, cheatZ = cheatZ, cheatType = cheatType)
                    #else:
                    RateTest = Rate_Fitter(datadir.format(simInd), dataname.format(simInd), simdir, simname,simgenfile, 2.1, zmin = 0.1, zmax = 1.2, cheatZ = cheatZ, cheatType = cheatType, cuts = cuts)# , MJDMin = 0, MJDMax = np.inf)
                    

                    simLoaded = True
                    RateTest.effCalc(nbins = nbins, fracContamCut = fcc)
                    #RateTest.effCalc(nbins = 20)
                    RateTest.fit_rate(fixK = fixK, fixBeta = fixBeta, trueBeta = trueBeta - 2.11)

                    print "AAA CC Scales"
                    if combinedScaleFit:
                        CCScale, CCScaleErr =  getCCScale(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBins = muresBins, Beta = RateTest.Beta, binList = RateTest.binList, fracCCData = RateTest.fracCCData, outfilePrefix =  dataname)
                        print 'CC Scale'
                        print CCScale
                        CCScaleStorage.append(CCScale[0])
                        CCScaleErrStorage.append(CCScaleErr[0])
                    else:

                        CCScaleLow, CCScaleHigh =  getCCScaleSplit(RateTest.simcat.Catalog, RateTest.realcat.Catalog, MURESWindow = (ScaleMuResCutLow, ScaleMuResCutHigh), zbins = scaleZBins, muresBinslow = muresBinsLow, muresBinshigh =muresBinsHigh, Beta = RateTest.Beta, binList = RateTest.binList, outfilePrefix =  dataname)
                        print CCScaleLow
                        print CCScaleHigh
                        CCScaleStorageLow.append(CCScaleLow)
                        CCScaleStorageHigh.append(CCScaleHigh)


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
        print "std. k"
        print np.std(ks)
        print "Mean beta"
        print np.mean(Betas)
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

        if combinedScaleFit:
            print "AAA CC Scale means (weighted, unweighted)"
            print np.average(ma.masked_invalid(np.array(CCScaleStorage)),weights = 1.0/ma.masked_invalid(CCScaleErrStorage)**2, axis = 0)
            print np.mean(ma.masked_invalid(np.array(CCScaleStorage)), axis = 0)
            print "AAA CC Scale stds"
            print np.nanstd(np.array(CCScaleStorage), axis = 0)
            CCScaleStorageGlobal.append(CCScaleStorage)

        else:
            print "AAA CC Scale means Low"
            print np.nanmean(np.array(CCScaleStorageLow), axis = 0)
            print "AAA CC Scale stds Low"
            print np.nanstd(np.array(CCScaleStorageLow), axis = 0)
            CCScaleStorageLowGlobal.append(CCScaleStorageLow)
            print "AAA CC Scale means High"
            print np.nanmean(np.array(CCScaleStorageHigh), axis = 0)
            print "AAA CC Scale stds High"
            print np.nanstd(np.array(CCScaleStorageHigh), axis = 0)
            CCScaleStorageHighGlobal.append(CCScaleStorageHigh)



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

if combinedScaleFit:
    print "Individual Scales"
    print CCScaleStorage
    print "Individual ScaleErrs"
    print CCScaleErrStorage
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

else:
    print "AAA CC Scale means Low"
    print np.nanmean(np.array(CCScaleStorageLow))
    print "AAA CC Scale stds Low"
    print np.nanstd(np.array(CCScaleStorageLow))

    print "AAA CC Scale means High"
    print np.nanmean(np.array(CCScaleStorageHigh))
    print "AAA CC Scale stds High"
    print np.nanstd(np.array(CCScaleStorageHigh))

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