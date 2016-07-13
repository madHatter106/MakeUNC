# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 15:09:44 2015
script to run analytic uncertainty analysis of rrs data
@author: ekarakoy
"""
import numpy as np
import netCDF4 as nc
import glob
import sys
import re
import argparse
import logging
import datetime as dt

class MakeUnc(object):
    """
    Class to get Rrs uncertainties for a given L2 granule. Includes methods to:
        * get baseline from L2 granule
        * calculate uncertainties (as rmse) from corresponding perturbed L2 files
        * save uncertainty variables in original unperturbed granule
    Inputs:
        *args:
            1- baselineFile
            2- noisyDir -- directory where noisy files are located
        **kwargs:
            1- fnum -- number of files to process [=None |1,2,...]
            2- doSaniCheck -- recover noise model by calculating Lt_unc [=False | True]

    """
    def __init__(self,silFile,noisyDir,**kwargs):
        self.verbose = kwargs.pop('verbose',False)
        self._SetLogger()
        self.silFile = silFile
        self.noisyDir = noisyDir
        self.fnum = kwargs.pop("fnum",None)
        self.doSaniCheck = kwargs.pop("doSaniCheck",False)
        self.doChla = kwargs.pop("doChla",False)
        self.doNflh = kwargs.pop("doNflh",False)
        self.rrsSilDict = dict.fromkeys(self.bands)
        self.attrRrsUncDict = dict.fromkeys(self.bands)
        self.dimsDict = dict.fromkeys(self.bands)
        self.dTypeDict = dict.fromkeys(self.bands)
        self.rrsUncArrDict = dict.fromkeys(self.bands)
        if self.doSaniCheck:
            self.ltUncArrDict = dict.fromkeys(self.bands)
            self.ltSilDict = dict.fromkeys(self.bands)
            self.attrLtUncDict = dict.fromkeys(self.bands)
        otherProdkeys = []
        if self.doChla:
            otherProdkeys.append('chlor_a')
            otherProdkeys.append('chlor_a_unc')
        if self.doNflh:
            otherProdkeys.append('nflh')
            otherProdkeys.append('nflh_unc')
        if len(otherProdkeys) > 0:
            self.otherProdsDict = dict.fromkeys(otherProdkeys)
            attrUncKeys = [x for x in otherProdkeys if re.search('_unc',x)]
            self.attrOtherProdUncDict = dict.fromkeys(attrUncKeys)
            self.dTypeDict.update(dict.fromkeys(attrUncKeys))
            self.dimsDict.update(dict.fromkeys(attrUncKeys))

    def _SetLogger(self):
        self.logger = logging.getLogger(__name__)
        logfilename = './MakeUnc' + str(dt.datetime.now())
        fh = logging.FileHandler(logfilename,mode='w')
        if self.verbose:
            formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(name)s, %(funcName)s, %(lineno)d: %(message)s')
            fh.setLevel(logging.DEBUG)
        else:
            formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s')
            fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        self.logger.info('Logger initialized')
        return None

    def WriteToSilent(self):
        # first create NC variables if necessary
        # save unc in corresponding NC variables.
        with nc.Dataset(self.silFile,'r+') as dsSil:
            geoGr = dsSil.groups['geophysical_data']
            geoVar = geoGr.variables
            for band in self.bands:
                rrsUncProdName  = 'Rrs_unc_' + band
                if rrsUncProdName not in geoVar:
                    varRrsUnc = geoGr.createVariable(rrsUncProdName,
                                                    self.dTypeDict[band],
                                                    self.dimsDict[band])
                    varRrsUnc.setncatts(self.attrRrsUncDict[band])
                else:
                    varRrsUnc = geoVar[rrsUncProdName]
                varRrsUnc[:] = self.rrsUncArrDict[band]
                if self.doSaniCheck:
                    ltUncProdName = 'Lt_unc_' + band
                    if ltUncProdName not in geoVar:
                        varLtUnc = geoGr.createVariable(ltUncProdName,
                                                        self.dTypeDict[band],
                                                        self.dimsDict[band])
                        varLtUnc.setncatts(self.attrLtUncDict[band])
                    else:
                        varLtUnc = geoVar[ltUncProdName]
                    varLtUnc[:] = self.ltUncArrDict[band]
            if self.doChla:
                if 'chlor_a_unc' not in geoVar:
                    varChloraUnc = geoGr.createVariable('chlor_a_unc',
                                                        self.dTypeDict['chlor_a_unc'],
                                                        self.dimsDict['chlor_a_unc'])
                    varChloraUnc.setncatts(self.attrOtherProdUncDict['chlor_a_unc'])
                else:
                    varChloraUnc = geoVar['chlor_a_unc']
                varChloraUnc[:] = self.otherProdsDict['chlor_a_unc']

            if self.doNflh:
                if 'nflh_unc' not in geoVar:
                    print('nflh_unc there; creating variable...',flush=True)
                    varNflhUnc = geoGr.createVariable('nflh_unc',
                                                      self.otherProdsDict['nflh_unc'].dtype,
                                                        self.dimsDict['nflh_unc'])
                    varNflhUnc.setncatts(self.attrOtherProdUncDict['nflh_unc'])
                else:
                    print('nflh_unc available, using existing variable...',
                          flush=True)
                    varNflhUnc = geoVar['nflh_unc']
                varNflhUnc[:] = self.otherProdsDict['nflh_unc']

    def BuildUncs(self,noisySfx,verbose=False):
        """"
        Calculates rrs uncertainty as st.dev of rrs. Note that to save time
            I use unperturbed rrs as the rrs average for the simulation
        """
        fBaseName = self.silFile.split('/')[-1].split('.')[0].split('_')[0]
        matchFilPatt = self.noisyDir +  fBaseName + '*' + noisySfx + '*'
        if verbose:
            print ("searching for " + matchFilPatt)
        firstPass = [True] * len(self.bands)
        flis = glob.glob(matchFilPatt)
        lflis = len(flis)
        rrsAggDataDict = dict.fromkeys(self.bands)
        if self.doSaniCheck:
            ltAggDataDict = dict.fromkeys(self.bands)
        if self.doChla:
            chlAggDataArr = np.array([])
        if self.doNflh:
            nflhAggDataArr = np.array([])
        #process noisy data
        for fcount,fname in enumerate(flis):
            prcDone = 100 * fcount / (lflis - 1)
            if verbose:
                print("\rLoading and reading %s -- %.1f%%" %
                                                    (fname,prcDone),end='',flush=True)
            else:
                print("\r%.1f%%" % prcDone,end='',flush=True)
            with nc.Dataset(fname) as nds:
                nGeoGr = nds.groups['geophysical_data']
                nGeoVar = nGeoGr.variables
                for i,band in enumerate(self.bands):
                    noisyRrs = nGeoVar['Rrs_'+band][:]
                    if self.doSaniCheck:
                        noisyLt = nGeoVar['Lt_'+band][:]
                    if self.doChla:
                        noisyChl = nGeoVar['chlor_a'][:]
                    if self.doNflh:
                        noisyNflh = nGeoVar['nflh'][:]
                    if firstPass[i]:
                        rrsAggDataDict[band] = (noisyRrs - self.rrsSilDict[band]) ** 2
                        if self.doSaniCheck:
                            ltAggDataDict[band] = (noisyLt - self.ltSilDict[band]) ** 2
                        if self.doChla:
                            chlAggDataArr = (noisyChl - self.otherProdsDict['chlor_a']) ** 2
                        if self.doNflh:
                            nflhAggDataArr = (noisyNflh - self.otherProdsDict['nflh']) ** 2
                        firstPass[i] = False
                        if verbose:
                            print("\nFirst pass complete for band %s" % band,end='...',flush=True)
                    else:
                        rrsAggDataDict[band] += (noisyRrs - self.rrsSilDict[band]) ** 2
                        if self.doSaniCheck:
                            ltAggDataDict[band] += (noisyLt - self.ltSilDict[band]) ** 2
                        if self.doChla:
                            chlAggDataArr += (noisyChl - self.otherProdsDict['chlor_a']) ** 2
                        if self.doNflh:
                            nflhAggDataArr += (noisyNflh - self.otherProdsDict['nflh']) ** 2
                if self.fnum is not None:
                    # if a number of files to process was specified & reached, break.
                    if (fcount + 1) > self.fnum:
                        break

        for band in self.bands:
            if verbose:
                print("\n...computing stdev for band",band,flush=True)
            self.rrsUncArrDict[band] = np.ma.sqrt(rrsAggDataDict[band] / lflis)
            if self.doSaniCheck:
                self.ltUncArrDict[band] = np.sqrt(ltAggDataDict[band] / lflis)
            if self.doChla:
                self.otherProdsDict['chlor_a_unc'] = np.ma.sqrt(chlAggDataArr / lflis)
            if self.doNflh:
                self.otherProdsDict['nflh_unc'] = np.ma.sqrt(nflhAggDataArr / lflis)
        if verbose:
            print("\nProcessed %d files " % lflis)

    def ReadFromSilent(self):
        '''Reads Baseline file
            Flags: l2bin default flags, namely ATMFAIL(1), LAND(2), HIGLINT(8),
                HILT(16), HISATZEN(32), STRAYLIGHT(256), CLDICE(512), COCCOLITH(1024),
                HISOLZEN(4096), LOWLW(16384), CHLFAIL(32768), NAVWARN(65536),
                MAXAERITER(524288), CHLWARN(2097152), ATMWARN(4194304),
                NAVFAIL(33554432), FILTER(67108864)
        '''
        #flagBits = 1 + 2 + 8 + 16 + 32 +  256 + 512 + 1024 + 4096 + 16384 + \
         #       32768 +  65536 + 524288 + 2097152 + 4194304 + 33554432 + 67108864

        #l2flags = geoVar['l2_flags'][:]
        #flagMaskArr = (l2flags & flagBits > 0)
        with nc.Dataset(self.silFile,'r') as dsSil:
            geoGr = dsSil.groups['geophysical_data']
            geoVar = geoGr.variables
            for band in self.bands:
                rrs = geoVar['Rrs_'+band]
                self.rrsSilDict[band] = rrs[:]
                self.attrRrsUncDict[band] = {'long_name' : 'Uncertainty in ' +
                                                          rrs.long_name,
                                             '_FillValue': rrs._FillValue,
                                             'units': rrs.units,
                                             'scale_factor':rrs.scale_factor,
                                             'add_offset':rrs.add_offset,}
                self.dimsDict[band] = rrs.dimensions
                self.dTypeDict[band] = rrs.dtype
                if self.doSaniCheck:
                    lt = geoVar['Lt_'+band]
                    self.ltSilDict[band] = lt[:]
                    self.attrLtUncDict[band] = {'long_name' : 'Uncertainty in ' +
                                                          lt.long_name,
                                             '_FillValue': lt._FillValue,
                                             'units': lt.units}
            if self.doChla:
                chla = geoVar['chlor_a']
                self.otherProdsDict['chlor_a'] = chla[:]
                self.attrOtherProdUncDict['chlor_a_unc'] = {'long_name' : 'Uncertainty in ' +
                                                            chla.long_name,
                                                            '_FillValue':chla._FillValue,
                                                            'units': chla.units,
                                                            'valid_min': 0}
                self.dTypeDict['chlor_a_unc'] = chla.dtype
                self.dimsDict['chlor_a_unc'] = chla.dimensions
            if self.doNflh:
                nflh = geoVar['nflh']
                self.otherProdsDict['nflh'] = nflh[:]
                self.attrOtherProdUncDict['nflh_unc'] = {'long_name': 'Uncertainty in ' +
                                                    nflh.long_name,
                                                    '_FillValue': nflh._FillValue,
                                                    'units': nflh.units,
                                                    'scale_factor':nflh.scale_factor,
                                                    'add_offset':nflh.add_offset}
                self.dimsDict['nflh_unc'] = nflh.dimensions
                self.dTypeDict['nflh_unc'] = nflh.dtype

class MakeSwfUnc(MakeUnc):
    """Uncertainty subclass for SeaWiFS"""
    def __init__(self,*args,**kwargs):
        self.sensor ='SeaWiFS'
        self.bands = kwargs.pop("bands",
                                ['412','443','490','510','555','670','765','865'])
        self.colDict = {'412':'#001166','443':'#004488','490':'#116688',
                        '510':'#228844','555':'#667722','670':'#aa2211',
                        '765':'#770500','865':'#440000'}
        super(MakeSwfUnc,self).__init__(*args,**kwargs)

class MakeHMA(MakeUnc):
    """Uncertainty engine for HMODISA"""
    def __init__(self,*args,**kwargs):
        self.sensor='HMODISA'
        self.bands = kwargs.pop("bands",
                                ['412','443','488','531','547','555','645','667',
                                 '678','748','859','869','1240','1640','2130'])
        super(MakeHMA,self).__init__(*args,**kwargs)

def Main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--ifile', help='Initial L2 file path.',
                        type=str, required='True')
    parser.add_argument('-n', '--npath', help='Path to noisy data directory.',
                        type=str, required='True')
    parser.add_argument('-s', '--nsfx',
                        help='Noisy file suffix for pattern matching.',
                        type=str, default='_noisy_')
    parser.add_argument('-c', '--dochl', help='Compute chloropyll uncertainty.',
                        action='store_true')
    parser.add_argument('-f', '--doflh',
                        help='Compute normalized fluorescence line height.',
                        action='store_true')
    parser.add_argument('-v','--verbose',help='Augment output verbosity',
                        action='store_true')
    pArgs = parser.parse_args(argv)
    baseLineFile = pArgs.ifile
    noisyDataDir = pArgs.npath
    noisySfx = pArgs.nsfx
    baseLineFname = baseLineFile.split('/')[-1]
    if noisyDataDir[-1] != '/':
        noisyDataDir += '/'
    if baseLineFname[0] == 'S':
        uncObj = MakeSwfUnc(baseLineFile,noisyDataDir,verbose=pArgs.verbose)
    elif baseLineFname[0] == 'A':
        uncObj = MakeHMA(baseLineFile, noisyDataDir, doChla=pArgs.dochl,
                        doNflh=pArgs.doflh,verbose=pArgs.verbose)
    uncObj.ReadFromSilent()
    uncObj.BuildUncs(noisySfx,verbose=pArgs.verbose)
    uncObj.WriteToSilent()
    print("DONE!")

if __name__ == "__main__":

    Main(sys.argv[1:])
