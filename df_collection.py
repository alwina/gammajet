import glob
import os
import time

import numpy as np
import pandas as pd

from params import ptcuttext, centralitycuttext
from template_fit import applyBkgWeights, TemplateFit, BackgroundFit
from utils import getHistAndErr, getNormHistAndErr


def getDfFromCsv(filename):
    start = time.time()
    df = pd.read_csv(filename, '\t')
    end = time.time()
    print('Processed {0} in {1:0.1f} seconds'.format(os.path.basename(filename), end - start))
    return df


def getDfFromCsvs(filenames):
    dfs = []
    for filename in filenames:
        dfs.append(getDfFromCsv(filename))
    return pd.concat(dfs).drop_duplicates().reset_index(drop=True)


def applyCuts(inputdf, cuts, verbose=True):
    df = inputdf
    if verbose:
        print('No cuts: {0}'.format(df.shape[0]))
    for cut in cuts:
        df = df.query(cut)
        if verbose:
            print('{0}: {1}'.format(cut, df.shape[0]))
    return df


# To be built only from CSVs; creating the CSVs from the ntuples is a separate task
class DataframeCollection:
    def __init__(self, system='', datafilepattern='', gjmcfilepattern='', jjmcfilepattern='', csvdir='csv'):
        self.system = system

        datafilenames = []
        gjmcfilenames = []
        jjmcfilenames = []

        if datafilepattern:
            datafilenames = glob.glob(os.path.join(csvdir, datafilepattern))
        if gjmcfilepattern:
            gjmcfilenames = glob.glob(os.path.join(csvdir, gjmcfilepattern))
        if jjmcfilepattern:
            jjmcfilenames = glob.glob(os.path.join(csvdir, jjmcfilepattern))

        if datafilenames:
            self.fulldatadf = getDfFromCsvs(datafilenames)
        if gjmcfilenames:
            self.fullgjmcdf = getDfFromCsvs(gjmcfilenames)
        if jjmcfilenames:
            self.fulljjmcdf = getDfFromCsvs(jjmcfilenames)

    def copy(self, inputDfCollection):
        self.system = inputDfCollection.system
        self.fulldatadf = inputDfCollection.fulldatadf
        self.fullgjmcdf = inputDfCollection.fullgjmcdf
        self.fulljjmcdf = inputDfCollection.fulljjmcdf

    def applyCuts(self, cuts):
        if 'data' in cuts:
            print('Data')
            self.fulldatadf = applyCuts(self.fulldatadf, cuts['data'])
        if 'gjmc' in cuts:
            print('GJ MC')
            self.fullgjmcdf = applyCuts(self.fullgjmcdf, cuts['gjmc'])
        if 'jjmc' in cuts:
            print('JJ MC')
            self.fulljjmcdf = applyCuts(self.fulljjmcdf, cuts['jjmc'])

    def getTemplateFit(self, ptrange, isoParams, ssParams, bkgWeights=[], verbosity=0, centrange=None):
        isodatadf = self.fulldatadf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        isogjmcdf = self.fullgjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodatadf = self.fulldatadf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        if centrange:
            isodatadf = isodatadf.query(centralitycuttext(centrange))
            isogjmcdf = isogjmcdf.query(centralitycuttext(centrange))
            antiisodatadf = antiisodatadf.query(centralitycuttext(centrange))

        data, dataerr = getHistAndErr(isodatadf, ssParams.ssvar, ssParams.binEdges)
        signal, signalerr = getNormHistAndErr(isogjmcdf, ssParams.ssvar, ssParams.binEdges)
        bkg, bkgerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.binEdges)

        if len(bkgWeights) > 0:
            bkg, bkgerr = applyBkgWeights(bkgWeights, bkg, bkgerr)

        return TemplateFit(data, dataerr, signal, signalerr, bkg, bkgerr, ssParams, verbosity)

    def getBackgroundFit(self, ptrange, isoParams, ssParams, bkgWeights=[], verbosity=0, centrange=None):
        isodatadf = self.fulldatadf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodatadf = self.fulldatadf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        if centrange:
            isodatadf = isodatadf.query(centralitycuttext(centrange))
            antiisodatadf = antiisodatadf.query(centralitycuttext(centrange))

        data, dataerr = getHistAndErr(isodatadf, ssParams.ssvar, ssParams.binEdges)
        bkg, bkgerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.binEdges)

        if len(bkgWeights) > 0:
            bkg, bkgerr = applyBkgWeights(bkgWeights, bkg, bkgerr)

        return BackgroundFit(data, dataerr, bkg, bkgerr, ssParams, verbosity)

    # much cleaner to split off the calculation of the background weights,
    # even though it means there will always be this extra call
    def getBkgWeights(self, ptrange, isoParams, ssParams, centrange=None):
        isodf = self.fulljjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodf = self.fulljjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        if centrange:
            isodf = isodf.query(centralitycuttext(centrange))
            antiisodf = antiisodf.query(centralitycuttext(centrange))

        isohist, _ = getNormHistAndErr(isodf, ssParams.ssvar, ssParams.binEdges)
        antiisohist, _ = getNormHistAndErr(antiisodf, ssParams.ssvar, ssParams.binEdges)

        weights = np.ones_like(isohist)
        np.divide(isohist, antiisohist, out=weights, where=antiisohist != 0)

        return weights
