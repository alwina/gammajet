import glob
import os
import time

import numpy as np
import pandas as pd

from params import ptcuttext
from template_fit import applyBkgWeights, TemplateFit, BackgroundFit
from utils import getHistAndErr, getNormHistAndErr


def getDfFromCsv(filename):
    start = time.time()
    df = pd.read_csv(filename, '\t')
    end = time.time()
    print('Processed {0} in {1} seconds'.format(os.path.basename(filename), end - start))
    return df


def getDfFromCsvs(filenames):
    dfs = []
    for filename in filenames:
        dfs.append(getDfFromCsv(filename))
    return pd.concat(dfs).drop_duplicates().reset_index(drop=True)


def applyCuts(inputdf, cuts, verbose=True):
    df = inputdf
    for cut in cuts:
        df.query(cut, inplace=True)
        if verbose:
            print('{0}: {1}'.format(cut, df.shape[0]))
    return df


# To be built only from CSVs; creating the CSVs from the ntuples is a separate task
class DataframeCollection:
    def __init__(self, system, datafilepattern, gjmcfilepattern, jjmcfilepattern, csvdir='csv'):
        self.system = system
        self.csvdir = csvdir

        datafilenames = glob.glob(os.path.join(self.csvdir, datafilepattern))
        gjmcfilenames = glob.glob(os.path.join(self.csvdir, gjmcfilepattern))
        jjmcfilenames = glob.glob(os.path.join(self.csvdir, jjmcfilepattern))

        self.fulldatadf = getDfFromCsvs(datafilenames)
        self.fullgjmcdf = getDfFromCsvs(gjmcfilenames)
        self.fulljjmcdf = getDfFromCsvs(jjmcfilenames)

    def applyCuts(self, cuts):
        self.fulldatadf = applyCuts(self.fulldatadf, cuts['data'])
        self.fullgjmcdf = applyCuts(self.fullgjmcdf, cuts['gjmc'])
        self.fulljjmcdf = applyCuts(self.fulljjmcdf, cuts['jjmc'])

    def getTemplateFit(self, ptrange, isoParams, ssParams, bkgWeights=[], verbosity=0):
        isodatadf = self.fulldatadf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        isogjmcdf = self.fullgjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodatadf = self.fulldatadf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        data, dataerr = getHistAndErr(isodatadf, ssParams.ssvar, ssParams.binEdges)
        signal, signalerr = getNormHistAndErr(isogjmcdf, ssParams.ssvar, ssParams.binEdges)
        bkg, bkgerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.binEdges)

        if len(bkgWeights) > 0:
            bkg, bkgerr = applyBkgWeights(bkgWeights, bkg, bkgerr)

        return TemplateFit(data, dataerr, signal, signalerr, bkg, bkgerr, ssParams, verbosity)

    def getBackgroundFit(self, ptrange, isoParams, ssParams, bkgWeights=[], verbosity=0):
        isodatadf = self.fulldatadf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodatadf = self.fulldatadf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        data, dataerr = getHistAndErr(isodatadf, ssParams.ssvar, ssParams.binEdges)
        bkg, bkgerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.binEdges)

        if len(bkgWeights) > 0:
            bkg, bkgerr = applyBkgWeights(bkgWeights, bkg, bkgerr)

        return BackgroundFit(data, dataerr, bkg, bkgerr, ssParams, verbosity)

    # much cleaner to split off the calculation of the background weights,
    # even though it means there will always be this extra call
    def getBkgWeights(self, ptrange, isoParams, ssParams):
        isodf = self.fulljjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodf = self.fulljjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        isohist, _ = getNormHistAndErr(isodf, ssParams.ssvar, ssParams.binEdges)
        antiisohist, _ = getNormHistAndErr(antiisodf, ssParams.ssvar, ssParams.binEdges)

        weights = np.ones_like(isohist)
        np.divide(isohist, antiisohist, out=weights, where=antiisohist != 0)

        return weights
