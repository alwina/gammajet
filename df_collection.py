import glob
import os
import time

import numpy as np
import pandas as pd

from params import ptcuttext, centralitycuttext
from template_fit import applyBkgWeights, TemplateFit, BackgroundFit
from utils import getHistAndErr, getNormHistAndErr, divideHistsAndErrs


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


def parseCut(cutvar, cuttype, cutval):
    if cuttype == 'min':
        cuttext = '{0} > {1}'.format(cutvar, cutval)
    elif cuttype == 'max':
        cuttext = '{0} < {1}'.format(cutvar, cutval)
    elif cuttype == 'incmin':
        cuttext = '{0} >= {1}'.format(cutvar, cutval)
    elif cuttype == 'incmax':
        cuttext = '{0} <= {1}'.format(cutvar, cutval)
    elif cuttype == 'equals':
        cuttext = '{0} == {1}'.format(cutvar, cutval)
    else:
        print('Cut type {0} not recognized'.format(cutval))
        cuttext = ''

    return cuttext


def getDataframeCollection(config):
    datacuts = []
    gjmccuts = []
    jjmccuts = []

    for cutvar, cutvals in config['clustercuts']['all'].iteritems():
        for cuttype, cutval in cutvals.iteritems():
            cuttext = parseCut(cutvar, cuttype, cutval)

            datacuts.append(cuttext)
            gjmccuts.append(cuttext)
            jjmccuts.append(cuttext)

    if 'data' in config['clustercuts']:
        for cutvar, cutvals in config['clustercuts']['data'].iteritems():
            for cuttype, cutval in cutvals.iteritems():
                cuttext = parseCut(cutvar, cuttype, cutval)
                datacuts.append(cuttext)

    if 'gjmc' in config['clustercuts']:
        for cutvar, cutvals in config['clustercuts']['gjmc'].iteritems():
            for cuttype, cutval in cutvals.iteritems():
                cuttext = parseCut(cutvar, cuttype, cutval)
                gjmccuts.append(cuttext)

    if 'jjmc' in config['clustercuts']:
        for cutvar, cutvals in config['clustercuts']['jjmc'].iteritems():
            for cuttype, cutval in cutvals.iteritems():
                cuttext = parseCut(cutvar, cuttype, cutval)
                jjmccuts.append(cuttext)

    cuts = {}
    cuts['data'] = datacuts
    cuts['gjmc'] = gjmccuts
    cuts['jjmc'] = jjmccuts

    dfs = DataframeCollection(config['longname'], config['filelists']['clustercsvs']['data'], config['filelists']['clustercsvs']['gjmc'], config['filelists']['clustercsvs']['jjmc'])
    dfs.applyCuts(cuts)

    return dfs


# To be built only from CSVs; creating the CSVs from the ntuples is a separate task
class DataframeCollection:
    def __init__(self, system='', datafilepattern='', gjmcfilepattern='', jjmcfilepattern='', csvdir=''):
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

    def getTemplateFit(self, ptrange, isoParams, ssParams, bkgWeights=[], verbosity=0, centrange=None, additionalCuts={}):
        isodatadf = self.fulldatadf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        isogjmcdf = self.fullgjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodatadf = self.fulldatadf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        if centrange:
            isodatadf = isodatadf.query(centralitycuttext(centrange))
            isogjmcdf = isogjmcdf.query(centralitycuttext(centrange))
            antiisodatadf = antiisodatadf.query(centralitycuttext(centrange))

        if 'data' in additionalCuts:
            for cut in additionalCuts['data']:
                isodatadf = isodatadf.query(cut)
                antiisodatadf = antiisodatadf.query(cut)

        if 'gjmc' in additionalCuts:
            for cut in additionalCuts['gjmc']:
                isogjmcdf = isogjmcdf.query(cut)

        data, dataerr = getHistAndErr(isodatadf, ssParams.ssvar, ssParams.binEdges)
        signal, signalerr = getNormHistAndErr(isogjmcdf, ssParams.ssvar, ssParams.binEdges)
        bkg, bkgerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.binEdges)

        if len(bkgWeights) > 0:
            bkg, bkgerr = applyBkgWeights(bkgWeights, bkg, bkgerr)

        return TemplateFit(data, dataerr, signal, signalerr, bkg, bkgerr, ssParams, verbosity)

    def getBackgroundFit(self, ptrange, isoParams, ssParams, bkgWeights=[], verbosity=0, centrange=None, additionalCuts={}):
        isodatadf = self.fulldatadf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodatadf = self.fulldatadf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        if centrange:
            isodatadf = isodatadf.query(centralitycuttext(centrange))
            antiisodatadf = antiisodatadf.query(centralitycuttext(centrange))

        if 'data' in additionalCuts:
            for cut in additionalCuts['data']:
                isodatadf = isodatadf.query(cut)
                antiisodatadf = antiisodatadf.query(cut)

        data, dataerr = getHistAndErr(isodatadf, ssParams.ssvar, ssParams.binEdges)
        bkg, bkgerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.binEdges)

        if len(bkgWeights) > 0:
            bkg, bkgerr = applyBkgWeights(bkgWeights, bkg, bkgerr)

        return BackgroundFit(data, dataerr, bkg, bkgerr, ssParams, verbosity)

    # much cleaner to split off the calculation of the background weights,
    # even though it means there will always be this extra call
    def getBkgWeights(self, ptrange, isoParams, ssParams, centrange=None, useraa=True, additionalCuts={}):
        isojjmcdf = self.fulljjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisogjmcdf = self.fullgjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))
        antiisojjmcdf = self.fulljjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        if centrange:
            isojjmcdf = isojjmcdf.query(centralitycuttext(centrange))
            antiisogjmcdf = antiisogjmcdf.query(centralitycuttext(centrange))
            antiisojjmcdf = antiisojjmcdf.query(centralitycuttext(centrange))

        if 'gjmc' in additionalCuts:
            for cut in additionalCuts['gjmc']:
                antiisogjmcdf = antiisogjmcdf.query(cut)

        if 'jjmc' in additionalCuts:
            for cut in additionalCuts['jjmc']:
                isojjmcdf = isojjmcdf.query(cut)
                antiisojjmcdf = antiisojjmcdf.query(cut)

        if useraa:
            antiisomcdf = pd.concat([antiisogjmcdf, antiisojjmcdf])
        else:
            antiisomcdf = antiisojjmcdf

        isohist, _ = getNormHistAndErr(isojjmcdf, ssParams.ssvar, ssParams.binEdges)
        antiisohist, _ = getNormHistAndErr(antiisomcdf, ssParams.ssvar, ssParams.binEdges)

        weights = np.ones_like(isohist)
        np.divide(isohist, antiisohist, out=weights, where=antiisohist != 0)

        return weights

    def getBkgWeightsAndErrs(self, ptrange, isoParams, ssParams, centrange=None, useraa=True, additionalCuts={}):
        isojjmcdf = self.fulljjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisogjmcdf = self.fullgjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))
        antiisojjmcdf = self.fulljjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        if centrange:
            isojjmcdf = isojjmcdf.query(centralitycuttext(centrange))
            antiisogjmcdf = antiisogjmcdf.query(centralitycuttext(centrange))
            antiisojjmcdf = antiisojjmcdf.query(centralitycuttext(centrange))

        if 'gjmc' in additionalCuts:
            for cut in additionalCuts['gjmc']:
                antiisogjmcdf = antiisogjmcdf.query(cut)

        if 'jjmc' in additionalCuts:
            for cut in additionalCuts['jjmc']:
                isojjmcdf = isojjmcdf.query(cut)
                antiisojjmcdf = antiisojjmcdf.query(cut)

        if useraa:
            antiisomcdf = pd.concat([antiisogjmcdf, antiisojjmcdf])
            weightvar = 'weightswithraa'
        else:
            antiisomcdf = antiisojjmcdf
            weightvar = 'weights'

        isohist, isoerr = getNormHistAndErr(isojjmcdf, ssParams.ssvar, ssParams.binEdges, weightvar=weightvar)
        antiisohist, antiisoerr = getNormHistAndErr(antiisomcdf, ssParams.ssvar, ssParams.binEdges, weightvar=weightvar)

        weights, weightserr = divideHistsAndErrs(isohist, isoerr, antiisohist, antiisoerr)

        return weights, weightserr

    def getDoubleRatioAndError(self, ptrange, isoParams, ssParams, centrange=None, useraa=True, additionalCuts={}):
        isodatadf = self.fulldatadf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisodatadf = self.fulldatadf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))
        isogjmcdf = self.fullgjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisogjmcdf = self.fullgjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))
        isojjmcdf = self.fulljjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
        antiisojjmcdf = self.fulljjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

        if centrange:
            isodatadf = isodatadf.query(centralitycuttext(centrange))
            antiisodatadf = antiisodatadf.query(centralitycuttext(centrange))
            isogjmcdf = isogjmcdf.query(centralitycuttext(centrange))
            antiisogjmcdf = antiisogjmcdf.query(centralitycuttext(centrange))
            isojjmcdf = isojjmcdf.query(centralitycuttext(centrange))
            antiisojjmcdf = antiisojjmcdf.query(centralitycuttext(centrange))

        if 'data' in additionalCuts:
            for cut in additionalCuts:
                isodatadf = isodatadf.query(cut)
                antiisodatadf = antiisodatadf.query(cut)

        if 'gjmc' in additionalCuts:
            for cut in additionalCuts:
                isogjmcdf = isogjmcdf.query(cut)
                antiisogjmcdf = antiisogjmcdf.query(cut)

        if 'jjmc' in additionalCuts:
            for cut in additionalCuts:
                isojjmcdf = isojjmcdf.query(cut)
                antiisojjmcdf = antiisojjmcdf.query(cut)

        if useraa:
            isomcdf = pd.concat([isojjmcdf, isogjmcdf])
            antiisomcdf = pd.concat([antiisojjmcdf, antiisogjmcdf])
            weightvar = 'weightswithraa'
        else:
            isomcdf = isojjmcdf
            antiisomcdf = antiisojjmcdf
            weightvar = 'weights'

        isomchist, isomcerr = getNormHistAndErr(isomcdf, ssParams.ssvar, ssParams.doubleRatioBinEdges, weightvar=weightvar)
        antiisomchist, antiisomcerr = getNormHistAndErr(antiisomcdf, ssParams.ssvar, ssParams.doubleRatioBinEdges, weightvar=weightvar)

        isodatahist, isodataerr = getNormHistAndErr(isodatadf, ssParams.ssvar, ssParams.doubleRatioBinEdges)
        antiisodatahist, antiisodataerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.doubleRatioBinEdges)

        with np.errstate(divide='ignore', invalid='ignore'):
            jjmcratio, jjmcerr = divideHistsAndErrs(isomchist, isomcerr, antiisomchist, antiisomcerr)
            dataratio, dataerr = divideHistsAndErrs(isodatahist, isodataerr, antiisodatahist, antiisodataerr)
            doubleratio, doubleratioerr = divideHistsAndErrs(dataratio, dataerr, jjmcratio, jjmcerr)

        return doubleratio, doubleratioerr
