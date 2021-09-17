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
            if 'Pb-Pb' in system:
                self.fullgjmcdf['weightswithraa'] = self.fullgjmcdf.apply(getWeightsWithRaa, axis=1)
        if jjmcfilenames:
            self.fulljjmcdf = getDfFromCsvs(jjmcfilenames)
            if 'Pb-Pb' in system:
                self.fulljjmcdf['weightswithraa'] = self.fulljjmcdf.apply(getWeightsWithRaa, axis=1)

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


raas = [[0.188494, 0.197833, 0.21251, 0.240708, 0.277338, 0.329186, 0.399066, 0.497627, 0.612843],
        [0.19801, 0.205901, 0.219408, 0.246369, 0.282159, 0.33265, 0.400528, 0.495272, 0.605806],
        [0.207792, 0.215552, 0.229416, 0.256647, 0.292653, 0.34283, 0.410216, 0.502372, 0.608506],
        [0.219327, 0.22716, 0.241571, 0.269387, 0.305382, 0.35592, 0.422556, 0.513105, 0.615657],
        [0.231176, 0.239203, 0.254035, 0.282325, 0.318783, 0.369641, 0.435871, 0.525221, 0.624577],
        [0.243387, 0.251714, 0.266929, 0.295954, 0.332856, 0.384008, 0.449979, 0.537574, 0.634146],
        [0.255471, 0.264083, 0.279983, 0.309645, 0.347274, 0.398806, 0.464663, 0.551295, 0.645121],
        [0.267528, 0.276849, 0.293235, 0.323968, 0.362188, 0.414036, 0.480073, 0.565651, 0.656529],
        [0.279677, 0.289361, 0.306376, 0.337932, 0.377019, 0.428968, 0.494085, 0.579171, 0.668245],
        [0.291486, 0.30143, 0.319046, 0.351662, 0.391104, 0.443572, 0.508926, 0.592154, 0.67925],
        [0.302886, 0.313222, 0.33148, 0.364769, 0.404731, 0.457459, 0.521857, 0.603834, 0.688718],
        [0.314375, 0.324846, 0.343726, 0.377508, 0.417996, 0.470983, 0.535021, 0.615801, 0.698341],
        [0.324738, 0.335729, 0.355079, 0.389842, 0.43077, 0.483672, 0.547258, 0.626229, 0.707458],
        [0.334996, 0.346594, 0.366417, 0.401538, 0.443014, 0.495936, 0.559099, 0.637065, 0.715108],
        [0.34458, 0.356343, 0.376649, 0.413197, 0.454417, 0.507353, 0.570194, 0.646609, 0.722773],
        [0.35326, 0.365369, 0.386718, 0.423485, 0.465517, 0.517916, 0.580194, 0.655166, 0.728586],
        [0.361728, 0.374677, 0.39586, 0.433284, 0.475294, 0.527615, 0.588766, 0.662121, 0.733835],
        [0.37243, 0.385905, 0.408434, 0.446277, 0.488134, 0.540125, 0.599592, 0.669964, 0.738654],
        [0.386406, 0.40084, 0.42425, 0.462491, 0.504584, 0.555529, 0.612746, 0.680579, 0.747146],
        [0.402806, 0.418223, 0.442969, 0.482594, 0.524538, 0.573984, 0.628629, 0.692367, 0.752749],
        [0.419814, 0.436492, 0.462543, 0.503165, 0.546383, 0.593386, 0.645866, 0.705326, 0.76207],
        [0.431153, 0.448903, 0.476458, 0.517984, 0.560118, 0.605492, 0.65457, 0.708915, 0.764124],
        [0.434408, 0.454063, 0.482069, 0.523636, 0.565777, 0.606436, 0.65387, 0.706258, 0.763671],
        [0.435141, 0.455392, 0.484958, 0.524784, 0.5656, 0.606101, 0.653692, 0.704346, 0.763553],
        [0.427501, 0.447811, 0.477647, 0.517433, 0.559215, 0.597926, 0.645374, 0.697378, 0.756974],
        [0.413536, 0.433489, 0.465074, 0.505704, 0.547201, 0.585693, 0.634513, 0.687034, 0.750963],
        [0.393387, 0.416557, 0.44601, 0.488185, 0.531349, 0.57093, 0.622151, 0.675746, 0.743261],
        [0.370042, 0.393088, 0.426042, 0.466699, 0.512259, 0.553493, 0.608739, 0.667277, 0.735648],
        [0.347612, 0.369306, 0.404447, 0.446322, 0.492533, 0.536444, 0.593985, 0.652662, 0.731987],
        [0.312928, 0.337618, 0.370519, 0.413696, 0.464106, 0.512405, 0.57695, 0.646534, 0.732126],
        [0.262523, 0.288506, 0.322048, 0.368801, 0.422067, 0.473717, 0.545127, 0.626515, 0.714316],
        [0.205135, 0.23011, 0.264207, 0.315144, 0.373652, 0.43509, 0.517932, 0.602702, 0.706452],
        [0.149629, 0.17454, 0.207704, 0.259074, 0.322556, 0.39306, 0.487361, 0.588318, 0.703458],
        [0.131015, 0.153942, 0.186519, 0.240468, 0.304813, 0.385128, 0.477856, 0.582286, 0.696775],
        [0.139377, 0.163514, 0.198945, 0.255147, 0.323724, 0.40547, 0.509364, 0.618491, 0.726998],
        [0.161794, 0.187576, 0.22437, 0.284281, 0.359381, 0.447792, 0.547126, 0.666463, 0.735024],
        [0.20074, 0.228634, 0.273301, 0.338241, 0.414265, 0.505078, 0.593209, 0.673702, 0.748473],
        [0.294104, 0.324423, 0.371913, 0.451159, 0.536475, 0.624313, 0.691564, 0.793201, 0.860148],
        [0.364925, 0.421362, 0.482606, 0.577559, 0.615381, 0.765339, 0.755291, 0.826902, 1.001153]]


def getWeightsWithRaa(row):
    mcweight = row['weights']
    parentpi0pt = row['cluster_mc_parentpi0pt']
    centrality = row['centrality_v0m']

    # if this is the dummy value, use the weight
    if parentpi0pt == -1:
        return mcweight

    # for now, use RAA of charged particles because I can't find pi0 RAA
    if 0 < centrality <= 5:
        centindex = 0
    elif 5 < centrality <= 80:
        centindex = int(math.ceil(centrality / 10.0))
    else:
        centindex = 8

    if 0.15 < parentpi0pt <= 1.00:
        ptindex = int(math.ceil(parentpi0pt / 0.05) - 4)
    elif 1.00 < parentpi0pt <= 1.10:
        ptindex = 17
    elif 1.10 < parentpi0pt <= 1.20:
        ptindex = 18
    elif 1.20 < parentpi0pt <= 3.20:
        ptindex = int(math.ceil(parentpi0pt / 0.20) + 12)
    elif 3.20 < parentpi0pt <= 3.60:
        ptindex = 29
    elif 3.60 < parentpi0pt <= 4.00:
        ptindex = 30
    elif 4.00 < parentpi0pt <= 5.00:
        ptindex = 31
    elif 5.00 < parentpi0pt <= 6.00:
        ptindex = 32
    elif 6.00 < parentpi0pt <= 8.00:
        ptindex = 33
    elif 8.00 < parentpi0pt <= 10.00:
        ptindex = 34
    elif 10 < parentpi0pt <= 13:
        ptindex = 35
    elif 13 < parentpi0pt <= 20:
        ptindex = 36
    elif 20 < parentpi0pt <= 30:
        ptindex = 37
    elif 30 < parentpi0pt <= 50:
        ptindex = 38
    else:
        if parentpi0pt < 0.15:
            ptindex = 0
        elif parentpi0pt > 50:
            ptindex = 38

    raa = raas[ptindex][centindex]

    return mcweight * raa
