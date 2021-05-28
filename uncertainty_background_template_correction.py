import matplotlib.pyplot as plt
import numpy as np

from params import ptcuttext, centralitycuttext
from template_fit import modifyBkgWeights
from utils import getNormHistAndErr, divideHistsAndErrs, getCenters, getBinRange


# call this first
def getDoubleRatioAndError(dfs, ptrange, isoParams, ssParams, centrange=None):
    isojjmcdf = dfs.fulljjmcdf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
    antiisojjmcdf = dfs.fulljjmcdf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))
    isodatadf = dfs.fulldatadf.query(isoParams.isocuttext()).query(ptcuttext(ptrange))
    antiisodatadf = dfs.fulldatadf.query(isoParams.antiisocuttext()).query(ptcuttext(ptrange))

    if centrange:
        isojjmcdf = isojjmcdf.query(centralitycuttext(centrange))
        antiisojjmcdf = antiisojjmcdf.query(centralitycuttext(centrange))
        isodatadf = isodatadf.query(centralitycuttext(centrange))
        antiisodatadf = antiisodatadf.query(centralitycuttext(centrange))

    isojjmchist, isojjmcerr = getNormHistAndErr(isojjmcdf, ssParams.ssvar, ssParams.doubleRatioBinEdges)
    antiisojjmchist, antiisojjmcerr = getNormHistAndErr(antiisojjmcdf, ssParams.ssvar, ssParams.doubleRatioBinEdges)
    isodatahist, isodataerr = getNormHistAndErr(isodatadf, ssParams.ssvar, ssParams.doubleRatioBinEdges)
    antiisodatahist, antiisodataerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.doubleRatioBinEdges)

    with np.errstate(divide='ignore', invalid='ignore'):
        jjmcratio, jjmcerr = divideHistsAndErrs(isojjmchist, isojjmcerr, antiisojjmchist, antiisojjmcerr)
        dataratio, dataerr = divideHistsAndErrs(isodatahist, isodataerr, antiisodatahist, antiisodataerr)
        doubleratio, doubleratioerr = divideHistsAndErrs(dataratio, dataerr, jjmcratio, jjmcerr)

    return doubleratio, doubleratioerr


# after getDoubleRatioAndError
# result gets saved into doubleRatioFit
# this performs the fit, which means it's slow!
def getDoubleRatioFitAndError(doubleratio, doubleratioerr, ssParams, doubleRatioFit, centrange=None):
    binCenters = getCenters(ssParams.doubleRatioBinEdges)
    fitRange = slice(*getBinRange(ssParams.doubleRatioBinEdges, *ssParams.doubleRatioFitRange))

    doubleRatioFit.getParamsAndErrors(doubleratio[fitRange], doubleratioerr[fitRange], binCenters[fitRange])


# after getDoubleRatioFitAndError has been called
def plotDoubleRatioAndFit(doubleratio, doubleratioerr, ssParams, doubleRatioFit, centrange=None):
    binCenters = getCenters(ssParams.doubleRatioBinEdges)
    fitRange = slice(*getBinRange(ssParams.doubleRatioBinEdges, *ssParams.doubleRatioFitRange))

    plt.errorbar(binCenters, doubleratio, yerr=doubleratioerr, fmt='ko')
    plt.plot(binCenters, list(map(doubleRatioFit.getFunctionFit(), binCenters)), 'r:')
    plt.plot(binCenters[fitRange], list(map(doubleRatioFit.getFunctionFit(), binCenters[fitRange])), 'r-')
    plt.annotate(doubleRatioFit.getResultText(), (0.95, 0.9), xycoords='axes fraction', ha='right', va='top', color='r')


# after getDoubleRatioFitAndError has been called
def getPuritiesWithDoubleRatio(dfs, ptrange, isoParams, ssParams, doubleRatioFit, centrange=None):
    baseBkgWeights = dfs.getBkgWeights(ptrange, isoParams, ssParams, centrange=centrange)
    purities = {}

    modFunction = doubleRatioFit.getFunctionFit()
    bkgWeights = modifyBkgWeights(baseBkgWeights, ssParams.binEdges, modFunction)
    tf = dfs.getTemplateFit(ptrange, isoParams, ssParams, bkgWeights=bkgWeights, centrange=centrange)
    purities['fit'] = tf.purity

    modFunction = doubleRatioFit.getFunctionErrHigh()
    bkgWeights = modifyBkgWeights(baseBkgWeights, ssParams.binEdges, modFunction)
    tf = dfs.getTemplateFit(ptrange, isoParams, ssParams, bkgWeights=bkgWeights, centrange=centrange)
    purities['high'] = tf.purity

    modFunction = doubleRatioFit.getFunctionErrLow()
    bkgWeights = modifyBkgWeights(baseBkgWeights, ssParams.binEdges, modFunction)
    tf = dfs.getTemplateFit(ptrange, isoParams, ssParams, bkgWeights=bkgWeights, centrange=centrange)
    purities['low'] = tf.purity

    return purities


def getBkgCorrectionUncertainty(dfs, ptrange, isoParams, ssParams, doubleRatioFit, centrange=None):
    getDoubleRatioFitAndError(dfs, ptrange, isoParams, ssParams, doubleRatioFit, centrange)
    purities = getPuritiesWithDoubleRatio(dfs, ptrange, isoParams, ssParams, doubleRatioFit, centrange)
    return np.ptp(list(purities.values())) / 2.0


def plotBkgCorrectionComp(dfs, ptrange, isoParams, ssParams, centrange=None):
    plt.subplot(121)
    tf = dfs.getTemplateFit(ptrange, isoParams, ssParams, centrange=centrange)
    handles = tf.plotFit()
    plt.legend(handles=handles)
    plt.title('No background template correction')

    plt.subplot(122)
    bkgWeights = dfs.getBkgWeights(ptrange, isoParams, ssParams)
    tf = dfs.getTemplateFit(ptrange, isoParams, ssParams, bkgWeights=bkgWeights, centrange=centrange)
    handles = tf.plotFit()
    plt.legend(handles=handles)
    plt.title('With background template correction')
