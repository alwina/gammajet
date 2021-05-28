import matplotlib.pyplot as plt
import numpy as np


def haveSameLength(*args):
    n = len(args[0])
    return all(len(array) == n for array in args)


def is1D(x):
    return len(np.array(x).shape) == 1


def quadSumAll(x):
    if not is1D(x):
        return TypeError('Input should be 1D')
    else:
        return np.sqrt(np.sum(np.square(x)))


def quadSumPairwise(x, y):
    if not is1D(x) or not is1D(y):
        raise TypeError('Both inputs should be 1D')
    elif not haveSameLength(x, y):
        raise TypeError('Inputs should have the same length')
    else:
        return np.sqrt(np.sum(np.square([x, y]), axis=0))


def getRanges(rangesOrEdges):
    if len(np.array(rangesOrEdges).shape) == 1:
        ranges = list(zip(rangesOrEdges[:-1], rangesOrEdges[1:]))
    elif len(np.array(rangesOrEdges).shape) == 2:
        ranges = rangesOrEdges
    else:
        raise TypeError('Input should be edges or ranges')
    return ranges


def getCenters(rangesOrEdges):
    return np.mean(getRanges(rangesOrEdges), axis=1)


def getWidths(rangesOrEdges):
    return np.ptp(getRanges(rangesOrEdges), axis=1)


def getXerrForPlot(rangesOrEdges):
    return getWidths(rangesOrEdges) / 2.0


def getBinRange(binEdges, valuemin, valuemax):
    binmin = min([i for i, edge in enumerate(binEdges) if edge >= valuemin])
    binmax = max([i for i, edge in enumerate(binEdges) if edge <= valuemax])
    return binmin, binmax + 1


def getErr(dist, disterr):
    if disterr is None:
        return np.sqrt(dist)
    else:
        return np.array(disterr, dtype='f')


def getHistAndErr(df, var, binEdges, weightvar='weights'):
    hist = np.histogram(df[var], bins=binEdges, weights=df[weightvar])[0]
    errs = []
    for left, right in zip(binEdges[:-1], binEdges[1:]):
        weights = [w for (x, w) in zip(df[var], df[weightvar]) if x >= left and x <= right]
        errs.append(quadSumAll(weights))
    binWidths = getWidths(binEdges)
    return np.divide(hist, binWidths), np.divide(errs, binWidths)


def normalizeHistAndErr(hist, err):
    norm = float(np.sum(hist))
    return np.divide(hist, norm), np.divide(err, norm)


def getNormHistAndErr(df, var, binEdges, weightvar='weights'):
    hist, err = getHistAndErr(df, var, binEdges, weightvar)
    return normalizeHistAndErr(hist, err)


def divideHistsAndErrs(hist1, err1, hist2, err2):
    hist = np.full_like(hist1, np.nan)
    np.divide(hist1, hist2, out=hist, where=hist2 != 0)

    relerr1 = np.zeros_like(hist1)
    np.divide(err1, hist1, out=relerr1, where=err1 != 0)

    relerr2 = np.zeros_like(hist2)
    np.divide(err2, hist2, out=relerr2, where=err2 != 0)

    relativeerr = quadSumPairwise(relerr1, relerr2)
    histerr = np.multiply(hist, relativeerr)

    return hist, histerr


def integrateHist(hist, binEdges):
    return np.sum(np.multiply(hist, getWidths(binEdges)))


def getUniformUncertainty(dist):
    return np.ptp(dist) / np.sqrt(12)


# skip the underflow and overflow bins for now -- we can change our minds later
# these assume equally-spaced bins

def th1ToArrays(th1):
    hist, err, binCenters, binWidths = [], [], [], []
    for binX in range(1, th1.GetNbinsX() + 1):
        hist.append(th1.GetBinContent(binX))
        err.append(th1.GetBinError(binX))
        binCenters.append(th1.GetBinCenter(binX))
        binWidths.append(th1.GetBinWidth(binX))
    return hist, err, binCenters, binWidths


def plotTH1(th1, **kwargs):
    hist, err, binCenters, binWidths = th1ToArrays(th1)
    plt.errorbar(binCenters, hist, yerr=err, xerr=binWidths, **kwargs)


def plotTH2(th2):
    nBinsX = th2.GetNbinsX()
    nBinsY = th2.GetNbinsY()

    hist2d = np.zeros((nBinsX + 1, nBinsY + 1))
    for binX in range(1, nBinsX + 1):
        for binY in range(1, nBinsY + 1):
            globalBinNumber = th2.GetBin(binX, binY)
            hist2d[binX][binY] = th2.GetBinContent(globalBinNumber)
    # set 0 values to be white when plotting
    hist2d[hist2d == 0] = np.nan

    minX = th2.GetXaxis().GetXmin()
    maxX = th2.GetXaxis().GetXmax()
    minY = th2.GetYaxis().GetXmin()
    maxY = th2.GetYaxis().GetXmax()
    plt.imshow(hist2d, origin='lower', extent=(minX, maxX, minY, maxY))


def sliceAndProjectTHnSparse(thnSparse, slices, *axesToProject):
    for (axisNumber, axisMin, axisMax) in slices:
        thnSparse.GetAxis(axisNumber).SetRangeUser(axisMin, axisMax)
    return thnSparse.Projection(*axesToProject)


def Chi2Histograms(df1, df2, var, bins):
    hist1, err1 = getNormHistAndErr(df1, var, bins)
    hist2, err2 = getNormHistAndErr(df2, var, bins)

    toterr = quadSumPairwise(err1, err2)

    numerator = np.square(np.subtract(hist1, hist2))
    denominator = np.square(toterr)

    chi2 = np.sum(np.divide(numerator, denominator))
    return chi2


def getTimeText(duration):
    if duration < 120:
        timetext = '{0:0.0f} seconds'.format(duration)
    elif duration < 3600:
        timetext = '{0:0.0f} minutes'.format(duration / 60.0)
    else:
        timetext = '{0:0.0f} hours {1:0.0f} minutes'.format(duration / 3600, (duration % 3600) / 60.0)
    return timetext
