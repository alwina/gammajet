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
    return binmin, binmax


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
