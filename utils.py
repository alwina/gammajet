import matplotlib.pyplot as plt
import numpy as np
import struct


def haveSameLength(*args):
    """
    Checks whether the input arrays have the same length
    Each array is its own argument, and there can (in principle) be an arbitrary number of them
    """
    n = len(args[0])
    return all(len(arr) == n for arr in args)


def is1D(x):
    """
    Checks whether the input array is one-dimensional
    """
    return len(np.array(x).shape) == 1


def quadSumAll(x):
    """
    Does the quadruture sum (square root of the sum of squares) of 1D input array
    """
    if not is1D(x):
        return TypeError('Input should be 1D')
    else:
        return np.sqrt(np.sum(np.square(x)))


def quadSumPairwise(x, y):
    """
    Does the quadruture sum pairwise (element-by-element) of two 1D arrays
    """
    if not is1D(x) or not is1D(y):
        raise TypeError('Both inputs should be 1D')
    elif not haveSameLength(x, y):
        raise TypeError('Inputs should have the same length')
    else:
        return np.sqrt(np.sum(np.square([x, y]), axis=0))


def getRanges(rangesOrEdges):
    """
    Generates ranges (list of tuples, e.g. [(1, 2), (2, 3)]) from either ranges or edges (1D list, e.g. [1, 2, 3])

    Mostly used as an intermediate for things like getCenters, getWidths, etc
    """
    if len(np.array(rangesOrEdges).shape) == 1:
        ranges = list(zip(rangesOrEdges[:-1], rangesOrEdges[1:]))
    elif len(np.array(rangesOrEdges).shape) == 2:
        ranges = rangesOrEdges
    else:
        raise TypeError('Input should be edges or ranges')
    return ranges


def getCenters(rangesOrEdges):
    """
    Gets the centers of the bins, either from ranges (list of tuple pairs) or edges
    """
    return np.mean(getRanges(rangesOrEdges), axis=1)


def getWidths(rangesOrEdges):
    """
    Gets the widths of the bins, either from ranges (list of tuple pairs) or edges
    """
    return np.ptp(getRanges(rangesOrEdges), axis=1)


def getXerrForPlot(rangesOrEdges):
    """
    Gets half the bin widths which is what pyplot wants for the xerr when plotting errorbars, either from ranges (list of tuple pairs) or edges
    """
    return getWidths(rangesOrEdges) / 2.0


def getBinRange(binEdges, valuemin, valuemax):
    """
    Gets the bin indices for a slice corresponding to the range (valuemin, valuemax). Takes only bins that are fully within the input range.
    """
    binmin = min([i for i, edge in enumerate(binEdges) if edge >= valuemin])
    binmax = max([i for i, edge in enumerate(binEdges) if edge <= valuemax])
    return binmin, binmax


def getErr(dist, disterr):
    """
    Calculates the sqrt(N) errors on a distribution if they don't exist; otherwise, converts them to a numpy array of floats
    """
    if disterr is None:
        return np.sqrt(dist)
    else:
        return np.array(disterr, dtype='f')


def getHistAndErr(inputdf, var, binEdges, weightvar='weights'):
    """
    Creates and returns a histogram and its errors. The bin widths are NOT divided out.

    inputdf: a pandas dataframe or equivalent structure
    var: the variable to histogram
    binEdges: the edges of the bins of the histogram
    weightvar: the variable to use for the weights of the histogram; this is not optional
    """
    # handle the situation where the weights might be NaN for whatever reason
    df = inputdf.query('{0}=={0}'.format(weightvar))
    hist = np.histogram(df[var], bins=binEdges, weights=df[weightvar])[0]
    errs = []
    for left, right in zip(binEdges[:-1], binEdges[1:]):
        weights = [w for (x, w) in zip(df[var], df[weightvar]) if x >= left and x < right]
        errs.append(quadSumAll(weights))
    return hist, errs


def normalizeHistAndErr(hist, err):
    """
    Normalizes a histogram and its errors such that it sums to 1
    """
    norm = float(np.sum(hist))
    return np.divide(hist, norm), np.divide(err, norm)


def getNormHistAndErr(df, var, binEdges, weightvar='weights'):
    """
    Creates and returns a histogram and its errors such that it sums to 1
    """
    hist, err = getHistAndErr(df, var, binEdges, weightvar)
    return normalizeHistAndErr(hist, err)


def divideHistsAndErrs(hist1, err1, hist2, err2):
    """
    Divides two histograms and combines their uncertainties correctly. If there's a divide-by-zero issue, the result in that bin will be NaN.
    """

    # start with NaN so that any 0s in hist2 give NaN
    hist = np.full_like(hist1, np.nan)
    np.divide(hist1, hist2, out=hist, where=hist2 != 0)

    # start with 0s
    relerr1 = np.zeros_like(hist1)
    np.divide(err1, hist1, out=relerr1, where=err1 != 0)

    # start with 0s
    relerr2 = np.zeros_like(hist2)
    np.divide(err2, hist2, out=relerr2, where=err2 != 0)

    relativeerr = quadSumPairwise(relerr1, relerr2)
    histerr = np.multiply(hist, relativeerr)

    return hist, histerr


def integrateHist(hist, binEdges):
    """
    Integrates a histogram, i.e. summing after multiplying by bin widths
    """
    return np.sum(np.multiply(hist, getWidths(binEdges)))


def getUniformUncertainty(dist):
    """
    Gets the uncertainty assuming a uniform distribution
    """
    return np.ptp(dist) / np.sqrt(12)


# skip the underflow and overflow bins for now -- we can change our minds later
# these assume equally-spaced bins
def th1ToArrays(th1):
    """
    Converts a ROOT TH1 to numpy arrays (hist, err, binCenters, binWidths)
    """
    hist, err, binCenters, binWidths = [], [], [], []
    for binX in range(1, th1.GetNbinsX() + 1):
        hist.append(th1.GetBinContent(binX))
        err.append(th1.GetBinError(binX))
        binCenters.append(th1.GetBinCenter(binX))
        binWidths.append(th1.GetBinWidth(binX))
    return np.array(hist), np.array(err), np.array(binCenters), np.array(binWidths)


def plotTH1(th1, **kwargs):
    """
    Plots a ROOT TH1 with pyplot
    """
    hist, err, binCenters, binWidths = th1ToArrays(th1)
    plt.errorbar(binCenters, hist, yerr=err, xerr=np.divide(binWidths, 2.0), **kwargs)


def plotTH2(th2, **kwargs):
    """
    Plots a ROOT TH2 with pyplot. No colorbar plotted; use plt.colorbar() after to show that. If colorbar is desired, recommend a figure size like (10, 8) to have a square plot.
    """
    nBinsX = th2.GetNbinsX()
    nBinsY = th2.GetNbinsY()

    hist2d = np.zeros((nBinsY, nBinsX))
    for binX in range(1, nBinsX + 1):
        for binY in range(1, nBinsY + 1):
            globalBinNumber = th2.GetBin(binX, binY)
            hist2d[binY - 1][binX - 1] = th2.GetBinContent(globalBinNumber)
    # set 0 values to be white when plotting
    hist2d[hist2d == 0] = np.nan

    minX = th2.GetXaxis().GetXmin()
    maxX = th2.GetXaxis().GetXmax()
    minY = th2.GetYaxis().GetXmin()
    maxY = th2.GetYaxis().GetXmax()
    plt.imshow(hist2d, origin='lower', extent=(minX, maxX, minY, maxY), **kwargs)


def sliceAndProjectTHnSparse(thnSparse, slices, *axesToProject):
    """
    Makes cuts on a THnSparse and returns a TH* projection

    thnSparse: the THnSparse to manipulate
    slices: cuts in the form of a list of tuples of the form (axisNumber, axisMin, axisMax)
    axesToProject: the axes numbers to pass to THnSparse::Projection

    the slices/cuts use the actual values, not the bin numbers. they should correspond to bin edges though
    the axes numbers should be passed in "reverse" order; (y, x) or (z, y, x), for example
    """
    for (axisNumber, axisMin, axisMax) in slices:
        # https://root.cern.ch/doc/master/TAxis_8cxx_source.html#l00962
        # SetRangeUser appears to have a non-inclusive upper bound,
        # but just in case, we'll subtract epsilon because it doesn't hurt to do so
        # and if things ever change, it'll still behave the way we expect
        # to avoid precision issues, we use nextafter, which takes the number you want to change
        # and the direction you want to change it. since we always want to go down,
        # we set the direction to be negative infinity
        thnSparse.GetAxis(axisNumber).SetRangeUser(axisMin, np.nextafter(axisMax, -np.inf))
    projection = thnSparse.Projection(*axesToProject)

    # reset the axes so that we aren't carrying around these cuts
    for axis in thnSparse.GetListOfAxes():
        # apparently axis.UnZoom doesn't actually work, because why would it
        axis.SetRangeUser(-np.inf, np.inf)

    return projection


# mapping of ROOT TBranch type code to python type code for struct
rootBranchTypeToStructType = {}
rootBranchTypeToStructType['O'] = '?'  # boolean (boolean)
rootBranchTypeToStructType['B'] = 'b'  # signed char (integer) 8-bit
rootBranchTypeToStructType['b'] = 'B'  # unsigned char (integer) 8-bit
rootBranchTypeToStructType['S'] = 'h'  # signed short (integer) 16-bit
rootBranchTypeToStructType['s'] = 'H'  # unsigned short (integer) 16-bit
rootBranchTypeToStructType['I'] = 'i'  # signed integer (integer) 32-bit
rootBranchTypeToStructType['i'] = 'I'  # unsigned integer (integer) 32-bit
rootBranchTypeToStructType['L'] = 'l'  # signed long (integer) 64-bit
rootBranchTypeToStructType['l'] = 'L'  # unsigned long (integer) 64-bit
rootBranchTypeToStructType['F'] = 'f'  # float (float) 32-bit
rootBranchTypeToStructType['D'] = 'd'  # double (float) 64-bit


def tBranchToArray(branch, branchType, arrayShape):
    """
    Takes a multi-dimensional TBranch and returns a numpy array.
    Also useful for anything python doesn't automatically handle, like unsigned integers

    branch: TBranch, usually retrieved via getattr(tree, 'branchName')
    branchType: the ROOT corresponding to this TBranch
    arrayShape: tuple with the dimensions of the output; for a 1D array, it's just the length

    For TBranch type codes: https://root.cern.ch/doc/master/classTBranch.html#ac0412c423e6c8388b42247e0410cf822
    numpy is smart enough to handle a non-tuple arrayShape, so no need to check for that
    """
    totalLength = int(np.product(arrayShape))
    structType = rootBranchTypeToStructType[branchType]
    structFormat = '{0}{1}'.format(totalLength, structType)
    unpacked = struct.unpack(structFormat, branch)
    return np.reshape(unpacked, arrayShape)


def Chi2Histograms(df1, df2, var, bins):
    """
    Calculates chi2 between for a given variable between 2 dataframes

    df1: dataframe for comparison
    df2: dataframe for comparison
    var: name of variable to compare
    bins: bin edges for histograms
    """
    hist1, err1 = getNormHistAndErr(df1, var, bins)
    hist2, err2 = getNormHistAndErr(df2, var, bins)

    toterr = quadSumPairwise(err1, err2)

    numerator = np.square(np.subtract(hist1, hist2))
    denominator = np.square(toterr)

    chi2 = np.sum(np.divide(numerator, denominator))
    return chi2


def getTimeText(duration):
    """
    Converts a duration into the preferred human-readable format
    """
    if duration < 120:
        timetext = '{0:0.0f} seconds'.format(duration)
    elif duration < 3600:
        timetext = '{0:0.0f} minutes'.format(duration / 60.0)
    else:
        timetext = '{0:0.0f} hours {1:0.0f} minutes'.format(duration / 3600, (duration % 3600) / 60.0)
    return timetext


def getSizeText(size):
    """
    Converts a number of bytes into the preferred human-readable format
    """
    if size < 1024:
        sizetext = '{0} B'.format(size)
    elif size < 1024 * 1024:
        sizetext = '{0:0.1f} kB'.format(float(size) / 1024)
    elif size < 1024 * 1024 * 1024:
        sizetext = '{0:0.1f} MB'.format(float(size) / (1024 * 1024))
    else:
        sizetext = '{0:0.1f} GB'.format(float(size) / (1024 * 1024 * 1024))
    return sizetext


def getMixedEventFilename(config, mixlabel):
    """
    Takes the gamma-jet configuration and a mixing label and does some lookup and concatenation to get a filename for that particular mixed event ROOT file.
    """
    basename = config['filelists']['correlations']['mixedevent']
    return basename.replace('.root', '_{0}.root'.format(mixlabel))
