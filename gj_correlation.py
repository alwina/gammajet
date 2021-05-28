import matplotlib.pyplot as plt
import numpy as np
import ROOT
import root_numpy as rnp

from unfolder import Unfolder
from utils import getCenters, getWidths, getXerrForPlot, quadSumPairwise, th1ToArrays


# instructions for use
# 1. initialize with the observable name
# 2. call setSameEvent, setMixedEvent, and setResponseMatrix (optional)
# 3. call useNbins, which also turns the histograms into numpy arrays and divides by the bin widths
# 4. call subtractMixedEvent
# 5. call getSignalCorrelation
# 6. optionally, do stuff with the unfolder
class GammaJetCorrelation:
    def __init__(self, observable):
        self.observable = observable
        self.unfolder = Unfolder()

        self.nBins = -1
        self.minBin = -1
        self.maxBin = -1

    def binsMatch(self, th1):
        if self.nBins < 0:
            self.nBins = th1.GetNbinsX()
            self.minBin = th1.GetXaxis().GetXmin()
            self.maxBin = th1.GetXaxis().GetXmax()
            return True
        else:
            if self.nBins != th1.GetNbinsX():
                print('Number of bins does not match')
                return False
            if self.minBin != th1.GetXaxis().GetXmin():
                print('Bin minimum does not match')
                return False
            if self.maxBin != th1.GetXaxis().GetXmax():
                print('Bin maximum does not match')
                return False
            return True

    def setBinEdges(self):
        self.binEdges = np.linspace(self.minBin, self.maxBin, self.nBins + 1)
        self.centers = getCenters(self.binEdges)
        self.xerr = getXerrForPlot(self.binEdges)
        self.widths = getWidths(self.binEdges)

    def setSameEvent(self, srth1, brth1, ntrigsr, ntrigbr):
        if self.binsMatch(srth1):
            self.sesrth1 = srth1
        if self.binsMatch(brth1):
            self.sebrth1 = brth1

        # scale by number of triggers
        self.sesrth1.Scale(1.0 / ntrigsr)
        self.sebrth1.Scale(1.0 / ntrigbr)
        self.sesrntrig = ntrigsr
        self.sebrntrig = ntrigbr

    def setMixedEvent(self, srth1, brth1, ntrigsr, ntrigbr, srnmix, brnmix):
        if self.binsMatch(srth1):
            self.mesrth1 = srth1
        if self.binsMatch(brth1):
            self.mebrth1 = brth1

        # scale by number of triggers and number of events
        self.mesrth1.Scale(1.0 / (ntrigsr * srnmix))
        self.mebrth1.Scale(1.0 / (ntrigbr * brnmix))
        self.mesrntrig = ntrigsr
        self.mebrntrig = ntrigbr
        self.srnmix = srnmix
        self.brnmix = brnmix

    def setResponseMatrix(self, rooUnfoldResponse):
        self.unfolder.setResponseMatrix(rooUnfoldResponse)

    def convertAndDivideWidths(self, th1):
        hist, err, centers, widths = th1ToArrays(th1)
        return np.divide(hist, widths), np.divide(err, widths)

    def useNbins(self, nBins):
        rebin = self.nBins / nBins
        self.sesrth1.Rebin(rebin)
        self.sebrth1.Rebin(rebin)
        self.mesrth1.Rebin(rebin)
        self.mebrth1.Rebin(rebin)

        # set bin info
        self.nBins = nBins
        self.minBin = self.sesrth1.GetXaxis().GetXmin()
        self.maxBin = self.sesrth1.GetXaxis().GetXmax()
        self.binEdges = np.linspace(self.minBin, self.maxBin, self.nBins + 1)
        self.centers = getCenters(self.binEdges)
        self.xerr = getXerrForPlot(self.binEdges)
        self.widths = getWidths(self.binEdges)

        # turn into numpy arrays and divide by bin widths here
        self.sesrhist, self.sesrerr = self.convertAndDivideWidths(self.sesrth1)
        self.sebrhist, self.sebrerr = self.convertAndDivideWidths(self.sebrth1)
        self.mesrhist, self.mesrerr = self.convertAndDivideWidths(self.mesrth1)
        self.mebrhist, self.mebrerr = self.convertAndDivideWidths(self.mebrth1)

    def subtractMixedEvent(self):
        self.srhist = np.subtract(self.sesrhist, self.mesrhist)
        self.srerr = quadSumPairwise(self.sesrerr, self.mesrerr)
        self.brhist = np.subtract(self.sebrhist, self.mebrhist)
        self.brerr = quadSumPairwise(self.sebrerr, self.mebrerr)

    # rebin before calling this
    def getSignalCorrelation(self):
        self.corrhist = np.subtract(self.srhist, self.brhist)
        self.correrr = quadSumPairwise(self.srerr, self.brerr)

        # convert to TH1 and feed to unfolder
        self.corrth1 = ROOT.TH1F("{0}_measured".format(self.observable), "{0} measured".format(self.observable), self.nBins, self.minBin, self.maxBin)
        rnp.array2hist(self.corrhist, self.corrth1, self.correrr)
        self.unfolder.setMeasuredTH1(self.corrth1)

    def plotCorr(self, corrprefix, **kwargs):
        # retrieve the indicated histogram and error
        hist = getattr(self, '{0}hist'.format(corrprefix))
        err = getattr(self, '{0}err'.format(corrprefix))

        # set zero values to nan for plotting
        hist[hist == 0] = np.nan

        # jet pT is 1/pT dN/dpT
        if self.observable == 'jetpt':
            plt.errorbar(self.centers, np.divide(hist, np.multiply(self.widths, self.centers)),
                         yerr=err, xerr=self.xerr, **kwargs)
        else:
            plt.errorbar(self.centers, np.divide(hist, self.widths), yerr=err, xerr=self.xerr, **kwargs)
