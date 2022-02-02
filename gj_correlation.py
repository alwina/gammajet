from enum import Enum
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import root_numpy as rnp

from unfolder import Unfolder
from utils import getCenters, getWidths, getXerrForPlot, quadSumPairwise, th1ToArrays, sliceAndProjectTHnSparse

# instructions for use
# 1. initialize with the observable name
# 2. call setSameEvent, setMixedEvent, and setResponseMatrix (optional)
# 3. call useNbins, which also turns the histograms into numpy arrays and divides by the bin widths
# 4. call subtractBkgRegion and/or subtractMixedEvent
# 5. call getSignalCorrelation
# 6. optionally, do stuff with the unfolder


class GammaJetCorrelation:
    """Manipulate ROOT objects for same- and mixed-event correlations to get various
    stages of subtracted correlations and generate plots"""
    def __init__(self, observable):
        self.observable = observable
        self.unfolder = Unfolder()

        self.nBins = -1
        self.minBin = np.inf
        self.maxBin = -np.inf

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

    def setSameEvent(self, srth1, brth1, ntrigsr, ntrigbr):
        if self.binsMatch(srth1):
            self.sesrth1 = srth1
        if self.binsMatch(brth1):
            self.sebrth1 = brth1

        # scale by number of triggers
        if ntrigsr != 0:
            self.sesrth1.Scale(1.0 / ntrigsr)
        if ntrigbr != 0:
            self.sebrth1.Scale(1.0 / ntrigbr)
        self.sesrntrig = ntrigsr
        self.sebrntrig = ntrigbr

    def setMixedEvent(self, srth1, brth1, srnmix, brnmix, mesrscale=1.0, mebrscale=1.0):
        if self.binsMatch(srth1):
            self.mesrth1 = srth1
        if self.binsMatch(brth1):
            self.mebrth1 = brth1

        # scale by number of trigger-event pairs and any other scaling
        if srnmix != 0:
            self.mesrth1.Scale(mesrscale / srnmix)
        if brnmix != 0:
            self.mebrth1.Scale(mebrscale / brnmix)
        self.srnmix = srnmix
        self.brnmix = brnmix

    def setResponseMatrix(self, rooUnfoldResponse):
        self.unfolder.setResponseMatrix(rooUnfoldResponse)

    def convertAndDivideWidths(self, th1):
        hist, err, centers, widths = th1ToArrays(th1)

        if not np.allclose(centers, self.centers):
            print('Warning: TH1 centers {0} do not match calculated centers {1}'.format(centers, self.centers))
        if not np.allclose(widths, self.widths):
            print('Warning: TH1 widths {0} do not match calculated widths {1}'.format(widths, self.widths))

        return np.divide(hist, widths), np.divide(err, widths)

    def useNbins(self, nBins):
        # calculate the mean of the distribution before rebinning
        # and also just keep around the TH1 with all the bins
        self.signalth1 = self.sesrth1.Clone()
        self.signalth1.Add(self.sebrth1, -1.0)
        self.signalth1.Add(self.mesrth1, -1.0)
        self.signalth1.Add(self.mebrth1, 1.0)
        self.mean = self.signalth1.GetMean()
        self.meanerr = self.signalth1.GetMeanError()
        self.rms = self.signalth1.GetRMS()
        self.rmserr = self.signalth1.GetRMSError()

        rebinFactor = self.nBins / nBins
        self.sesrth1.Rebin(rebinFactor)
        self.sebrth1.Rebin(rebinFactor)
        self.mesrth1.Rebin(rebinFactor)
        self.mebrth1.Rebin(rebinFactor)

        # set bin info
        self.nBins = nBins
        self.binEdges = np.linspace(self.minBin, self.maxBin, self.nBins + 1)
        self.centers = getCenters(self.binEdges)
        self.xerr = getXerrForPlot(self.binEdges)
        self.widths = getWidths(self.binEdges)

        # check for issues with rebinning
        if not (self.binsMatch(self.sesrth1) and self.binsMatch(self.sebrth1) and self.binsMatch(self.mesrth1) and self.binsMatch(self.mebrth1)):
            print('Warning: issue with rebinning')

        # turn into numpy arrays and divide by bin widths here
        self.sesrhist, self.sesrerr = self.convertAndDivideWidths(self.sesrth1)
        self.sebrhist, self.sebrerr = self.convertAndDivideWidths(self.sebrth1)
        self.mesrhist, self.mesrerr = self.convertAndDivideWidths(self.mesrth1)
        self.mebrhist, self.mebrerr = self.convertAndDivideWidths(self.mebrth1)

    def getSignalCorrelation(self):
        # it absolutely should not matter if we do (SESR-SEBR)-(MESR-MEBR)
        # or (SESR-MESR)-(SEBR-MEBR); just do SE - ME
        # totalbkg = MESR + SEBR - MEBR = SEBR + MESR - MEBR
        try:
            self.sehist = np.subtract(self.sesrhist, self.sebrhist)
            self.seerr = quadSumPairwise(self.sesrerr, self.sebrerr)
            self.mehist = np.subtract(self.mesrhist, self.mebrhist)
            self.meerr = quadSumPairwise(self.mesrerr, self.mebrerr)

            self.srhist = np.subtract(self.sesrhist, self.mesrhist)
            self.srerr = quadSumPairwise(self.sesrerr, self.mesrerr)
            self.brhist = np.subtract(self.sebrhist, self.mebrhist)
            self.brerr = quadSumPairwise(self.sebrerr, self.mebrerr)
        except AttributeError:
            print('Must call useNbins before calling getSignalCorrelation')
            raise

        self.corrhist = np.subtract(self.sehist, self.mehist)
        self.correrr = quadSumPairwise(self.seerr, self.meerr)
        self.totalbkghist = np.add(self.sebrhist, self.mehist)
        self.totalbkgerr = quadSumPairwise(self.sebrerr, self.meerr)

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

        # plot
        plt.errorbar(self.centers, hist, yerr=err, xerr=self.xerr, **kwargs)


class AxisNum(Enum):
    """Enum describing what each THnSparse axis represents"""
    centrality = 0
    clusterpt = 1
    deltaphi = 2
    jetpt = 3
    ptratio = 4
    jetarea = 5
    jetmultiplicity = 6


def getAllCorr(centranges, photonptranges, observableInfo, rootfileSE, rootfileME, rootfileRM='', **kwargs):
    """
    Parse ROOT files and generate GammaJetCorrelation objects for each centrality range, ptrange, and observable

    centranges: list of tuples of centrality ranges
    photonptranges: list of tuples of photon pT ranges
    observableInfo: dictionary of observable info {name: {info}}, generally read from config yaml
    rootfileSE: name of ROOT file for same-event correlation
    rootfileME: name of ROOT file for mixed-event correlation
    rootfileRM: name of ROOT file for RooUnfoldResponse object (optional)
    """
    # set up the dictionary for all correlation objects
    allCorr = {}

    # parsing ROOT files - extracting THnSparses
    rootfile = ROOT.TFile.Open(rootfileSE)
    sehTrigSR = rootfile.Get(kwargs.get('sehTrigSR', 'hTrigSR'))
    sehTrigBR = rootfile.Get(kwargs.get('sehTrigBR', 'hTrigBR'))
    sehCorrSR = rootfile.Get(kwargs.get('sehCorrSR', 'hCorrSR'))
    sehCorrBR = rootfile.Get(kwargs.get('sehCorrBR', 'hCorrBR'))
    sehCorr1ptSR = rootfile.Get(kwargs.get('sehCorr1ptSR', 'hCorr1ptSR'))
    sehCorr1ptBR = rootfile.Get(kwargs.get('sehCorr1ptBR', 'hCorr1ptBR'))
    rootfile.Close()

    rootfile = ROOT.TFile.Open(rootfileME)
    mehnMixSR = rootfile.Get(kwargs.get('mehnMixSR', 'hnMixSR'))
    mehnMixBR = rootfile.Get(kwargs.get('mehnMixBR', 'hnMixBR'))
    mehCorrSR = rootfile.Get(kwargs.get('mehCorrSR', 'hCorrSR'))
    mehCorrBR = rootfile.Get(kwargs.get('mehCorrBR', 'hCorrBR'))
    mehCorr1ptSR = rootfile.Get(kwargs.get('mehCorr1ptSR', 'hCorr1ptSR'))
    mehCorr1ptBR = rootfile.Get(kwargs.get('mehCorr1ptBR', 'hCorr1ptBR'))
    rootfile.Close()

    if rootfileRM:
        rootfile = ROOT.TFile.Open(rootfileRM)
        # TO-DO: split by pT
        rooUnfoldResponses = {}
        for observable in observableInfo:
            rooUnfoldResponses[observable] = {}
            for icent, centrange in enumerate(centranges):
                rooUnfoldResponses[observable][centrange] = rootfile.Get('{0}Response{1}'.format(observable, icent))
        rootfile.Close()

    # make all correlation objects
    for centrange in centranges:
        allCorr[centrange] = {}
        for photonptrange in photonptranges:
            allCorr[centrange][photonptrange] = {}
            slices = []
            slices.append((AxisNum.centrality.value, centrange[0], centrange[1]))
            slices.append((AxisNum.clusterpt.value, photonptrange[0], photonptrange[1]))

            # get the number of triggers for each bin in the same-event
            hTrigSESR = sliceAndProjectTHnSparse(sehTrigSR, slices, 0)
            nTrigSESR = hTrigSESR.Integral()

            hTrigSEBR = sliceAndProjectTHnSparse(sehTrigBR, slices, 0)
            nTrigSEBR = hTrigSEBR.Integral()

            # get the number of trigger-MB event pairs for each bin in the mixed-event
            hnMixSR = sliceAndProjectTHnSparse(mehnMixSR, slices, 0)
            nMixSR = hnMixSR.Integral()

            hnMixBR = sliceAndProjectTHnSparse(mehnMixBR, slices, 0)
            nMixBR = hnMixBR.Integral()

            # look for more mixed-event scaling
            try:
                mesrscale = kwargs['mescaling'][centrange][photonptrange]['sr']
                mebrscale = kwargs['mescaling'][centrange][photonptrange]['br']
            except KeyError:
                mesrscale = 1.0
                mebrscale = 1.0

            for observable in observableInfo:
                # set up the correlation object
                gjCorr = GammaJetCorrelation(observable)

                # add any additional cuts for the particular observable
                additionalCuts = []
                for cutvar in observableInfo[observable]['cuts']:
                    cut = observableInfo[observable]['cuts'][cutvar]
                    additionalCuts.append((AxisNum[cutvar].value, cut['min'], cut['max']))

                # project and scale the same-event
                if observable == 'jetpt':
                    try:
                        srTH1 = sliceAndProjectTHnSparse(sehCorr1ptSR, slices + additionalCuts, AxisNum[observable].value)
                        brTH1 = sliceAndProjectTHnSparse(sehCorr1ptBR, slices + additionalCuts, AxisNum[observable].value)
                    except AttributeError:
                        srTH1 = sliceAndProjectTHnSparse(sehCorrSR, slices + additionalCuts, AxisNum[observable].value)
                        brTH1 = sliceAndProjectTHnSparse(sehCorrBR, slices + additionalCuts, AxisNum[observable].value)
                        print('Warning: SE jetpt distribution does not have 1/pT factor')
                else:
                    srTH1 = sliceAndProjectTHnSparse(sehCorrSR, slices + additionalCuts, AxisNum[observable].value)
                    brTH1 = sliceAndProjectTHnSparse(sehCorrBR, slices + additionalCuts, AxisNum[observable].value)

                gjCorr.setSameEvent(srTH1, brTH1, nTrigSESR, nTrigSEBR)

                # project and scale the mixed-event
                if observable == 'jetpt':
                    try:
                        srTH1 = sliceAndProjectTHnSparse(mehCorr1ptSR, slices + additionalCuts, AxisNum[observable].value)
                        brTH1 = sliceAndProjectTHnSparse(mehCorr1ptBR, slices + additionalCuts, AxisNum[observable].value)
                    except AttributeError:
                        srTH1 = sliceAndProjectTHnSparse(mehCorrSR, slices + additionalCuts, AxisNum[observable].value)
                        brTH1 = sliceAndProjectTHnSparse(mehCorrBR, slices + additionalCuts, AxisNum[observable].value)
                        print('Warning: ME jetpt distribution does not have 1/pT factor')
                else:
                    srTH1 = sliceAndProjectTHnSparse(mehCorrSR, slices + additionalCuts, AxisNum[observable].value)
                    brTH1 = sliceAndProjectTHnSparse(mehCorrBR, slices + additionalCuts, AxisNum[observable].value)

                gjCorr.setMixedEvent(srTH1, brTH1, nMixSR, nMixBR, mesrscale, mebrscale)

                # rebin
                gjCorr.useNbins(observableInfo[observable]['nbins'])

                # subtract backgrounds, which includes both subtracting BR in each of SE and ME
                # and also subtracting ME in each of SR and BR
                gjCorr.getSignalCorrelation()

                # set response matrix if it exists
                # TO-DO: split by pT
                gjCorr.setResponseMatrix(rooUnfoldResponses[observable][centrange])

                # add to dictionary
                allCorr[centrange][photonptrange][observable] = gjCorr

    return allCorr
