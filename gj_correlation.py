from enum import Enum
import numpy as np
import ROOT

from utils import sliceAndProjectTHnSparse, getTH1MeanBinValue

# instructions for use
# 1. initialize with the observable name
# 2. call setSameEvent, setMixedEvent, and setResponseMatrix (optional)
# 3. call useNbins, which also turns the histograms into numpy arrays and divides by the bin widths
# 4. call subtractBkgRegion and/or subtractMixedEvent
# 5. call getSignalCorrelation
# 6. optionally, do stuff with the unfolder


class GammaJetCorrelation2D:
    """Manipulate ROOT objects for same- and mixed-event correlations to get
    2D histogram that can be unfolded

    2D response matrix axes: observable reco, observable truth, jetpt reco, jetpt truth
    """
    def __init__(self, observableX, observableY, observableInfo):
        self.observableX = observableX
        self.observableY = observableY
        self.observableInfo = observableInfo

        self.nBinsX = observableInfo[observableX]['nbins']
        self.minBinX = observableInfo[observableX]['min']
        self.maxBinX = observableInfo[observableX]['max']

        self.nBinsY = observableInfo[observableY]['nbins']
        self.minBinY = observableInfo[observableY]['min']
        self.maxBinY = observableInfo[observableY]['max']

        self.rebinFactorX = 120 / self.nBinsX
        self.rebinFactorY = 120 / self.nBinsY

    def setSameEvent(self, srth2, brth2, ntrigsr, ntrigbr):
        self.sesrth2 = srth2.Clone()
        self.sebrth2 = brth2.Clone()

        # scale by number of triggers
        if ntrigsr != 0:
            self.sesrth2.Scale(1.0 / ntrigsr)
        if ntrigbr != 0:
            self.sebrth2.Scale(1.0 / ntrigbr)
        self.sesrntrig = ntrigsr
        self.sebrntrig = ntrigbr

        self.sesrth2finebins = self.sesrth2.Clone()
        self.sebrth2finebins = self.sebrth2.Clone()

        # rebin
        self.sesrth2.Rebin2D(self.rebinFactorX, self.rebinFactorY)
        self.sebrth2.Rebin2D(self.rebinFactorX, self.rebinFactorY)

        # subtract
        self.seth2 = self.sesrth2.Clone()
        self.seth2.Add(self.sebrth2, -1.0)
        self.seth2finebins = self.sesrth2finebins.Clone()
        self.seth2finebins.Add(self.sebrth2finebins, -1.0)

    def setMixedEvent(self, srth2, brth2, srnmix, brnmix):
        self.srnmix = srnmix
        self.brnmix = brnmix

        self.mesrth2 = srth2.Clone()
        self.mebrth2 = brth2.Clone()

        self.mesrth2finebins = self.mesrth2.Clone()
        self.mebrth2finebins = self.mebrth2.Clone()

        # rebin
        self.mesrth2.Rebin2D(self.rebinFactorX, self.rebinFactorY)
        self.mebrth2.Rebin2D(self.rebinFactorX, self.rebinFactorY)

        # subtract
        self.meth2 = self.mesrth2.Clone()
        self.meth2.Add(self.mebrth2, -1.0)
        self.meth2finebins = self.mesrth2finebins.Clone()
        self.meth2finebins.Add(self.mebrth2finebins, -1.0)

    def doCorrelationSubtraction(self):
        self.corrth2 = self.sesrth2.Clone()
        self.corrth2.Add(self.sebrth2, -1.0)
        self.corrth2.Add(self.mesrth2, -1.0)
        self.corrth2.Add(self.mebrth2, 1.0)

        self.corrth2finebins = self.sesrth2finebins.Clone()
        self.corrth2finebins.Add(self.sebrth2finebins, -1.0)
        self.corrth2finebins.Add(self.mesrth2finebins, -1.0)
        self.corrth2finebins.Add(self.mebrth2finebins, 1.0)

        self.totalbkgth2 = self.sebrth2.Clone()
        self.totalbkgth2.Add(self.mesrth2, 1.0)
        self.totalbkgth2.Add(self.mebrth2, -1.0)

        self.totalbkgth2finebins = self.sebrth2finebins.Clone()
        self.totalbkgth2finebins.Add(self.mesrth2finebins, 1.0)
        self.totalbkgth2finebins.Add(self.mebrth2finebins, -1.0)

    def getSignalCorrelation(self, specialBins=[]):
        self.unfolder.setMeasuredTH2(self.corrth2)

        # make sure we always divide by bin widths! use this instead of corrth2.ProjectionX()!
        self.corrth1regbins = self.corrth2.ProjectionX()
        self.corrth1finebins = self.corrth2finebins.ProjectionX()

        if specialBins:
            self.corrth1 = self.corrth1regbins.Rebin(5, '', np.array(specialBins))
        else:
            self.corrth1 = self.corrth1regbins.Clone()

        self.corrth1.Scale(1.0, 'width')
        self.corrth1regbins.Scale(1.0, 'width')
        self.corrth1finebins.Scale(1.0, 'width')


class AxisNum(Enum):
    """Enum describing what each THnSparse axis represents"""
    centrality = 0
    clusterpt = 1
    deltaphi = 2
    jetpt = 3
    ptratio = 4
    jetkt = 5


def getAll2DCorr(centranges, photonptranges, observables, observableInfo, rootfileSE, rootfileME, allMEScales, **kwargs):
    # set up the dictionary for all correlation objects
    allCorr = {}

    # parsing ROOT files - extracting THnSparses
    rootfile = ROOT.TFile.Open(rootfileSE)
    sehTrigSR = rootfile.Get('hTrigSR')
    sehTrigBR = rootfile.Get('hTrigBR')
    sehCorrSR = rootfile.Get('hCorrSR')
    sehCorrBR = rootfile.Get('hCorrBR')
    rootfile.Close()

    rootfile = ROOT.TFile.Open(rootfileME)
    mehnMixSR = rootfile.Get('hnMixSR')
    mehnMixBR = rootfile.Get('hnMixBR')
    mehCorrSR = rootfile.Get('hCorrSR')
    mehCorrBR = rootfile.Get('hCorrBR')
    rootfile.Close()

    for i, centrange in enumerate(centranges):
        allCorr[centrange] = {}
        for j, ptrange in enumerate(photonptranges):
            allCorr[centrange][ptrange] = {}

            slices = []
            slices.append((AxisNum.centrality.value, centrange[0], centrange[1]))
            slices.append((AxisNum.clusterpt.value, ptrange[0], ptrange[1]))

            # get the number of triggers for each bin in the same-event
            nTrigSESR = sliceAndProjectTHnSparse(sehTrigSR, slices, 0).Integral()
            nTrigSEBR = sliceAndProjectTHnSparse(sehTrigBR, slices, 0).Integral()

            # get the number of trigger-MB event pairs for each bin in the mixed-event
            nMixSR = sliceAndProjectTHnSparse(mehnMixSR, slices, 0).Integral()
            nMixBR = sliceAndProjectTHnSparse(mehnMixBR, slices, 0).Integral()

            # get the scaling for ME
            mescales = allMEScales[centrange][ptrange]
            mesrscale = mescales['sr']
            mebrscale = mescales['br']

            for observable in observables:
                gjCorr = GammaJetCorrelation2D(observable, 'jetpt', observableInfo)

                additionalCuts = []
                for cutvar in observableInfo[observable]['cuts']:
                    cut = observableInfo[observable]['cuts'][cutvar]
                    additionalCuts.append((AxisNum[cutvar].value, cut['min'], cut['max']))

                srth2 = sliceAndProjectTHnSparse(sehCorrSR, slices + additionalCuts, AxisNum['jetpt'].value, AxisNum[observable].value)
                brth2 = sliceAndProjectTHnSparse(sehCorrBR, slices + additionalCuts, AxisNum['jetpt'].value, AxisNum[observable].value)
                gjCorr.setSameEvent(srth2, brth2, nTrigSESR, nTrigSEBR)

                srth2 = combineME(centrange, ptrange, ['jetpt', observable], nTrigSESR, sehTrigSR, mehnMixSR, mehCorrSR, additionalCuts)
                brth2 = combineME(centrange, ptrange, ['jetpt', observable], nTrigSEBR, sehTrigBR, mehnMixBR, mehCorrBR, additionalCuts)
                srth2.Scale(mesrscale)
                brth2.Scale(mebrscale)
                gjCorr.setMixedEvent(srth2, brth2, nMixSR, nMixBR)

                gjCorr.doCorrelationSubtraction()
                specialBins = []
                if observable == 'deltaphi':
                    if centrange in [(0, 10), (0, 30)]:
                        # combine first two bins in certain cases
                        tempth1 = gjCorr.corrth2.ProjectionX()
                        specialBins = [tempth1.GetBinLowEdge(i) for i in [1, 3, 4, 5, 6, 7]]
                gjCorr.getSignalCorrelation(specialBins)

                allCorr[centrange][ptrange][observable] = gjCorr

    return allCorr


def getAllMEScales(centranges, photonptranges, observableInfo, rootfileSE, rootfileME, jetptranges, deltaphirange):
    """
    get scales for MESR and MEBR based on a negative jet pT range and the ratio SESR/MESR or SEBR/MEBR in some deltaphi range
    """
    allMEScales = {}

    # parsing ROOT files - extracting THnSparses
    rootfile = ROOT.TFile.Open(rootfileSE)
    sehTrigSR = rootfile.Get('hTrigSR')
    sehTrigBR = rootfile.Get('hTrigBR')
    sehCorrSR = rootfile.Get('hCorrSR')
    sehCorrBR = rootfile.Get('hCorrBR')
    rootfile.Close()

    rootfile = ROOT.TFile.Open(rootfileME)
    mehnMixSR = rootfile.Get('hnMixSR')
    mehnMixBR = rootfile.Get('hnMixBR')
    mehCorrSR = rootfile.Get('hCorrSR')
    mehCorrBR = rootfile.Get('hCorrBR')
    rootfile.Close()

    for i, centrange in enumerate(centranges):
        allMEScales[centrange] = {}
        jetptrange = jetptranges[centrange]
        for j, ptrange in enumerate(photonptranges):
            allMEScales[centrange][ptrange] = {}

            slices = []
            slices.append((AxisNum.centrality.value, centrange[0], centrange[1]))
            slices.append((AxisNum.clusterpt.value, ptrange[0], ptrange[1]))

            # get the number of triggers for each bin in the same-event
            nTrigSESR = sliceAndProjectTHnSparse(sehTrigSR, slices, 0).Integral()
            nTrigSEBR = sliceAndProjectTHnSparse(sehTrigBR, slices, 0).Integral()

            jetptslice = (AxisNum.jetpt.value, jetptrange[0], jetptrange[1])
            rebin = 120 / observableInfo['deltaphi']['nbins']

            # get per-trigger yields for SE
            sesrh1 = sliceAndProjectTHnSparse(sehCorrSR, slices + [jetptslice], AxisNum.deltaphi.value)
            sebrh1 = sliceAndProjectTHnSparse(sehCorrBR, slices + [jetptslice], AxisNum.deltaphi.value)
            sesrh1.Scale(1.0 / nTrigSESR)
            sebrh1.Scale(1.0 / nTrigSEBR)

            # get per-trigger yields for ME
            mesrh1 = combineME(centrange, ptrange, ['deltaphi'], nTrigSESR, sehTrigSR, mehnMixSR, mehCorrSR, [jetptslice])
            mebrh1 = combineME(centrange, ptrange, ['deltaphi'], nTrigSEBR, sehTrigBR, mehnMixBR, mehCorrBR, [jetptslice])

            # rebin
            sesrh1.Rebin(rebin)
            sebrh1.Rebin(rebin)
            mesrh1.Rebin(rebin)
            mebrh1.Rebin(rebin)

            # get ratios
            ratiosr = sesrh1.Clone()
            ratiosr.Divide(mesrh1)
            ratiobr = sebrh1.Clone()
            ratiobr.Divide(mebrh1)

            # restrict to deltaphi range
            ratiosr.GetXaxis().SetRangeUser(deltaphirange[0], deltaphirange[1])
            ratiobr.GetXaxis().SetRangeUser(deltaphirange[0], deltaphirange[1])

            allMEScales[centrange][ptrange]['sr'] = {}
            avgscale = getTH1MeanBinValue(ratiosr, deltaphirange[0], deltaphirange[1])
            minscale = ratiosr.GetBinContent(ratiosr.GetMinimumBin())
            maxscale = ratiosr.GetBinContent(ratiosr.GetMaximumBin())
            allMEScales[centrange][ptrange]['sr']['avg'] = avgscale
            allMEScales[centrange][ptrange]['sr']['min'] = minscale
            allMEScales[centrange][ptrange]['sr']['max'] = maxscale

            allMEScales[centrange][ptrange]['br'] = {}
            avgscale = getTH1MeanBinValue(ratiobr, deltaphirange[0], deltaphirange[1])
            minscale = ratiobr.GetBinContent(ratiobr.GetMinimumBin())
            maxscale = ratiobr.GetBinContent(ratiobr.GetMaximumBin())
            allMEScales[centrange][ptrange]['br']['avg'] = avgscale
            allMEScales[centrange][ptrange]['br']['min'] = minscale
            allMEScales[centrange][ptrange]['br']['max'] = maxscale

    return allMEScales


def combineME(centrange, ptrange, observables, nTrig, hTrig, hnMix, hCorr, additionalCuts=[]):
    """
    combine ME per-trigger yields such that the centrality distributions match SE in 10% bins
    """
    smallmehists = []
    ntrigs = []

    for centmin in range(centrange[0], centrange[1], 10):
        otherslices = []
        otherslices.append((AxisNum.centrality.value, centmin, centmin + 10))
        otherslices.append((AxisNum.clusterpt.value, ptrange[0], ptrange[1]))
        nTrigsmall = sliceAndProjectTHnSparse(hTrig, otherslices, 0).Integral()
        nMixsmall = sliceAndProjectTHnSparse(hnMix, otherslices, 0).Integral()
        axesToProject = [AxisNum[observable].value for observable in observables]
        mesmall = sliceAndProjectTHnSparse(hCorr, otherslices + additionalCuts, *axesToProject)
        mesmall.Scale(1.0 / nMixsmall)

        smallmehists.append(mesmall)
        ntrigs.append(nTrigsmall)

    mehist = smallmehists[0].Clone()
    mehist.Reset()
    for h1, trig, in zip(smallmehists, ntrigs):
        mehist.Add(h1, float(trig) / nTrig)

    return mehist
