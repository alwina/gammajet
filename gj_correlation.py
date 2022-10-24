from enum import Enum
import ROOT

from unfolder import Unfolder2D
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
    def __init__(self, observableX, observableY, observableInfo, unfoldVerbosity=1):
        self.observableX = observableX
        self.observableY = observableY
        self.observableInfo = observableInfo
        self.unfolder = Unfolder2D(unfoldVerbosity)

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

    def setMixedEvent(self, srth2, brth2, srnmix, brnmix, mesrscale=1.0, mebrscale=1.0):
        self.mesrth2 = srth2.Clone()
        self.mebrth2 = brth2.Clone()

        # scale by number of trigger-event pairs and any other scaling
        if srnmix != 0:
            self.mesrth2.Scale(mesrscale / srnmix)
        if brnmix != 0:
            self.mebrth2.Scale(mebrscale / brnmix)
        self.srnmix = srnmix
        self.brnmix = brnmix

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

    def getSignalCorrelation(self):
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

        self.unfolder.setMeasuredTH2(self.corrth2)

        # make sure we always divide by bin widths! use this instead of corrth2.ProjectionX()!
        self.corrth1 = self.corrth2.ProjectionX()
        self.corrth1.Scale(1.0, 'width')

        self.corrth1finebins = self.corrth2finebins.ProjectionX()
        self.corrth1finebins.Scale(1.0, 'width')


class GJMCCorrelation2D:
    def __init__(self, centrange, ptrange, observable, observableInfo, minjetrecopt=10, minjettruthpt=0):
        self.centrange = centrange
        self.ptrange = ptrange
        self.observable = observable
        self.observableInfo = observableInfo
        self.minjetrecopt = minjetrecopt
        self.minjettruthpt = minjettruthpt

        self.nBinsX = observableInfo[observable]['nbins']
        self.minBinX = observableInfo[observable]['min']
        self.maxBinX = observableInfo[observable]['max']

        self.nBinsY = observableInfo['jetpt']['nbins']
        self.minBinY = observableInfo['jetpt']['min']
        self.maxBinY = observableInfo['jetpt']['max']

        self.rebinFactorX = 120 / self.nBinsX
        self.rebinFactorY = 120 / self.nBinsY

    def getCorrelations(self, hTrigSR, hCorrSRTruth, hCorrSRAll, th4):
        # these jets are not required to be truth-reco matched
        slices = []
        slices.append((AxisNum.centrality.value, self.centrange[0], self.centrange[1]))
        slices.append((AxisNum.clusterpt.value, self.ptrange[0], self.ptrange[1]))

        self.nTrig = sliceAndProjectTHnSparse(hTrigSR, slices, 0).Integral()

        for cutvar in self.observableInfo[self.observable]['cuts']:
            # we handle this one separately
            if cutvar == 'jetpt':
                continue
            cut = self.observableInfo[self.observable]['cuts'][cutvar]
            slices.append((AxisNum[cutvar].value, cut['min'], cut['max']))

        self.unmatchedrecoth2 = sliceAndProjectTHnSparse(hCorrSRAll, slices + [(AxisNum.jetpt.value, self.minjetrecopt, self.maxBinY)], AxisNum.jetpt.value, AxisNum[self.observable].value)
        self.unmatchedtruthth2 = sliceAndProjectTHnSparse(hCorrSRTruth, slices + [(AxisNum.jetpt.value, self.minjettruthpt, self.maxBinY)], AxisNum.jetpt.value, AxisNum[self.observable].value)

        # TH4 axes: reco obs, truth obs, reco jetpt, truth jetpt
        # in principle, this is the same as the response matrix
        # these are already split by centrality and pT
        th4slices = []
        th4slices.append((2, self.minjetrecopt, self.maxBinY))
        th4slices.append((3, self.minjettruthpt, self.maxBinY))

        self.th4 = th4
        self.matchedrecoth2 = sliceAndProjectTHnSparse(th4, th4slices, 2, 0)
        self.matchedtruthth2 = sliceAndProjectTHnSparse(th4, th4slices, 3, 1)

        # scale by number of triggers and rebin
        self.unmatchedrecoth2.Scale(1.0 / self.nTrig)
        self.unmatchedtruthth2.Scale(1.0 / self.nTrig)
        self.matchedrecoth2.Scale(1.0 / self.nTrig)
        self.matchedtruthth2.Scale(1.0 / self.nTrig)

        self.unmatchedrecoth2.Rebin2D(self.rebinFactorX, self.rebinFactorY)
        self.unmatchedtruthth2.Rebin2D(self.rebinFactorX, self.rebinFactorY)
        self.matchedrecoth2.Rebin2D(self.rebinFactorX, self.rebinFactorY)
        self.matchedtruthth2.Rebin2D(self.rebinFactorX, self.rebinFactorY)


class AxisNum(Enum):
    """Enum describing what each THnSparse axis represents"""
    centrality = 0
    clusterpt = 1
    deltaphi = 2
    jetpt = 3
    ptratio = 4
    jetkt = 5


def getAll2DCorr(centranges, photonptranges, observables, observableInfo, rootfileSE, rootfileME, rmfilename, allMEScales, unfoldVerbosity=1, **kwargs):
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

    rmfile = ROOT.TFile.Open(rmfilename)

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
                gjCorr = GammaJetCorrelation2D(observable, 'jetpt', observableInfo, unfoldVerbosity)

                additionalCuts = []
                for cutvar in observableInfo[observable]['cuts']:
                    cut = observableInfo[observable]['cuts'][cutvar]
                    additionalCuts.append((AxisNum[cutvar].value, cut['min'], cut['max']))

                srth2 = sliceAndProjectTHnSparse(sehCorrSR, slices + additionalCuts, AxisNum['jetpt'].value, AxisNum[observable].value)
                brth2 = sliceAndProjectTHnSparse(sehCorrBR, slices + additionalCuts, AxisNum['jetpt'].value, AxisNum[observable].value)
                gjCorr.setSameEvent(srth2, brth2, nTrigSESR, nTrigSEBR)

                if centrange == (0, 30):
                    # separate 0-10% and 10-30, calculate ME correlation
                    # then recombine with weights from SE centrality distribution

                    # calculate per-trigger yield for 0-10%
                    otherslices = []
                    otherslices.append((AxisNum.centrality.value, 0, 10))
                    otherslices.append((AxisNum.clusterpt.value, ptrange[0], ptrange[1]))
                    nTrigsr010 = sliceAndProjectTHnSparse(sehTrigSR, otherslices, 0).Integral()
                    nTrigbr010 = sliceAndProjectTHnSparse(sehTrigBR, otherslices, 0).Integral()
                    nMixsr010 = sliceAndProjectTHnSparse(mehnMixSR, otherslices, 0).Integral()
                    nMixbr010 = sliceAndProjectTHnSparse(mehnMixBR, otherslices, 0).Integral()
                    mesr010 = sliceAndProjectTHnSparse(mehCorrSR, otherslices + additionalCuts, AxisNum.jetpt.value, AxisNum[observable].value)
                    mebr010 = sliceAndProjectTHnSparse(mehCorrBR, otherslices + additionalCuts, AxisNum.jetpt.value, AxisNum[observable].value)
                    mesr010.Scale(mesrscale / nMixsr010)
                    mebr010.Scale(mebrscale / nMixbr010)

                    # calculate per-trigger yield for 10-30%
                    otherslices = []
                    otherslices.append((AxisNum.centrality.value, 10, 30))
                    otherslices.append((AxisNum.clusterpt.value, ptrange[0], ptrange[1]))
                    nTrigsr1030 = sliceAndProjectTHnSparse(sehTrigSR, otherslices, 0).Integral()
                    nTrigbr1030 = sliceAndProjectTHnSparse(sehTrigBR, otherslices, 0).Integral()
                    nMixsr1030 = sliceAndProjectTHnSparse(mehnMixSR, otherslices, 0).Integral()
                    nMixbr1030 = sliceAndProjectTHnSparse(mehnMixBR, otherslices, 0).Integral()
                    mesr1030 = sliceAndProjectTHnSparse(mehCorrSR, otherslices + additionalCuts, AxisNum.jetpt.value, AxisNum[observable].value)
                    mebr1030 = sliceAndProjectTHnSparse(mehCorrBR, otherslices + additionalCuts, AxisNum.jetpt.value, AxisNum[observable].value)
                    mesr1030.Scale(mesrscale / nMixsr1030)
                    mebr1030.Scale(mebrscale / nMixbr1030)

                    # combine 0-10% and 10-30%
                    srth2 = mesr010.Clone()
                    srth2.Scale(float(nTrigsr010) / nTrigSESR)
                    srth2.Add(mesr1030, float(nTrigsr1030) / nTrigSESR)
                    brth2 = mebr010.Clone()
                    brth2.Scale(float(nTrigbr010) / nTrigSEBR)
                    brth2.Add(mebr1030, float(nTrigbr1030) / nTrigSEBR)

                    # modified version of setMixedEvent
                    gjCorr.mesrth2finebins = srth2.Clone()
                    gjCorr.mebrth2finebins = brth2.Clone()
                    gjCorr.srnmix = nMixSR
                    gjCorr.brnmix = nMixBR

                    srth2.Rebin2D(gjCorr.rebinFactorX, gjCorr.rebinFactorY)
                    brth2.Rebin2D(gjCorr.rebinFactorX, gjCorr.rebinFactorY)
                    gjCorr.mesrth2 = srth2
                    gjCorr.mebrth2 = brth2

                    gjCorr.meth2 = gjCorr.mesrth2.Clone()
                    gjCorr.meth2.Add(gjCorr.mebrth2, -1.0)
                    gjCorr.meth2finebins = gjCorr.mesrth2finebins.Clone()
                    gjCorr.meth2finebins.Add(gjCorr.mebrth2finebins, -1.0)
                else:
                    srth2 = sliceAndProjectTHnSparse(mehCorrSR, slices + additionalCuts, AxisNum['jetpt'].value, AxisNum[observable].value)
                    brth2 = sliceAndProjectTHnSparse(mehCorrBR, slices + additionalCuts, AxisNum['jetpt'].value, AxisNum[observable].value)
                    gjCorr.setMixedEvent(srth2, brth2, nMixSR, nMixBR, mesrscale, mebrscale)

                gjCorr.unfolder.setResponseMatrix(rmfile.Get('{0}{1}Response{2}{3}'.format(observable, 'jetpt', i, j)))
                gjCorr.unfolder.setTH4(rmfile.Get('{0}{1}Hist{2}{3}'.format(observable, 'jetpt', i, j)))

                gjCorr.getSignalCorrelation()

                allCorr[centrange][ptrange][observable] = gjCorr

    rmfile.Close()
    return allCorr


def getAllMEScales(centranges, photonptranges, observables, observableInfo, rootfileSE, rootfileME, jetptranges, deltaphirange):
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

            # get the number of trigger-MB event pairs for each bin in the mixed-event
            nMixSR = sliceAndProjectTHnSparse(mehnMixSR, slices, 0).Integral()
            nMixBR = sliceAndProjectTHnSparse(mehnMixBR, slices, 0).Integral()

            jetptslice = (AxisNum.jetpt.value, jetptrange[0], jetptrange[1])
            rebin = 120 / observableInfo['deltaphi']['nbins']

            # get per-trigger yields for SE
            sesrh1 = sliceAndProjectTHnSparse(sehCorrSR, slices + [jetptslice], AxisNum.deltaphi.value)
            sebrh1 = sliceAndProjectTHnSparse(sehCorrBR, slices + [jetptslice], AxisNum.deltaphi.value)
            sesrh1.Scale(1.0 / nTrigSESR)
            sebrh1.Scale(1.0 / nTrigSEBR)

            # get per-trigger yields for ME
            if centrange == (0, 30):
                # calculate per-trigger yield for 0-10%
                otherslices = []
                otherslices.append((AxisNum.centrality.value, 0, 10))
                otherslices.append((AxisNum.clusterpt.value, ptrange[0], ptrange[1]))
                nTrigsr010 = sliceAndProjectTHnSparse(sehTrigSR, otherslices, 0).Integral()
                nTrigbr010 = sliceAndProjectTHnSparse(sehTrigBR, otherslices, 0).Integral()
                nMixsr010 = sliceAndProjectTHnSparse(mehnMixSR, otherslices, 0).Integral()
                nMixbr010 = sliceAndProjectTHnSparse(mehnMixBR, otherslices, 0).Integral()
                mesr010 = sliceAndProjectTHnSparse(mehCorrSR, otherslices + [jetptslice], AxisNum.deltaphi.value)
                mebr010 = sliceAndProjectTHnSparse(mehCorrBR, otherslices + [jetptslice], AxisNum.deltaphi.value)
                mesr010.Scale(1.0 / nMixsr010)
                mebr010.Scale(1.0 / nMixbr010)

                # calculate per-trigger yield for 10-30%
                otherslices = []
                otherslices.append((AxisNum.centrality.value, 10, 30))
                otherslices.append((AxisNum.clusterpt.value, ptrange[0], ptrange[1]))
                nTrigsr1030 = sliceAndProjectTHnSparse(sehTrigSR, otherslices, 0).Integral()
                nTrigbr1030 = sliceAndProjectTHnSparse(sehTrigBR, otherslices, 0).Integral()
                nMixsr1030 = sliceAndProjectTHnSparse(mehnMixSR, otherslices, 0).Integral()
                nMixbr1030 = sliceAndProjectTHnSparse(mehnMixBR, otherslices, 0).Integral()
                mesr1030 = sliceAndProjectTHnSparse(mehCorrSR, otherslices + [jetptslice], AxisNum.deltaphi.value)
                mebr1030 = sliceAndProjectTHnSparse(mehCorrBR, otherslices + [jetptslice], AxisNum.deltaphi.value)
                mesr1030.Scale(1.0 / nMixsr1030)
                mebr1030.Scale(1.0 / nMixbr1030)

                # combine 0-10% and 10-30%
                mesrh1 = mesr010.Clone()
                mesrh1.Scale(float(nTrigsr010) / nTrigSESR)
                mesrh1.Add(mesr1030, float(nTrigsr1030) / nTrigSESR)
                mebrh1 = mebr010.Clone()
                mebrh1.Scale(float(nTrigbr010) / nTrigSEBR)
                mebrh1.Add(mebr1030, float(nTrigbr1030) / nTrigSEBR)
            else:
                mesrh1 = sliceAndProjectTHnSparse(mehCorrSR, slices + [jetptslice], AxisNum.deltaphi.value)
                mebrh1 = sliceAndProjectTHnSparse(mehCorrBR, slices + [jetptslice], AxisNum.deltaphi.value)
                mesrh1.Scale(1.0 / nMixSR)
                mebrh1.Scale(1.0 / nMixBR)

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
