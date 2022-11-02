import matplotlib.pyplot as plt
import numpy as np
import ROOT

from utils import plotTH1, plotTH2, th1ToArrays, sliceAndProjectTHnSparse

ROOT.gSystem.Load("/home/alwina/RooUnfold/libRooUnfold")


class Unfolder2D:
    def __init__(self, verbosity=1):
        # keys are regularization parameter values
        # values are TH1s
        self.unfoldedBayesTH2 = {}
        self.unfoldedSvdTH2 = {}
        self.verbosity = verbosity

    def setMeasuredTH2(self, measuredTH2):
        self.measuredTH2 = measuredTH2

    def setResponseMatrix(self, rooUnfoldResponse):
        self.responseMatrix = rooUnfoldResponse

    def setTH4(self, th4):
        # this is a THnSparse, not a TH4, but whatever
        self.TH4 = th4

    def unfoldBayes(self, *args):
        for nIterations in args:
            unfold = ROOT.RooUnfoldBayes(self.responseMatrix, self.measuredTH2, nIterations)
            unfold.SetVerbose(self.verbosity)
            self.unfoldedBayesTH2[nIterations] = unfold.Hreco()

    def unfoldSvd(self, *args):
        # in the 1D case, k should be between 2 and the number of bins
        # not sure what that means for 2D
        for k in args:
            unfold = ROOT.RooUnfoldSvd(self.responseMatrix, self.measuredTH2, k)
            self.unfoldedSvdTH2[k] = unfold.Hreco()

    def plotAllUnfolded(self, bayesKeys=None, svdKeys=None, axis=1):
        if axis not in (1, 2):
            print('axis should be 1 (X) or 2 (Y), not {0}'.format(axis))
            return

        if axis == 1:
            plotTH1(self.measuredTH2.ProjectionX(), fmt='ko', label='Measured', divideBinWidths=True, skipZeros=True)
        else:
            plotTH1(self.measuredTH2.ProjectionY(), fmt='ko', label='Measured', divideBinWidths=True, skipZeros=True)

        if bayesKeys is None:
            bayesKeys = sorted(self.unfoldedBayesTH2.keys())
        for key in bayesKeys:
            if axis == 1:
                plotTH1(self.unfoldedBayesTH2[key].ProjectionX(), fmt='o-', label='Bayesian {0}'.format(key), divideBinWidths=True, skipZeros=True)
            else:
                plotTH1(self.unfoldedBayesTH2[key].ProjectionY(), fmt='o-', label='Bayesian {0}'.format(key), divideBinWidths=True, skipZeros=True)

        if svdKeys is None:
            svdKeys = sorted(self.unfoldedSvdTH2.keys())
        for key in svdKeys:
            if axis == 1:
                plotTH1(self.unfoldedSvdTH2[key].ProjectionX(), fmt='o-', label='SVD {0}'.format(key), divideBinWidths=True, skipZeros=True)
            else:
                plotTH1(self.unfoldedSvdTH2[key].ProjectionY(), fmt='o-', label='SVD {0}'.format(key), divideBinWidths=True, skipZeros=True)

    def plotBayesConvergence(self, bayesKeys=None, axis=1):
        if axis not in (1, 2):
            print('axis should be 1 (X) or 2 (Y), not {0}'.format(axis))
            return

        if bayesKeys is None:
            bayesKeys = self.unfoldedBayesTH2.keys()

        sortedkeys = sorted(bayesKeys)
        if axis == 1:
            denth1 = self.measuredTH2.ProjectionX()
        else:
            denth1 = self.measuredTH2.ProjectionY()
        denkey = 0
        for key in sortedkeys:
            numkey = key
            if axis == 1:
                numth1 = self.unfoldedBayesTH2[key].ProjectionX()
            else:
                numth1 = self.unfoldedBayesTH2[key].ProjectionY()
            numth1.Divide(denth1)

            hist, err, binCenters, binWidths = th1ToArrays(numth1)
            plt.plot(binCenters, hist, 'o-', label='{0}/{1}'.format(numkey, denkey))

            if axis == 1:
                denth1 = self.unfoldedBayesTH2[key].ProjectionX()
            else:
                denth1 = self.unfoldedBayesTH2[key].ProjectionY()

            denkey = key

    def shapeClosureTestMCDataRatio(self, bayesKeys, axis=1):
        if axis not in (1, 2):
            print('axis should be 1 (X) or 2 (Y), not {0}'.format(axis))
            return

        recoth2 = sliceAndProjectTHnSparse(self.TH4, [], 2, 0)
        truthth2 = sliceAndProjectTHnSparse(self.TH4, [], 3, 1)

        recoth2.RebinX(int(120 / self.measuredTH2.GetNbinsX()))
        recoth2.RebinY(int(120 / self.measuredTH2.GetNbinsY()))
        truthth2.RebinX(int(120 / self.measuredTH2.GetNbinsX()))
        truthth2.RebinY(int(120 / self.measuredTH2.GetNbinsY()))

        # scale factor: MC reco / data
        scaleth2 = recoth2.Clone()
        # don't propagate the error from the measured in this case
        measurednoerror = self.measuredTH2.Clone()
        for binx in range(1, measurednoerror.GetNbinsX() + 1):
            for biny in range(1, measurednoerror.GetNbinsY() + 1):
                measurednoerror.SetBinError(binx, biny, 0)
        scaleth2.Divide(measurednoerror)

        # scale MC truth and reco
        recoth2.Multiply(scaleth2)
        truthth2.Multiply(scaleth2)

        # make integral match measured
        measuredintegral = abs(self.measuredTH2.Integral())
        recointegral = recoth2.Integral()
        truthintegral = truthth2.Integral()

        recoth2.Scale(measuredintegral / recointegral)
        truthth2.Scale(measuredintegral / truthintegral)

        # smear by statistical errors on measured data
        for binx in range(1, self.measuredTH2.GetNbinsX() + 1):
            for biny in range(1, self.measuredTH2.GetNbinsY() + 1):
                if self.measuredTH2.GetBinContent(binx, biny) != 0:
                    relerr = self.measuredTH2.GetBinError(binx, biny) / self.measuredTH2.GetBinContent(binx, biny)
                    smearedbin = recoth2.GetBinContent(binx, biny) * np.random.normal(1.0, abs(relerr))
                    recoth2.SetBinContent(binx, biny, smearedbin)
                    # make the relative error match the measured
                    recoth2.SetBinError(binx, biny, smearedbin * relerr)

        # unfold this new thing - this shouldn't be infinitely recursive, right???
        unfolder = Unfolder2D(0)
        unfolder.setMeasuredTH2(recoth2)
        unfolder.setResponseMatrix(self.responseMatrix)
        unfolder.unfoldBayes(*bayesKeys)

        # compare to scaled truth
        if axis == 1:
            plotTH1(truthth2.ProjectionX(), fmt='ks', ms=8, label='Scaled MC truth')
            for key in bayesKeys:
                plotTH1(unfolder.unfoldedBayesTH2[key].ProjectionX(), fmt='o-', label='Unfolded {0}'.format(key))
        else:
            plotTH1(truthth2.ProjectionY(), fmt='ks', ms=8, label='Scaled MC truth')
            for key in bayesKeys:
                plotTH1(unfolder.unfoldedBayesTH2[key].ProjectionY(), fmt='o-', label='Unfolded {0}'.format(key))

        return unfolder, truthth2

    def statisticalClosureTest(self, bayesKeys, axis=1):
        if axis not in (1, 2):
            print('axis should be 1 (X) or 2 (Y), not {0}'.format(axis))
            return

        recoth2 = sliceAndProjectTHnSparse(self.TH4, [], 2, 0)
        truthth2 = sliceAndProjectTHnSparse(self.TH4, [], 3, 1)

        recoth2.RebinX(int(120 / self.measuredTH2.GetNbinsX()))
        recoth2.RebinY(int(120 / self.measuredTH2.GetNbinsY()))
        truthth2.RebinX(int(120 / self.measuredTH2.GetNbinsX()))
        truthth2.RebinY(int(120 / self.measuredTH2.GetNbinsY()))

        # make integral match measured
        measuredintegral = abs(self.measuredTH2.Integral())
        recointegral = recoth2.Integral()
        truthintegral = truthth2.Integral()

        recoth2.Scale(measuredintegral / recointegral)
        truthth2.Scale(measuredintegral / truthintegral)
        unsmearedrecoth2 = recoth2.Clone()

        # smear by statistical errors on measured data
        for binx in range(1, self.measuredTH2.GetNbinsX() + 1):
            for biny in range(1, self.measuredTH2.GetNbinsY() + 1):
                if self.measuredTH2.GetBinContent(binx, biny) != 0:
                    relerr = self.measuredTH2.GetBinError(binx, biny) / self.measuredTH2.GetBinContent(binx, biny)
                    smearedbin = recoth2.GetBinContent(binx, biny) * np.random.normal(1.0, abs(relerr))
                    recoth2.SetBinContent(binx, biny, smearedbin)
                    # make the relative error match the measured
                    recoth2.SetBinError(binx, biny, smearedbin * relerr)

        # unfold this new thing - this shouldn't be infinitely recursive, right???
        unfolder = Unfolder2D(0)
        unfolder.setMeasuredTH2(recoth2)
        unfolder.setResponseMatrix(self.responseMatrix)
        unfolder.unfoldBayes(*bayesKeys)

        # compare to truth
        if axis == 1:
            plotTH1(truthth2.ProjectionX(), fmt='ko-', linewidth=2, ms=10, label='MC truth')
            plotTH1(self.measuredTH2.ProjectionX(), fmt='s', mfc='None', mew=2, label='Measured')
            plotTH1(unsmearedrecoth2.ProjectionX(), fmt='d-', label='Unsmeared MC reco')
            plotTH1(recoth2.ProjectionX(), fmt='d-', label='Smeared MC reco')
            for key in bayesKeys:
                plotTH1(unfolder.unfoldedBayesTH2[key].ProjectionX(), fmt='o-', ms=10, mfc='None', label='Unfolded {0}'.format(key))
        else:
            plotTH1(truthth2.ProjectionY(), fmt='ko-', linewidth=2, ms=10, label='MC truth')
            plotTH1(self.measuredTH2.ProjectionY(), fmt='s', mfc='None', mew=2, label='Measured')
            plotTH1(unsmearedrecoth2.ProjectionY(), fmt='d-', label='Unsmeared MC reco')
            plotTH1(recoth2.ProjectionY(), fmt='d-', label='Smeared MC reco')
            for key in bayesKeys:
                plotTH1(unfolder.unfoldedBayesTH2[key].ProjectionY(), fmt='o-', ms=10, mfc='None', label='Unfolded {0}'.format(key))


class Unfolder:
    def __init__(self):
        # keys are regularization parameter values
        # values are TH1s
        self.unfoldedBayesTH1 = {}
        self.unfoldedSvdTH1 = {}

    def setMeasuredTH1(self, measuredTH1):
        self.measuredTH1 = measuredTH1

    def setResponseMatrix(self, rooUnfoldResponse):
        self.responseMatrix = rooUnfoldResponse

    def plotResponseMatrix(self):
        responseTH2 = self.responseMatrix.Hresponse()
        plotTH2(responseTH2)
        plt.xlabel('Measured')
        plt.ylabel('Truth')
        plt.colorbar()

    def unfoldBayes(self, *args):
        for nIterations in args:
            unfold = ROOT.RooUnfoldBayes(self.responseMatrix, self.measuredTH1, nIterations)
            self.unfoldedBayesTH1[nIterations] = unfold.Hreco()

    def unfoldSvd(self, *args):
        # k should be between 2 and the number of bins
        for k in args:
            unfold = ROOT.RooUnfoldSvd(self.responseMatrix, self.measuredTH1, k)
            self.unfoldedSvdTH1[k] = unfold.Hreco()

    def plotAllUnfolded(self, bayesKeys=None, svdKeys=None):
        plotTH1(self.measuredTH1, fmt='ko', label='Measured')

        if bayesKeys is None:
            bayesKeys = self.unfoldedBayesTH1.keys()
        for key in bayesKeys:
            plotTH1(self.unfoldedBayesTH1[key], label='Bayesian {0}'.format(key), linestyle='None', marker='.')

        if svdKeys is None:
            svdKeys = self.unfoldedSvdTH1.keys()
        for key in svdKeys:
            plotTH1(self.unfoldedSvdTH1[key], label='SVD {0}'.format(key), linestyle='None', marker='.')
