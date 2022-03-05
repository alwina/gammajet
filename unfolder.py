import matplotlib.pyplot as plt
import ROOT

from utils import plotTH1, plotTH2, th1ToArrays

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
        if axis == 1:
            plotTH1(self.measuredTH2.ProjectionX(), fmt='ko', label='Measured')

            if bayesKeys is None:
                bayesKeys = sorted(self.unfoldedBayesTH2.keys())
            for key in bayesKeys:
                plotTH1(self.unfoldedBayesTH2[key].ProjectionX(), fmt='o-', label='Bayesian {0}'.format(key))

            if svdKeys is None:
                svdKeys = sorted(self.unfoldedSvdTH2.keys())
            for key in svdKeys:
                plotTH1(self.unfoldedSvdTH2[key].ProjectionX(), fmt='o-', label='Bayesian {0}'.format(key))
        elif axis == 2:
            plotTH1(self.measuredTH2.ProjectionY(), fmt='ko', label='Measured')

            if bayesKeys is None:
                bayesKeys = sorted(self.unfoldedBayesTH2.keys())
            for key in bayesKeys:
                plotTH1(self.unfoldedBayesTH2[key].ProjectionY(), fmt='o-', label='Bayesian {0}'.format(key))

            if svdKeys is None:
                svdKeys = sorted(self.unfoldedSvdTH2.keys())
            for key in svdKeys:
                plotTH1(self.unfoldedSvdTH2[key].ProjectionY(), fmt='o-', label='Bayesian {0}'.format(key))
        else:
            print('axis should be 1 (X) or 2 (Y), not {0}'.format(axis))

    def plotBayesConvergence(self, bayesKeys=None, axis=1):
        if bayesKeys is None:
            bayesKeys = self.unfoldedBayesTH2.keys()

        sortedkeys = sorted(bayesKeys)
        if axis == 1:
            denth1 = self.measuredTH2.ProjectionX()
            denkey = 0
            for key in sortedkeys:
                numkey = key
                numth1 = self.unfoldedBayesTH2[key].ProjectionX()
                numth1.Divide(denth1)

                hist, err, binCenters, binWidths = th1ToArrays(numth1)
                plt.plot(binCenters, hist, 'o-', label='{0}/{1}'.format(numkey, denkey))

                denth1 = self.unfoldedBayesTH2[key].ProjectionX()
                denkey = key
        elif axis == 2:
            denth1 = self.measuredTH2.ProjectionY()
            denkey = 0
            for key in sortedkeys:
                numkey = key
                numth1 = self.unfoldedBayesTH2[key].ProjectionY()
                numth1.Divide(denth1)

                hist, err, binCenters, binWidths = th1ToArrays(numth1)
                plt.plot(binCenters, hist, 'o-', label='{0}/{1}'.format(numkey, denkey))

                denth1 = self.unfoldedBayesTH2[key].ProjectionY()
                denkey = key
        else:
            print('axis should be 1 (X) or 2 (Y), not {0}'.format(axis))


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
