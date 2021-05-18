import matplotlib.pyplot as plt
import ROOT

from utils import plotTH1, plotTH2

ROOT.gSystem.Load("/home/alwina/RooUnfold/libRooUnfold")


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
        plt.xlabel('Truth')
        plt.ylabel('Measured')
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
