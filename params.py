import numpy as np


class IsolationParams:
    # set defaults here or override via keywords
    def __init__(self, **kwargs):
        self.isovar = kwargs.get('isovar', 'cluster_iso_its_04_sub')
        self.isocut = kwargs.get('isocut', 1.5)
        self.antiisocutlow = kwargs.get('antiisocutlow', 5.0)
        self.antiisocuthigh = kwargs.get('antiisocuthigh', 10.0)

    def isocuttext(self):
        return '{0} < {1}'.format(self.isovar, self.isocut)

    def antiisocuttext(self):
        return '{0} < {1} and {1} < {2}'.format(self.antiisocutlow, self.isovar, self.antiisocuthigh)

    def antiisocutrange(self):
        return (self.antiisocutlow, self.antiisocuthigh)

    def copy(self, inputIsoParams):
        self.isovar = inputIsoParams.isovar
        self.isocut = inputIsoParams.isocut
        self.antiisocutlow = inputIsoParams.antiisocutlow
        self.antiisocuthigh = inputIsoParams.antiisocuthigh


class ShowerShapeParams:
    def __init__(self):
        self.signalColor = '#3B7EA1'
        self.bkgColor = '#FDB515'

    def setDefaultLambda(self):
        self.ssvar = 'cluster_Lambda'
        self.binEdges = np.linspace(0, 2, 101)
        self.purityRange = (0.0, 0.3)
        self.bkgFitRange = (0.4, 1.5)
        self.tfFitRange = None
        self.axisLabel = '$\mathrm{\sigma^2_{long}}$'
        self.legendLabel = '$\mathrm{\sigma^2_{long}}$'
        self.doubleRatioBinEdges = np.linspace(0, 2, 21)
        self.doubleRatioFitRange = (0.5, 1.75)

    def setDefaultDNN(self):
        self.ssvar = 'cluster_NN1'
        self.binEdges = np.linspace(0, 1, 101)
        self.purityRange = (0.7, 0.9)
        self.bkgFitRange = (0.0, 0.3)
        self.tfFitRange = None
        self.axisLabel = 'Deep Neural Network Output'
        self.legendLabel = 'DNN'
        self.doubleRatioBinEdges = np.linspace(0, 1, 21)
        self.doubleRatioFitRange = (0.0, 0.6)


def cuttext(cutvar, cutrange):
    return '{1} < {0} and {0} < {2}'.format(cutvar, *cutrange)


def ptcuttext(ptrange):
    return cuttext('cluster_pt', ptrange)


def centralitycuttext(centrality):
    return cuttext('centrality_v0m', centrality)
