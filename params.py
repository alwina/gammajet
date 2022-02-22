import numpy as np


class IsolationParams:
    # set defaults here or override via keywords
    def __init__(self, isovar='cluster_iso_tpc_02_sub'):
        self.isovar = isovar

    def setFromConfig(self, config):
        if 'isovar' in config:
            self.isovar = config['isovar']
        self.isocut = config['isocut']
        self.antiisocutlow = config['antiisocutlow']
        self.antiisocuthigh = config['antiisocuthigh']

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

    def setFromConfig(self, config):
        self.ssvar = config['ssvar']
        self.binEdges = config['tfbins']
        self.purityRange = (config['srmin'], config['srmax'])
        if 'bkgonlyfitrange' in config:
            self.bkgFitRange = config['bkgonlyfitrange']
        else:
            self.bkgFitRange = (config['brmin'], config['brmax'])
        self.tfFitRange = config.get('tffitrange', None)
        self.axisLabel = config.get('axislabel', self.ssvar)
        self.legendLabel = config.get('legendlabel', self.ssvar)
        self.doubleRatioBinEdges = config['doubleratiobinedges']
        self.doubleRatioFitRange = config['doubleratiofitrange']

    def setDefaultLambda(self):
        self.ssvar = 'cluster_Lambda'
        self.binEdges = np.linspace(0.1, 2, 96)
        self.purityRange = (0.1, 0.3)
        self.bkgFitRange = (0.6, 1.5)
        self.tfFitRange = None
        self.axisLabel = '$\mathrm{\sigma^2_{long}}$'
        self.legendLabel = '$\mathrm{\sigma^2_{long}}$'
        self.doubleRatioBinEdges = np.linspace(0.1, 2, 20)
        self.doubleRatioFitRange = (0.5, 1.75)

    def setVarBinLambda(self):
        self.setDefaultLambda()
        self.binEdges = np.concatenate((np.linspace(0.1, 0.64, 28), np.linspace(0.72, 1.12, 6), np.linspace(1.2, 2, 6)))

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


def ptcuttext(ptrange):
    return '{1} <= {0} and {0} < {2}'.format('cluster_pt', *ptrange)


def centralitycuttext(centrange):
    return '{1} <= {0} and {0} < {2}'.format('centrality_v0m', *centrange)


def parseCut(cutvar, cuttype, cutval):
    if cuttype == 'min':
        cuttext = '{0} > {1}'.format(cutvar, cutval)
    elif cuttype == 'max':
        cuttext = '{0} < {1}'.format(cutvar, cutval)
    elif cuttype == 'incmin':
        cuttext = '{0} >= {1}'.format(cutvar, cutval)
    elif cuttype == 'incmax':
        cuttext = '{0} <= {1}'.format(cutvar, cutval)
    elif cuttype == 'equals':
        cuttext = '{0} == {1}'.format(cutvar, cutval)
    elif cuttype == 'absmax':
        cuttext = 'abs({0}) < {1}'.format(cutvar, cutval)
    else:
        print('Cut type {0} not recognized'.format(cutval))
        cuttext = ''

    return cuttext
