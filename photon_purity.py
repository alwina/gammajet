import numpy as np

from fit_functions import SingleParameterLinearFit
from uncertainty_background_template_correction import getBkgCorrectionUncertainty
from uncertainty_sideband_selection import getSidebandSelectionUncertainty
from uncertainty_signal_template import getBfTfDiff
from utils import quadSumAll


class PhotonPurity:
    def __init__(self, inputValues=None):
        if inputValues:
            for key in inputValues:
                setattr(self, key, inputValues[key])

    def computePurityAndUncertainties(self, dfs, ptrange, isoParams, ssParams, computeSystematics=True, verbose=False):
        bkgWeights = dfs.getBkgWeights(ptrange, isoParams, ssParams)
        tf = dfs.getTemplateFit(ptrange, isoParams, ssParams, bkgWeights=bkgWeights)
        self.purity = tf.purity
        self.staterr = tf.purityerr

        if verbose:
            print 'Purity = {0:2.1f}% +/- {1:2.1f}%'.format(100 * self.purity, 100 * self.staterr)

        if computeSystematics:
            self.signal = getBfTfDiff(dfs, ptrange, isoParams, ssParams)
            if verbose:
                print 'Signal template uncertainty: {0:2.2f}%'.format(100 * self.signal)

            sidebandWindows = [(x, x + 2.0) for x in np.arange(-3, 25, 0.5)]
            self.sideband = getSidebandSelectionUncertainty(dfs, ptrange, isoParams, ssParams, sidebandWindows)
            if verbose:
                print 'Sideband selection uncertainty: {0:2.2f}%'.format(100 * self.sideband)

            doubleRatioFit = SingleParameterLinearFit()
            self.bkgcorr = getBkgCorrectionUncertainty(dfs, ptrange, isoParams, ssParams, doubleRatioFit)
            if verbose:
                print 'Background correction uncertainty: {0:2.2f}%'.format(100 * self.bkgcorr)

            self.systerr = quadSumAll([self.signal, self.sideband, self.bkgcorr])

    # mostly useful for pickling
    def getDictionary(self, verbose=False):
        d = {}
        d['purity'] = self.purity
        d['staterr'] = self.staterr
        d['signal'] = self.signal
        d['sideband'] = self.sideband
        d['bkgcorr'] = self.bkgcorr
        d['systerr'] = self.systerr
        if verbose:
            print d
        return d
