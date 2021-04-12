import matplotlib.pyplot as plt
import numpy as np

import iminuit

from utils import getBinRange, getCenters, getWidths, normalizeHistAndErr, quadSumPairwise


def getPurity(normsignal, normbkg, binEdges, frac, purityMin, purityMax):
    pmin, pmax = getBinRange(binEdges, purityMin, purityMax)
    signal = frac * normsignal[pmin:pmax]
    bkg = (1 - frac) * normbkg[pmin:pmax]

    return np.sum(signal) / (np.sum(signal) + np.sum(bkg))


def modifyBkgWeights(weights, binEdges, modFunction):
    binCenters = getCenters(binEdges)
    mods = list(map(modFunction, binCenters))
    return np.multiply(weights, mods)


def applyBkgWeights(weights, hist, histerr):
    newhist = np.multiply(hist, weights)
    newhisterr = np.multiply(histerr, weights)
    return normalizeHistAndErr(newhist, newhisterr)


# Requires all histograms as normalized numpy arrays
# In particular, the background weights must already be applied
class TemplateFit:
    def __init__(self, data, dataerr, signal, signalerr, bkg, bkgerr, ssParams, verbosity=0):
        self.data = np.array(data, dtype='f')
        self.dataerr = np.array(dataerr, dtype='f')
        self.signal = np.array(signal, dtype='f')
        self.signalerr = np.array(signalerr, dtype='f')
        self.bkg = np.array(bkg, dtype='f')
        self.bkgerr = np.array(bkgerr, dtype='f')
        self.binEdges = ssParams.binEdges
        self.verbosity = verbosity

        if ssParams.tfFitRange:
            self.fitRange = slice(*getBinRange(self.binEdges, *ssParams.tfFitRange))
        else:
            self.fitRange = slice(0, None)

        self.N = np.sum(self.data)
        self.binCenters = getCenters(self.binEdges)
        self.binWidths = getWidths(self.binEdges)
        self.purityMin = ssParams.purityRange[0]
        self.purityMax = ssParams.purityRange[1]
        self.xlabel = ssParams.axisLabel
        self.signalColor = ssParams.signalColor
        self.bkgColor = ssParams.bkgColor

        self.doFit()
        self.getPurity()

    def doFit(self):
        def Chi2(f):
            model = self.N * (f * self.signal + (1 - f) * self.bkg)
            totalerr = quadSumPairwise(self.dataerr, self.N * (1 - f) * self.bkgerr)
            return np.sum(np.power(np.divide(self.data - model, totalerr, where=np.array(totalerr) != 0), 2.0)[self.fitRange])

        for initf in np.arange(0.10, 1.00, 0.10):
            mt = iminuit.Minuit(Chi2, f=initf, error_f=0.01, limit_f=(0.0, 1.0), errordef=1, print_level=self.verbosity)
            mt.migrad()
            if mt.migrad_ok():
                break

        if not mt.migrad_ok():
            print('Warning: template fit did not converge')

        self.fitf = mt.values['f']
        self.fitferr = mt.errors['f']

        self.fitSignal = self.N * self.fitf * self.signal
        self.fitSignalerr = self.N * self.fitf * self.signalerr
        self.fitBkg = self.N * (1 - self.fitf) * self.bkg
        self.fitBkgerr = self.N * (1 - self.fitf) * self.bkgerr

        fitTotal = self.fitSignal + self.fitBkg
        totalerr = quadSumPairwise(self.dataerr, self.fitBkgerr)
        self.residuals = np.divide(fitTotal - self.data, totalerr, where=np.array(totalerr) != 0)
        self.chi2 = Chi2(self.fitf)
        self.dof = len(self.data[self.fitRange]) - 1  # number of fit parameters

    def getChi2DofInRange(self, rangeMin, rangeMax):
        binmin, binmax = getBinRange(self.binEdges, rangeMin, rangeMax)
        dof = len(self.residuals[binmin:binmax]) - 1
        chi2 = np.sum(np.power(self.residuals[binmin:binmax], 2.0))

        return chi2 / dof

    def getPurity(self):
        self.purity = getPurity(self.signal, self.bkg, self.binEdges, self.fitf, self.purityMin, self.purityMax)
        puritylow = getPurity(self.signal, self.bkg, self.binEdges, self.fitf - self.fitferr, self.purityMin, self.purityMax)
        purityhigh = getPurity(self.signal, self.bkg, self.binEdges, self.fitf + self.fitferr, self.purityMin, self.purityMax)
        self.purityerrlow = self.purity - puritylow
        self.purityerrhigh = purityhigh - self.purity
        self.purityerr = np.mean([self.purityerrlow, self.purityerrhigh])

        if self.verbosity:
            print('Purity = {0:2.2f}, +{1:2.2f}, -{2:2.2f}'.format(self.purity, self.purityerrhigh, self.purityerrlow))

    # returns list of handles: data, bkg, signal, chi2, purity
    def plotFit(self, dataLabel='Data, iso', signalLabel='Signal', bkgLabel='Bkg'):
        norm = self.N * (self.binCenters[1] - self.binCenters[0])
        totalerr = quadSumPairwise(self.dataerr, self.fitBkgerr)

        dataplot = plt.errorbar(self.binCenters, self.data / norm, yerr=self.dataerr / norm, label=dataLabel, fmt='ko')
        bkgplot = plt.bar(self.binCenters, self.fitBkg / norm, width=self.binWidths, align='center', label=bkgLabel, capsize=0, color=self.bkgColor, ec=self.bkgColor)
        signalplot = plt.bar(self.binCenters, self.fitSignal / norm, yerr=totalerr / norm, bottom=self.fitBkg / norm, width=self.binWidths, align='center', label='{0} + {1}'.format(signalLabel, 'Bkg'), capsize=3, color=self.signalColor, ec=self.signalColor, ecolor='gray')
        # signalplot = plt.bar(self.binCenters, self.fitSignal / norm, yerr=totalerr / norm, bottom=self.fitBkg / norm, width=self.binWidths, align='center', label='{0} + {1}'.format(signalLabel, bkgLabel), capsize=3, color=self.signalColor, ec=self.signalColor, ecolor='gray')

        chi2text, = plt.plot([], [], ' ', label='Chi2/dof = {0:2.2f}/{1:2.0f}'.format(self.chi2, self.dof))
        puritytext, = plt.plot([], [], ' ', label='Purity = ${0:2.1f}\pm{1:2.1f}$%'.format(100 * self.purity, 100 * self.purityerr))

        ax = plt.gca()
        pmin, pmax = getBinRange(self.binEdges, self.purityMin, self.purityMax)
        ax.axvspan(self.binCenters[pmin] - self.binWidths[pmin] / 2.0, self.binCenters[pmax] - self.binWidths[pmax] / 2.0, color='black', alpha=0.1)

        return dataplot, bkgplot, signalplot, chi2text, puritytext

    def plotResiduals(self, ylim=[-8.9, 8.9]):
        plt.plot(self.binCenters, self.residuals, 'ko', alpha=0.65)
        plt.ylabel('(Fit - Data)/Dataerr', fontsize=26)
        plt.ylim(ylim)

    def plotNormalizedTemplates(self, signalLabel='Signal', bkgLabel='Bkg'):
        plt.bar(self.binCenters, self.signal, width=self.binWidths, align='center', label=signalLabel, color=self.signalColor, ec=self.signalColor, alpha=0.4)
        plt.bar(self.binCenters, self.bkg, width=self.binWidths, align='center', label=bkgLabel, color=self.bkgColor, ec=self.bkgColor, alpha=0.4)
        plt.xlabel(self.xlabel)
        plt.ylabel('Arb. units')
        plt.legend(loc='best')

    def plotTemplates(self, signalLabel='Signal', bkgLabel='Bkg'):
        plt.bar(self.binCenters, self.fitSignal, width=self.binWidths, align='center', label=signalLabel, color=self.signalColor, ec=self.signalColor, alpha=0.4)
        plt.bar(self.binCenters, self.fitBkg, width=self.binWidths, align='center', label=bkgLabel, color=self.bkgColor, ec=self.bkgColor, alpha=0.4)
        plt.xlabel(self.xlabel)
        plt.ylabel('Entries')
        plt.legend(loc='best')

    def plotFitAndResiduals(self, ptrange, centrange=None, figfilename=None, dataLabel='Data, iso', signalLabel='Signal', bkgLabel='Bkg', system='Pb-Pb', ylim=[-8.9, 8.9]):
        fig = plt.figure()

        fig.add_axes((0.1, 0.3, 0.88, 0.6))

        ax = plt.gca()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='major', length=8, width=2, direction='in')
        ax.tick_params(axis='both', which='minor', length=4, width=1, direction='in')

        fig.add_axes((0.1, 0.3, 0.88, 0.6))
        dataplot, bkgplot, signalplot, chi2text, puritytext = self.plotFit(dataLabel, signalLabel, bkgLabel)

        datapoint, = plt.plot([], [], 'ko', label=dataLabel)
        plt.legend(handles=[datapoint, signalplot, bkgplot, chi2text, puritytext], ncol=1, numpoints=1, loc=1, fontsize=22, frameon=False)

        pttext = '{0} < pT < {1} GeV/$c$'.format(ptrange[0], ptrange[1])
        if centrange:
            centtext = '{0}-{1}% {2}'.format(centrange[0], centrange[1], system)
        infotext = pttext
        if centrange:
            infotext = infotext + '\n' + centtext
        # centtext = '{0}-{1}% Pb-Pb\n$\sqrt{{s_{{NN}}}}=5.02$ TeV'.format(centrange[0], centrange[1])
        plt.annotate(infotext, xy=(0.95, 0.4), xycoords='axes fraction', va='top', ha='right', fontsize=22)
        plt.ylabel('Arbitrary units', fontsize=26, y=1.0, ha='right')

        fig.add_axes((0.1, 0.1, 0.88, 0.2), sharex=ax)

        ax = plt.gca()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='major', length=8, width=2, direction='in')
        ax.tick_params(axis='both', which='minor', length=4, width=1, direction='in')

        self.plotResiduals(ylim=ylim)
        average = np.average(self.residuals)
        plt.axhline(y=average, color='r', label='Average')
        plt.legend(loc=1, frameon=False, fontsize=22)

        plt.xlabel(self.xlabel, fontsize=30, x=1.0, ha='right')

        if figfilename:
            fig.savefig(figfilename)


# Requires all histograms as normalized numpy arrays
# In particular, the background weights must already be applied
class BackgroundFit:
    def __init__(self, data, dataerr, bkg, bkgerr, ssParams, verbosity=0):
        self.data = np.array(data, dtype='f')
        self.dataerr = np.array(dataerr, dtype='f')
        self.bkg = np.array(bkg, dtype='f')
        self.bkgerr = np.array(bkgerr, dtype='f')
        self.binEdges = ssParams.binEdges
        self.fitRange = slice(*getBinRange(self.binEdges, *ssParams.bkgFitRange))
        self.verbosity = verbosity

        self.binCenters = getCenters(self.binEdges)
        self.binWidths = getWidths(self.binEdges)
        self.purityMin = ssParams.purityRange[0]
        self.purityMax = ssParams.purityRange[1]
        self.xlabel = ssParams.axisLabel
        self.bkgColor = ssParams.bkgColor

        self.doFit()
        self.getPurity()

    def doFit(self):
        def Chi2(N):
            model = N * self.bkg
            totalerr = quadSumPairwise(self.dataerr, N * self.bkgerr)
            return np.sum(np.power(np.divide(self.data - model, totalerr, where=np.array(totalerr) != 0), 2.0)[self.fitRange])

        mt = iminuit.Minuit(Chi2, N=np.sum(self.data) / np.sum(self.bkg), error_N=1, errordef=1, print_level=self.verbosity)
        mt.migrad()

        if not mt.migrad_ok():
            print('Warning: template fit did not converge')

        self.fitN = mt.values['N']
        self.fitNerr = mt.errors['N']

        self.fitBkg = self.fitN * self.bkg
        self.fitBkgerr = self.fitN * self.bkgerr

        totalerr = quadSumPairwise(self.dataerr, self.fitBkgerr)

        self.residuals = np.divide(self.fitBkg - self.data, totalerr, where=np.array(totalerr) != 0)
        self.chi2 = Chi2(self.fitN)
        self.dof = len(self.data[self.fitRange]) - 1  # number of fit parameters

    def getChi2DofInRange(self, rangeMin, rangeMax):
        binmin, binmax = getBinRange(self.binEdges, rangeMin, rangeMax)
        dof = len(self.residuals[binmin:binmax]) - 1
        chi2 = np.sum(np.power(self.residuals[binmin:binmax], 2.0))

        return chi2 / dof

    def getPurity(self):
        pmin, pmax = getBinRange(self.binEdges, self.purityMin, self.purityMax)
        fitSignal = np.subtract(self.data, self.fitBkg)
        lowSignal = np.subtract(self.data, (self.fitN + self.fitNerr) * self.bkg)
        highSignal = np.subtract(self.data, (self.fitN - self.fitNerr) * self.bkg)

        self.purity = np.sum(fitSignal[pmin:pmax]) / np.sum(self.data[pmin:pmax])
        puritylow = np.sum(lowSignal[pmin:pmax]) / np.sum(self.data[pmin:pmax])
        purityhigh = np.sum(highSignal[pmin:pmax]) / np.sum(self.data[pmin:pmax])
        self.purityerrlow = self.purity - puritylow
        self.purityerrhigh = purityhigh - self.purity
        self.purityerr = np.mean([self.purityerrlow, self.purityerrhigh])

        if self.verbosity:
            print('Purity = {0:2.2f}, +{1:2.2f}, -{2:2.2f}'.format(self.purity, self.purityerrhigh, self.purityerrlow))

    # returns list of handles: data, bkg, chi2, purity
    def plotFit(self, dataLabel='Data, iso', bkgLabel='Bkg (data, anti-iso)'):
        dataplot = plt.errorbar(self.binCenters, self.data, yerr=self.dataerr, label=dataLabel, fmt='ko')
        bkgplot = plt.bar(self.binCenters, self.fitBkg, yerr=self.fitBkgerr, width=self.binWidths, align='center', label=bkgLabel, capsize=0, color=self.bkgColor, ec=self.bkgColor, ecolor=self.bkgColor)
        chi2text, = plt.plot([], [], ' ', label='Chi2/dof = {0:2.2f}'.format(self.chi2 / self.dof))
        puritytext, = plt.plot([], [], ' ', label='Purity = ${0:2.1f}\pm{1:2.1f}$%'.format(100 * self.purity, 100 * self.purityerr))

        ax = plt.gca()
        pmin, pmax = getBinRange(self.binEdges, self.purityMin, self.purityMax)
        ax.axvspan(self.binEdges[pmin], self.binEdges[pmax], color='black', alpha=0.2)
        ax.axvspan(self.binEdges[self.fitRange.start], self.binEdges[self.fitRange.stop], color='red', alpha=0.2)

        plt.xlabel(self.xlabel)
        plt.ylabel('Entries')

        return [dataplot, bkgplot, chi2text, puritytext]

    def plotResiduals(self, ylim=[-8.9, 8.9]):
        plt.plot(self.binCenters[self.fitRange], self.residuals[self.fitRange], 'ko', alpha=0.65)
        plt.xlabel(self.xlabel)
        plt.ylabel('(Fit - Data)/Dataerr')
        plt.xlim([self.binCenters[0], self.binCenters[-1]])
        plt.ylim(ylim)
