import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import iminuit

from utils import getBinRange, getCenters, getWidths, normalizeHistAndErr, quadSumPairwise, getHistAndErr, getNormHistAndErr, divideHistsAndErrs
from plotstyle import upperRightText
from params import centralitycuttext, ptcuttext
from uncertainty_background_template_correction import getDoubleRatioAndError, getDoubleRatioFitAndError, plotDoubleRatioAndFit
from fit_functions import SingleParameterLinearFit


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

        for initf in [0.1, 0.2, 0.4, 0.8, 0.05, 0.02, 0.01, 0.005]:
            mt = iminuit.Minuit(Chi2, f=initf, error_f=0.01, limit_f=(0.0, 1.0), errordef=1, print_level=self.verbosity)
            mt.migrad()

            if not np.isnan(mt.values['f']) and mt.migrad_ok():
                break

        if np.isnan(mt.values['f']) or not mt.migrad_ok():
            print('Warning: template fit did not converge')
            print('Data:', self.data)
            print('Signal:', self.signal)
            print('Bkg:', self.bkg)
            print('Data err:', self.dataerr)
            print('Bkg err:', self.bkgerr)

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
        puritytext, = plt.plot([], [], ' ', label='Purity = ${0:2.1f} \pm {1:2.1f}$%'.format(100 * self.purity, 100 * self.purityerr))

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
        puritytext, = plt.plot([], [], ' ', label='Purity = ${0:2.1f} \pm {1:2.1f}$%'.format(100 * self.purity, 100 * self.purityerr))

        ax = plt.gca()
        pmin, pmax = getBinRange(self.binEdges, self.purityMin, self.purityMax)
        ax.axvspan(self.binEdges[pmin], self.binEdges[pmax], color='black', alpha=0.2)
        ax.axvspan(self.binEdges[self.fitRange.start], self.binEdges[self.fitRange.stop], color='red', alpha=0.2)

        plt.xlabel(self.xlabel)
        plt.ylabel('Entries')

        return dataplot, bkgplot, chi2text, puritytext

    def plotResiduals(self, ylim=[-8.9, 8.9]):
        plt.plot(self.binCenters[self.fitRange], self.residuals[self.fitRange], 'ko', alpha=0.65)
        plt.xlabel(self.xlabel)
        plt.ylabel('(Fit - Data)/Dataerr')
        plt.xlim([self.binCenters[0], self.binCenters[-1]])
        plt.ylim(ylim)

    def plotFitAndResiduals(self, ptrange, centrange=None, figfilename=None, dataLabel='Data, iso', bkgLabel='Bkg', system='Pb-Pb', ylim=[-8.9, 8.9]):
        fig = plt.figure()

        fig.add_axes((0.1, 0.3, 0.88, 0.6))

        ax = plt.gca()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='major', length=8, width=2, direction='in')
        ax.tick_params(axis='both', which='minor', length=4, width=1, direction='in')

        fig.add_axes((0.1, 0.3, 0.88, 0.6))
        dataplot, bkgplot, chi2text, puritytext = self.plotFit(dataLabel, bkgLabel)

        datapoint, = plt.plot([], [], 'ko', label=dataLabel)
        plt.legend(handles=[datapoint, bkgplot, chi2text, puritytext], ncol=1, numpoints=1, loc=1, fontsize=22, frameon=False)

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


class ExtendedTemplateFit:
    def __init__(self, dfs, ptrange, centrange, isoParams, ssParams, usezeta=False, useraa=True, usedoubleratio=False, verbosity=0):
        bincut = '{0} and {1}'.format(centralitycuttext(centrange), ptcuttext(ptrange))
        self.ptrange = ptrange
        self.centrange = centrange
        self.isoParams = isoParams
        self.ssParams = ssParams

        self.dfs = dfs
        isodatadf = dfs.fulldatadf.query(isoParams.isocuttext()).query(bincut)
        antiisodatadf = dfs.fulldatadf.query(isoParams.antiisocuttext()).query(bincut)
        isogjmcdf = dfs.fullgjmcdf.query(isoParams.isocuttext()).query(bincut)
        antiisogjmcdf = dfs.fullgjmcdf.query(isoParams.antiisocuttext()).query(bincut)
        isojjmcdf = dfs.fulljjmcdf.query(isoParams.isocuttext()).query(bincut)
        antiisojjmcdf = dfs.fulljjmcdf.query(isoParams.antiisocuttext()).query(bincut)

        self.isodatahist, self.isodataerr = getHistAndErr(isodatadf, ssParams.ssvar, ssParams.binEdges)
        self.isogjmchist, self.isogjmcerr = getNormHistAndErr(isogjmcdf, ssParams.ssvar, ssParams.binEdges)
        self.isojjmchist, self.isojjmcerr = getNormHistAndErr(isojjmcdf, ssParams.ssvar, ssParams.binEdges)
        self.antiisodatahist, self.antiisodataerr = getNormHistAndErr(antiisodatadf, ssParams.ssvar, ssParams.binEdges)
        self.antiisogjmchist, self.antiisogjmcerr = getNormHistAndErr(antiisogjmcdf, ssParams.ssvar, ssParams.binEdges)
        self.antiisojjmchist, self.antiisojjmcerr = getNormHistAndErr(antiisojjmcdf, ssParams.ssvar, ssParams.binEdges)

        if useraa:
            # make combined histogram with anti-isolated MCs using RAA-scaled weights
            antiisomcdf = pd.concat([antiisogjmcdf, antiisojjmcdf])
            self.antiisomchist, self.antiisomcerr = getNormHistAndErr(antiisomcdf, ssParams.ssvar, ssParams.binEdges, weightvar='weightswithraa')
            # also scale the isolated dijet MC with RAA
            self.isojjmchist, self.isojjmcerr = getNormHistAndErr(isojjmcdf, ssParams.ssvar, ssParams.binEdges, weightvar='weightswithraa')

        if ssParams.tfFitRange:
            self.fitRange = slice(*getBinRange(self.binEdges, *ssParams.tfFitRange))
        else:
            self.fitRange = slice(0, None)

        self.N = np.sum(self.isodatahist)
        self.binEdges = ssParams.binEdges
        self.binCenters = getCenters(self.binEdges)
        self.binWidths = getWidths(self.binEdges)
        self.purityMin = ssParams.purityRange[0]
        self.purityMax = ssParams.purityRange[1]
        self.doubleRatioFitRange = ssParams.doubleRatioFitRange
        self.xlabel = ssParams.axisLabel
        self.signalColor = ssParams.signalColor
        self.bkgColor = ssParams.bkgColor
        self.verbosity = verbosity
        self.drfit = None

        # flags for how to correct the background template
        self.usezeta = usezeta
        self.useraa = useraa
        self.usedoubleratio = usedoubleratio

        if usedoubleratio:
            self.getDoubleRatioFit()

        self.doFit()
        self.getPurity()

    def getDoubleRatioFit(self):
        self.doubleratio, self.doubleratioerr = getDoubleRatioAndError(self.dfs, self.ptrange, self.isoParams, self.ssParams, self.centrange, useraa=self.useraa)
        self.drfit = SingleParameterLinearFit()
        getDoubleRatioFitAndError(self.doubleratio, self.doubleratioerr, self.ssParams, self.drfit)

    def getBkgCorrectionWeights(self):
        if self.useraa:
            bkgweights, bkgweightserr = divideHistsAndErrs(self.isojjmchist, self.isojjmcerr, self.antiisomchist, self.antiisomcerr)
        else:
            bkgweights, bkgweightserr = divideHistsAndErrs(self.isojjmchist, self.isojjmcerr, self.antiisojjmchist, self.antiisojjmcerr)

        if self.usedoubleratio:
            bkgweights = modifyBkgWeights(bkgweights, self.binEdges, self.drfit.getFunctionFit())
            bkgweightserr = modifyBkgWeights(bkgweightserr, self.binEdges, self.drfit.getFunctionFit())

        return bkgweights, bkgweightserr

    def getBkgTemplate(self, zeta):
        # background template starts with the anti-isolated data
        bkgtemplate = self.antiisodatahist
        bkgtemplateerr = self.antiisodataerr

        # then subtract the anti-isolated prompt photons, scaled by zeta
        bkgtemplate = np.subtract(bkgtemplate, zeta * self.antiisogjmchist)
        # bkgtemplateerr = quadSumPairwise(bkgtemplateerr, zeta * self.antiisogjmcerr)

        # renormalize to 1
        bkgtemplate, bkgtemplateerr = normalizeHistAndErr(bkgtemplate, bkgtemplateerr)

        # get the background weights
        bkgweights, bkgweightserr = self.getBkgCorrectionWeights()

        # apply background weights
        bkgtemplate, bkgtemplateerr = applyBkgWeights(bkgweights, bkgtemplate, bkgtemplateerr)

        return bkgtemplate, bkgtemplateerr

    def doFit(self):
        if self.usezeta:
            def Chi2(f, zeta):
                sigtemplate = self.isogjmchist
                bkgtemplate, bkgtemplateerr = self.getBkgTemplate(zeta)

                model = self.N * (f * sigtemplate + (1 - f) * bkgtemplate)
                totalerr = quadSumPairwise(self.isodataerr, self.N * (1 - f) * bkgtemplateerr)
                return np.sum(np.power(np.divide(self.isodatahist - model, totalerr, where=np.array(totalerr) != 0), 2.0)[self.fitRange])

            for initf in np.arange(0.1, 1.0, 0.1):
                initzeta = 0.1
                mt = iminuit.Minuit(Chi2, f=initf, error_f=0.01, limit_f=(0.0, 1.0), zeta=initzeta, error_zeta=initzeta / 10, limit_zeta=(0.0, 1.0), errordef=1, print_level=self.verbosity)
                mt.migrad()

                if mt.migrad_ok() and not np.isnan(mt.values['f']):
                    break
        else:
            def Chi2(f):
                sigtemplate = self.isogjmchist
                bkgtemplate, bkgtemplateerr = self.getBkgTemplate(0)

                model = self.N * (f * sigtemplate + (1 - f) * bkgtemplate)
                totalerr = quadSumPairwise(self.isodataerr, self.N * (1 - f) * bkgtemplateerr)
                return np.sum(np.power(np.divide(self.isodatahist - model, totalerr, where=np.array(totalerr) != 0), 2.0)[self.fitRange])

            for initf in [0.1, 0.2, 0.4, 0.8, 0.05, 0.02, 0.01, 0.005]:
                mt = iminuit.Minuit(Chi2, f=initf, error_f=0.01, limit_f=(0.0, 1.0), errordef=1, print_level=self.verbosity)
                mt.migrad()

                if not np.isnan(mt.values['f']) and mt.migrad_ok():
                    break

        if np.isnan(mt.values['f']) or not mt.migrad_ok():
            print('Warning: template fit did not converge')
            print('Data:', self.isodatahist)
            print('Signal:', self.isogjmchist)
            print('Anti-iso data:', self.antiisodatahist)
            print('Iso dijet MC:', self.isojjmchist)
            print('Anti-iso MC:', self.antiisomchist)

        self.fitf = mt.values['f']
        self.fitferr = mt.errors['f']
        if self.usezeta:
            self.fitzeta = mt.values['zeta']
            self.fitzetaerr = mt.values['zeta']
        else:
            self.fitzeta = 0
            self.fitzetaerr = 0

        bkgtemplate, bkgtemplateerr = self.getBkgTemplate(self.fitzeta)

        self.fitSignal = self.N * self.fitf * self.isogjmchist
        self.fitSignalerr = self.N * self.fitf * self.isogjmcerr
        self.fitBkg = self.N * (1 - self.fitf) * bkgtemplate
        self.fitBkgerr = self.N * (1 - self.fitf) * bkgtemplateerr

        fitTotal = self.fitSignal + self.fitBkg
        totalerr = quadSumPairwise(self.isodataerr, self.fitBkgerr)
        self.residuals = np.divide(fitTotal - self.isodatahist, totalerr, where=np.array(totalerr) != 0)
        if self.usezeta:
            self.chi2 = Chi2(self.fitf, self.fitzeta)
            self.dof = len(self.isodatahist[self.fitRange]) - 2
        else:
            self.chi2 = Chi2(self.fitf)
            self.dof = len(self.isodatahist[self.fitRange]) - 1

    def getPurity(self):
        bkgtemplate, _ = self.getBkgTemplate(self.fitzeta)
        self.purity = getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf, self.purityMin, self.purityMax)

        # vary fit parameters and take absolute max/mins
        purityvars = []
        # central zeta
        purityvars.append(getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf - self.fitferr, self.purityMin, self.purityMax))
        purityvars.append(getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf + self.fitferr, self.purityMin, self.purityMax))
        # zeta down
        bkgtemplate, _ = self.getBkgTemplate(self.fitzeta - self.fitzetaerr)
        purityvars.append(getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf, self.purityMin, self.purityMax))
        purityvars.append(getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf - self.fitferr, self.purityMin, self.purityMax))
        purityvars.append(getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf + self.fitferr, self.purityMin, self.purityMax))
        # zeta up
        bkgtemplate, _ = self.getBkgTemplate(self.fitzeta + self.fitzetaerr)
        purityvars.append(getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf, self.purityMin, self.purityMax))
        purityvars.append(getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf - self.fitferr, self.purityMin, self.purityMax))
        purityvars.append(getPurity(self.isogjmchist, bkgtemplate, self.binEdges, self.fitf + self.fitferr, self.purityMin, self.purityMax))

        self.purityerrlow = self.purity - min(purityvars)
        self.purityerrhigh = max(purityvars) - self.purity
        self.purityerr = np.mean([self.purityerrlow, self.purityerrhigh])
        self.purityvars = purityvars

    def plotFit(self, dataLabel='Data, iso', signalLabel='Signal', bkgLabel='Bkg'):
        norm = self.N
        totalerr = quadSumPairwise(self.isodataerr, self.fitBkgerr)

        dataplot = plt.errorbar(self.binCenters, np.divide(self.isodatahist, norm), yerr=np.divide(self.isodataerr, norm), label=dataLabel, fmt='ko')
        bkgplot = plt.bar(self.binCenters, np.divide(self.fitBkg, norm), width=self.binWidths, align='center', label=bkgLabel, capsize=0, color=self.bkgColor, ec=self.bkgColor)
        sigplot = plt.bar(self.binCenters, np.divide(self.fitSignal, norm), yerr=np.divide(totalerr, norm), bottom=np.divide(self.fitBkg, norm), width=self.binWidths, align='center', label='{0} + {1}'.format(signalLabel, bkgLabel), capsize=3, color=self.signalColor, ec=self.signalColor, ecolor='gray')

        chi2text, = plt.plot([], [], ' ', label='Chi2/dof = {0:2.2f}/{1:2.0f}'.format(self.chi2, self.dof))
        puritytext, = plt.plot([], [], ' ', label='Purity = ${0:2.1f} \pm {1:2.1f}$%'.format(100 * self.purity, 100 * self.purityerr))

        ax = plt.gca()
        pmin, pmax = getBinRange(self.binEdges, self.purityMin, self.purityMax)
        ax.axvspan(self.binCenters[pmin] - self.binWidths[pmin] / 2.0, self.binCenters[pmax] - self.binWidths[pmax] / 2.0, color='black', alpha=0.1)

        return dataplot, bkgplot, sigplot, chi2text, puritytext

    def plotResiduals(self, ylim=[-8.9, 8.9]):
        plt.plot(self.binCenters, self.residuals, 'ko', alpha=0.65)
        plt.ylabel('(Fit- Data)/Dataerr', fontsize=26)
        plt.ylim(ylim)

    def plotBkgTemplateComponents(self):
        plt.figure(figsize=(8, 12))

        # contribution from the anti-isolated GJ MC
        plt.subplot(311)
        if self.usezeta:
            plt.bar(self.binCenters, np.subtract(self.antiisodatahist, np.multiply(self.fitzeta, self.antiisogjmchist)), width=self.binWidths, align='center',
                    label='Anti-isolated data - scaled GJ MC', capsize=0, color=self.bkgColor, ec=self.bkgColor)
            plt.bar(self.binCenters, np.multiply(self.fitzeta, self.antiisogjmchist), bottom=np.subtract(self.antiisodatahist, np.multiply(self.fitzeta, self.antiisogjmchist)), width=self.binWidths, align='center',
                    label='Scaled anti-iso GJ MC ($\zeta$={0:2.2})'.format(self.fitzeta), capsize=0, color=self.signalColor, ec=self.signalColor)
            plt.legend()
        elif self.useraa:
            plt.bar(self.binCenters, self.antiisomchist, width=self.binWidths, align='center',
                    label=r'$R_{AA}\times$ JJ MC + GJ MC', capsize=0, color=self.bkgColor, ec=self.bkgColor)
            plt.plot(self.binCenters, self.antiisojjmchist, 'k-', drawstyle='steps-mid', label='JJ MC')
            plt.legend()
        else:
            plt.bar(self.binCenters, self.antiisodatahist, width=self.binWidths, align='center', label='Anti-isolated data', capsize=0, color=self.bkgColor, ec=self.bkgColor)
            plt.legend()

        # background template correction
        plt.subplot(312)
        bkgweights, bkgweighterrs = self.getBkgCorrectionWeights()
        plt.errorbar(self.binCenters, bkgweights, yerr=bkgweighterrs, fmt='ko')
        plt.ylabel('Bkg template\ncorrection')
        if self.useraa:
            bkgtempcorrdesc = r'$R_\mathrm{AA} \times$ Iso JJ MC / ($R_\mathrm{AA} \times$ Anti-iso JJ MC + Anti-iso GJ MC)'
        else:
            bkgtempcorrdesc = 'Iso JJ MC / Anti-iso JJ MC'
        upperRightText(bkgtempcorrdesc)

        # double ratio
        plt.subplot(313)
        if self.drfit is not None:
            self.doubleratio[self.doubleratioerr == 0] = np.nan
            plotDoubleRatioAndFit(self.doubleratio, self.doubleratioerr, self.ssParams, self.drfit)
            if self.useraa:
                plt.ylabel(r'$\frac{\mathrm{(iso/antiiso)}_\mathrm{data}}{\mathrm{(iso/antiiso)}_\mathrm{MC}}$')
            else:
                plt.ylabel(r'$\frac{\mathrm{(iso/antiiso)}_\mathrm{data}}{\mathrm{(iso/antiiso)}_\mathrm{dijet MC}}$')
            plt.axvspan(self.doubleRatioFitRange[0], self.doubleRatioFitRange[1], color=self.bkgColor, alpha=0.2)
            plt.xlabel(self.xlabel)

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
