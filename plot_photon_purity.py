import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from utils import getCenters, getWidths


class PurityPlotParams:
    def __init__(self, **kwargs):
        self.color = 'C0'
        self.marker = 'o'
        self.boxFill = False
        self.boxAlpha = 0.1
        # self.boxHatch = '/'
        self.boxLinewidth = 1.25
        self.linestyle = '--'
        self.linewidth = 1.5

        for key in kwargs:
            setattr(self, key, kwargs[key])

    def getLineParams(self):
        p = {}
        p['color'] = self.color
        p['linestyle'] = self.linestyle
        p['linewidth'] = self.linewidth
        return p

    def getMarkerParams(self):
        p = {}
        p['color'] = self.color
        p['marker'] = self.marker
        p['linestyle'] = 'None'
        return p

    def getBoxParams(self):
        p = {}
        p['color'] = self.color
        p['fill'] = self.boxFill
        p['alpha'] = self.boxAlpha
        p['hatch'] = self.boxHatch
        p['linewidth'] = self.boxLinewidth
        return p


def plotErfFit(erfFit, xmin, xmax, label, plotParams):
    x = np.linspace(xmin, xmax, 200)
    erf, = plt.plot(x, np.multiply(map(erfFit.getFunctionFit(), x), 100), label=label, **plotParams)
    return erf


def plotSystematicBoxes(ptranges, photonPurities, ax, plotParams):
    purities = [100 * p.purity for p in photonPurities]
    systerrs = [100 * p.systerr for p in photonPurities]

    widths = getWidths(ptranges)
    heights = np.multiply(systerrs, 2)
    xs = [ledge for (ledge, hedge) in ptranges]
    ys = np.subtract(purities, systerrs)

    for x, y, width, height in zip(xs, ys, widths, heights):
        ax.add_patch(mpl.patches.Rectangle((x, y), width, height, **plotParams))


def plotPurityAndFit(ptranges, photonPurities, erfFit, system, plotParams):
    ptcenters = getCenters(ptranges)

    purities = [100 * p.purity for p in photonPurities]
    staterrs = [100 * p.staterr for p in photonPurities]

    ax = plt.gca()

    plt.errorbar(ptcenters, purities, yerr=staterrs, **plotParams.getMarkerParams())
    plotSystematicBoxes(ptranges, photonPurities, ax, plotParams.getBoxParams())
    erf = plotErfFit(erfFit, np.min(ptranges), np.max(ptranges), '{0} Erf fit'.format(system), plotParams.getLineParams())
    point, = plt.plot([], [], label=system, **plotParams.getMarkerParams())

    handles = [point, erf]
    return handles


def purityPlotLabels(handles):
    plt.xlabel('$p_{\mathrm{T}}$ (GeV/$c$)', x=1.0, ha='right')
    plt.ylabel('Purity (%)', y=1.0, ha='right')
    plt.ylim([0, 80])

    plt.legend(handles=handles, loc=2, fontsize=20)
