import matplotlib.pyplot as plt

from params import IsolationParams
from utils import getCenters, getUniformUncertainty, getXerrForPlot


def computeSidebandWindows(dfs, ptrange, baseIsoParams, ssParams, sidebands):
    purities = []
    chi2dofs = []

    for sideband in sidebands:
        isoParams = IsolationParams()
        isoParams.isovar = baseIsoParams.isovar
        isoParams.isocut = baseIsoParams.isocut
        isoParams.antiisocutlow = sideband[0]
        isoParams.antiisocuthigh = sideband[1]

        tf = dfs.getTemplateFit(ptrange, isoParams, ssParams)
        purities.append(tf.purity)
        chi2dofs.append(tf.chi2 / tf.dof)

    return purities, chi2dofs


def isSidebandInRange(sideband, isoParams):
    return sideband[0] >= isoParams.antiisocutlow and sideband[1] <= isoParams.antiisocuthigh


def getValuesInRange(values, sidebands, isoParams):
    rangeValues = [value for (value, sideband) in zip(values, sidebands) if isSidebandInRange(sideband, isoParams)]

    return rangeValues


def getSidebandSelectionUncertainty(dfs, ptrange, baseIsoParams, ssParams, sidebands):
    purities, _ = computeSidebandWindows(dfs, ptrange, baseIsoParams, ssParams, sidebands)

    return getUniformUncertainty(getValuesInRange(purities, sidebands, baseIsoParams))


def plotSidebandWindows(dfs, ptrange, baseIsoParams, ssParams, sidebands):
    purities, chi2dofs = computeSidebandWindows(dfs, ptrange, baseIsoParams, ssParams, sidebands)
    rangePurities = getValuesInRange(purities, sidebands, baseIsoParams)
    rangeChi2dofs = getValuesInRange(chi2dofs, sidebands, baseIsoParams)

    # fig = plt.figure(figsize=(16, 8))

    plt.subplot(121)
    ax = plt.gca()
    plt.errorbar(getCenters(sidebands), purities, xerr=getXerrForPlot(sidebands), fmt='b.')
    ax.axhspan(min(rangePurities), max(rangePurities), color='b', alpha=0.2)
    ax.axvspan(baseIsoParams.antiisocutlow, baseIsoParams.antiisocuthigh, color='green', alpha=0.2)
    plt.xlabel('Isolation energy (GeV)')
    plt.ylabel('Purities')

    plt.subplot(122)
    ax = plt.gca()
    plt.errorbar(getCenters(sidebands), chi2dofs, xerr=getXerrForPlot(sidebands), fmt='r.')
    ax.axhspan(min(rangeChi2dofs), max(rangeChi2dofs), color='r', alpha=0.2)
    ax.axvspan(baseIsoParams.antiisocutlow, baseIsoParams.antiisocuthigh, color='green', alpha=0.2)
    plt.xlabel('Isolation energy (GeV)')
    plt.ylabel('Chi2/dof')

    # fig.savefig('plots/sideband-windows-{0}-{1}-{2:2.0f}-{3:2.0f}.pdf'.format(dfs.system, ssParams.ssvar, ptrange[0], ptrange[1]))
