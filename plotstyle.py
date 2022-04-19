import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.rc('figure', figsize=(8, 8), titlesize=20)  # repeat this line in an attempt to get it to actually work
mpl.rc('figure', figsize=(8, 8), titlesize=20)
mpl.rc('savefig', bbox='tight', facecolor='white')
mpl.rc('legend', edgecolor='white', framealpha=0.0, fontsize=20)
mpl.rc('axes', labelsize=20, titlesize=20)
mpl.rc('xtick', labelsize=16, top=True, direction='in')
mpl.rc('ytick', labelsize=16, right=True, direction='in')


def upperRightText(text, xy=(0.96, 0.96), fontsize=24, **kwargs):
    plt.annotate(text, xy=xy, xycoords='axes fraction', ha='right', va='top', fontsize=fontsize, **kwargs)


def upperLeftText(text, xy=(0.04, 0.96), fontsize=24, **kwargs):
    plt.annotate(text, xy=xy, xycoords='axes fraction', ha='left', va='top', fontsize=fontsize, **kwargs)


def lowerRightText(text, xy=(0.96, 0.04), fontsize=24, **kwargs):
    plt.annotate(text, xy=xy, xycoords='axes fraction', ha='right', va='bottom', fontsize=fontsize, **kwargs)


def lowerLeftText(text, xy=(0.04, 0.04), fontsize=24, **kwargs):
    plt.annotate(text, xy=xy, xycoords='axes fraction', ha='left', va='bottom', fontsize=fontsize, **kwargs)


def setTicks(ax):
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', length=8, width=2, direction='in')
    ax.tick_params(axis='both', which='minor', length=4, width=1, direction='in')
