import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.rc('figure', figsize=(8, 8), titlesize=20)
mpl.rc('savefig', bbox='tight', facecolor='white')
mpl.rc('legend', edgecolor='white', framealpha=0.0, fontsize=20)
mpl.rc('axes', labelsize=20, titlesize=20)
mpl.rc('xtick', labelsize=16, top=True, direction='in')
mpl.rc('ytick', labelsize=16, right=True, direction='in')


def upperRightText(text, fontsize=16):
    plt.annotate(text, xy=(0.98, 0.98), xycoords='axes fraction', ha='right', va='top', fontsize=fontsize)


def upperLeftText(text, fontsize=16):
    plt.annotate(text, xy=(0.02, 0.98), xycoords='axes fraction', ha='left', va='top', fontsize=fontsize)


def lowerRightText(text, fontsize=16):
    plt.annotate(text, xy=(0.98, 0.02), xycoords='axes fraction', ha='right', va='bottom', fontsize=fontsize)


def lowerLeftText(text, fontsize=16):
    plt.annotate(text, xy=(0.02, 0.02), xycoords='axes fraction', ha='left', va='bottom', fontsize=fontsize)


def setTicks(ax):
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', length=8, width=2, direction='in')
    ax.tick_params(axis='both', which='minor', length=4, width=1, direction='in')
