import matplotlib.pyplot as plt


def plotBfTfComp(dfs, ptrange, isoParams, ssParams):
    # fig = plt.figure(figsize=(16, 8))

    plt.subplot(121)
    bf = dfs.getBackgroundFit(ptrange, isoParams, ssParams)
    handles = bf.plotFit()
    plt.legend(handles=handles)
    plt.title('Background-only fit')

    plt.subplot(122)
    tf = dfs.getTemplateFit(ptrange, isoParams, ssParams)
    handles = tf.plotFit()
    plt.legend(handles=handles)
    plt.title('Template fit\n(no background template correction)')

#     fig.savefig('plots/bf-vs-tf-{0}-{1}-{2:2.0f}-{3:2.0f}.pdf'.format(dfs.system, ssParams.ssvar, *ptrange))


def getBfTfDiff(dfs, ptrange, isoParams, ssParams):
    bf = dfs.getBackgroundFit(ptrange, isoParams, ssParams)
    tf = dfs.getTemplateFit(ptrange, isoParams, ssParams)
    return tf.purity - bf.purity


# TO-DO: ADD BACK IN THE SIGNAL TEMPLATE SMEARING STUDY
