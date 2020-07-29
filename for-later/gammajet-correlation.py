import numpy as np
import math


npzfiles = []
centralities = [(0, 10), (10, 30), (30, 50), (50, 100)]
ptranges = [(15, 20), (20, 30)]
sigrange = (0.1, 0.3)
bkgrange = (0.4, 1.5)

# correlation distributions of interest, split first by centrality, then by pt
sigdeltaphi = {}
sigjetpt = {}
bkgdeltaphi = {}
bkgjetpt = {}

for centrality in centralities:
    sigdeltaphi[centrality] = {}
    sigjetpt[centrality] = {}
    bkgdeltaphi[centrality] = {}
    bkgjetpt[centrality] = {}

    for ptrange in ptranges:
        sigdeltaphi[centrality][ptrange] = []
        sigjetpt[centrality][ptrange] = []
        bkgdeltaphi[centrality][ptrange] = []
        bkgjetpt[centrality][ptrange] = []


def deltaphi(jet_phi, cluster_phi):
    # this is just a placeholder; figure out the better version later
    return jet_phi - cluster_phi


for npzfile in npzfiles:
    with np.load(npzfile, allow_pickle=True) as arrays:
        # get relevant branches
        acluster_pt = arrays['cluster_pt']
        acluster_eta = arrays['cluster_eta']
        acluster_phi = arrays['cluster_phi']
        acluster_e_cross = arrays['cluster_e_cross']
        acluster_e = arrays['cluster_e']
        acluster_e_max = arrays['cluster_e_max']
        acluster_ncell = arrays['cluster_ncell']
        acluster_distance_to_bad_channel = arrays['cluster_distance_to_bad_channel']
        acluster_lambda_square = arrays['cluster_lambda_square']
        acluster_iso_tpc_04 = arrays['cluster_iso_tpc_04']
        acluster_iso_tpc_04_ue = arrays['cluster_iso_tpc_04_ue']
        aue_estimate_tpc_const = arrays['ue_estimate_tpc_const']

        ajet_ak04tpc_pt_raw = arrays['jet_ak04tpc_pt_raw']
        ajet_ak04tpc_phi = arrays['jet_ak04tpc_phi']

        acentrality_v0m = arrays['centrality_v0m']

        # loop through events
        for ievent in range(len(acentrality_v0m)):
            centrality_v0m = acentrality_v0m[ievent]
            ue_estimate_tpc_const = aue_estimate_tpc_const[ievent]

            # if the centrality is not in the range of interest, skip this event
            if centrality_v0m < np.min(centralities) or centrality_v0m > np.max(centralities):
                continue

            clusterptmax = 0
            iclusterptmax = -1

            for icluster in range(len(acluster_pt[ievent])):
                # first, check if the cluster passes all the cluster cuts
                if abs(acluster_eta[ievent][icluster]) > 0.67:
                    continue
                if (acluster_e_cross[ievent][icluster] / acluster_e[ievent][icluster]) < 0.05:
                    continue
                if acluster_ncell[ievent][icluster] < 2:
                    continue
                if acluster_distance_to_bad_channel[ievent][icluster] < 1.0:
                    continue

                # check if the cluster passes the isolation cut
                cluster_iso_tpc_04_sub = acluster_iso_tpc_04[ievent][icluster] + acluster_iso_tpc_04_ue[ievent][icluster] - (ue_estimate_tpc_const * 0.4 * 0.4 * math.pi)
                if cluster_iso_tpc_04_sub > 1.5:
                    continue

                # check if the cluster is the highest-pt cluster yet found
                cluster_pt = acluster_pt[ievent][icluster]
                if cluster_pt > clusterptmax:
                    clusterptmax = cluster_pt
                    iclusterptmax = icluster

            # if the max-pt cluster is out of the ptrange, skip this event
            if clusterptmax < np.min(ptranges) or clusterptmax > np.max(ptranges):
                continue

            # find the centrality range of interest
            thiscentrality = None
            for centrality in centralities:
                if centrality_v0m > centrality[0] and centrality_v0m < centrality[1]:
                    thiscentrality = centrality

            # find the pT range of interest
            thisptrange = None
            for ptrange in ptranges:
                if clusterptmax > ptrange[0] and clusterptmax < ptrange[1]:
                    thisptrange = ptrange

            # figure out if this is signal or background
            cluster_lambda0 = acluster_lambda_square[ievent][iclusterptmax][0]
            if cluster_lambda0 > sigrange[0] and cluster_lambda0 < sigrange[1]:
                isSignal = True
            elif cluster_lambda0 > bkgrange[0] and cluster_lambda0 < bkgrange[1]:
                isSignal = False
            else:
                continue

            cluster_phi = acluster_phi[ievent][iclusterptmax]

            # loop through jets
            for ijet in range(len(ajet_ak04tpc_pt_raw[ievent])):
                # ignore jets with pT < 10
                jet_pt = ajet_ak04tpc_pt_raw[ievent][ijet]
                if jet_pt < 10:
                    continue

                # fill histograms
                jet_phi = ajet_ak04tpc_phi[ievent][ijet]
                if isSignal:
                    sigdeltaphi[thiscentrality][thisptrange].append(deltaphi(cluster_phi, jet_phi))
                    sigjetpt[thiscentrality][thisptrange].append(jet_pt)
                else:
                    bkgdeltaphi[thiscentrality][thisptrange].append(deltaphi(cluster_phi, jet_phi))
                    bkgjetpt[thiscentrality][thisptrange].append(jet_pt)
