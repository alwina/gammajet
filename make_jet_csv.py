from __future__ import print_function

import csv
import datetime
import numpy as np
import os
import ROOT
import sys
import time
import warnings
import yaml

from utils import tBranchToArray, getTimeText, getSizeText


# for now, this works only for triggered data
# with exactly the cluster quality selections we have in the config file
def main(config):
    warnings.simplefilter('ignore', RuntimeWarning)
    start = time.time()

    ntuplefilenames = config['filelists']['ntuples']['data']
    csvfilename = config['filelists']['jetcsvs']['data']

    clustercuts = config['clustercuts']['all']
    clustercuts.update(config['clustercuts']['data'])

    jetcuts = config['jetcuts']

    with open(csvfilename, 'wb') as csvfile:
        fieldnames = ['cluster_pt', 'cluster_eta', 'cluster_phi', 'cluster_5x5all', 'cluster_iso_tpc_02_sub',
                      'centrality_v0m', 'ue_estimate_tpc_const', 'icluster',
                      'jet_ak02tpc_pt_raw', 'jet_ak02tpc_eta', 'jet_ak02tpc_phi',
                      'jet_ak02tpc_area', 'jet_ak02tpc_multiplicity_raw']
        csvwriter = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        csvwriter.writeheader()
        clustercount = 0
        jetcount = 0
        totalevents = 0

        ntuplefilenames.sort()
        for (ifile, ntuplefilename) in enumerate(ntuplefilenames):
            print('{0} Processing {1} ({2}/{3})'.format(datetime.datetime.now(), os.path.basename(ntuplefilename), ifile + 1, len(ntuplefilenames)))

            rootfile = ROOT.TFile.Open(ntuplefilename, 'READ')
            tree = rootfile.Get('AliAnalysisTaskNTGJ/_tree_event')
            nevents = tree.GetEntries()
            totalevents = totalevents + nevents

            auxfile = ROOT.TFile.Open(ntuplefilename.replace('.root', '_AUX.root'))
            auxtree = auxfile.Get('ntupleaux')

            # loop through events
            for ievent in range(nevents):
                tree.GetEntry(ievent)
                auxtree.GetEntry(ievent)

                ncluster = getattr(tree, 'ncluster')
                cluster_pt = getattr(tree, 'cluster_pt')
                cluster_eta = getattr(tree, 'cluster_eta')
                cluster_phi = getattr(tree, 'cluster_phi')
                cluster_e_cross = getattr(tree, 'cluster_e_cross')
                cluster_e_max = getattr(tree, 'cluster_e_max')
                cluster_ncell = getattr(tree, 'cluster_ncell')
                cluster_tof = getattr(tree, 'cluster_tof')
                cluster_distance_to_bad_channel = getattr(tree, 'cluster_distance_to_bad_channel')
                cluster_iso_tpc_02 = getattr(tree, 'cluster_iso_tpc_02')
                cluster_iso_tpc_02_ue = getattr(tree, 'cluster_iso_tpc_02_ue')
                cluster_lambda_square = tBranchToArray(getattr(tree, 'cluster_lambda_square'), 'F', (ncluster, 2))
                cluster_5x5all = getattr(auxtree, 'cluster_5x5all')

                njet = getattr(auxtree, 'njet_ak02tpc')
                jet_pt = getattr(auxtree, 'jet_ak02tpc_pt_raw')
                jet_eta = getattr(auxtree, 'jet_ak02tpc_eta')
                jet_phi = getattr(auxtree, 'jet_ak02tpc_phi')
                jet_area = getattr(auxtree, 'jet_ak02tpc_area')
                jet_multiplicity = tBranchToArray(getattr(auxtree, 'jet_ak02tpc_multiplicity_raw'), 's', njet)

                centrality_v0m = getattr(tree, 'centrality_v0m')
                ue_estimate_tpc_const = getattr(tree, 'ue_estimate_tpc_const')

                for icluster in range(ncluster):
                    # cluster-based values
                    pt = cluster_pt[icluster]
                    eta = cluster_eta[icluster]
                    phi = cluster_phi[icluster]
                    e_cross = cluster_e_cross[icluster]
                    e_max = cluster_e_max[icluster]
                    ncell = cluster_ncell[icluster]
                    tof = cluster_tof[icluster]
                    dbc = cluster_distance_to_bad_channel[icluster]
                    iso_tpc_02 = cluster_iso_tpc_02[icluster]
                    iso_tpc_02_ue = cluster_iso_tpc_02_ue[icluster]
                    iso_tpc_02_sub = iso_tpc_02 + iso_tpc_02_ue - (ue_estimate_tpc_const * 0.2 * 0.2 * np.pi)
                    lambda0 = cluster_lambda_square[icluster][0]
                    c5x5all = cluster_5x5all[icluster]

                    # make cluster cuts
                    if not(pt > clustercuts['cluster_pt']['min']):
                        continue
                    if not(pt < clustercuts['cluster_pt']['max']):
                        continue
                    if not(abs(eta) < clustercuts['cluster_eta']['absmax']):
                        continue
                    if not(e_cross / e_max > clustercuts['cluster_ecross_emax']['min']):
                        continue
                    if not(ncell >= clustercuts['cluster_ncell']['incmin']):
                        continue
                    if not(dbc >= clustercuts['cluster_distance_to_bad_channel']['incmin']):
                        continue
                    if not(lambda0 > clustercuts['cluster_Lambda']['min']):
                        continue
                    if not(abs(tof) < clustercuts['cluster_tof']['absmax']):
                        continue
                    # got lazy here...should probably read these from config too
                    if not(iso_tpc_02_sub < 1.5):
                        continue
                    if not((c5x5all > 0.1 and c5x5all < 0.3) or (c5x5all > 0.6 and c5x5all < 1.5)):
                        continue

                    # now loop through jets
                    for ijet in range(njet):
                        jpt = jet_pt[ijet]
                        jeta = jet_eta[ijet]
                        jphi = jet_phi[ijet]
                        jarea = jet_area[ijet]
                        jmult = jet_multiplicity[ijet]

                        # make jet cuts
                        if not(jpt > jetcuts['jet_pt_raw']['min']):
                            continue
                        if not(jpt < jetcuts['jet_pt_raw']['max']):
                            continue
                        if not(abs(jeta) < jetcuts['jet_eta']['max']):
                            continue

                        # write out to the CSV
                        row = {}
                        row['cluster_pt'] = pt
                        row['cluster_eta'] = eta
                        row['cluster_phi'] = phi
                        row['cluster_iso_tpc_02_sub'] = iso_tpc_02_sub
                        row['cluster_5x5all'] = c5x5all
                        row['ue_estimate_tpc_const'] = ue_estimate_tpc_const
                        row['centrality_v0m'] = centrality_v0m
                        row['icluster'] = clustercount
                        row['jet_ak02tpc_pt_raw'] = jpt
                        row['jet_ak02tpc_eta'] = jeta
                        row['jet_ak02tpc_phi'] = jphi
                        row['jet_ak02tpc_area'] = jarea
                        row['jet_ak02tpc_multiplicity_raw'] = jmult

                        csvwriter.writerow(row)
                        jetcount = jetcount + 1

                    clustercount = clustercount + 1

            auxfile.Close()
            rootfile.Close()

    end = time.time()
    timetext = getTimeText(end - start)

    csvsize = float(os.path.getsize(csvfilename))
    sizetext = getSizeText(csvsize)

    print('Took {0} to produce {1} ({2}, {3:.0f} jets, {4:.0f} clusters, {5:.0f} events)'.format(timetext, os.path.basename(csvfilename), sizetext, jetcount, clustercount, totalevents))


if __name__ == '__main__':
    configfilename = sys.argv[1]

    with open(configfilename) as runconfigfile:
        runconfig = yaml.safe_load(runconfigfile)
        
    with open(runconfig['systemconfig']) as sysconfigfile:
        sysconfig = yaml.safe_load(sysconfigfile)
        
    fullconfig = {}
    fullconfig.update(sysconfig)
    fullconfig.update(runconfig)

    main(fullconfig)
