from __future__ import print_function

import csv
import datetime
import numpy as np
import os
import ROOT
import root_numpy as rnp
import struct
import sys
import time
import warnings
import yaml

from alice_emcal import get5x5Cells, getCrossCells, calculateShowerShapeFromCells
from alice_mc_weights import fixedWeightProductions, getFixedWeight
from alice_raa import getWeightWithRaa502charged
from alice_triggers import getINT7TriggerIds, getCentralTriggerIds, getSemiCentralTriggerIds, getEMCEGATriggerIds, isEventSelected


def calculateShowerShapes5x5(icluster, cluster_e, cellMaxId, cell_e, cell_cluster_index):
    # first, get the cell Ids of the surrounding 5x5 grid
    cells5x5 = get5x5Cells(cellMaxId)

    # this is not supposed to happen, so return a nonsensical number to see how there are
    # and also to not count them
    if cell_e[cellMaxId] < 0.1:
        return {
            '5x5contiguouscluster': -1.0,
            '5x5contiguous': -1.0,
            '5x5cluster': -1.0,
            '5x5all': -1.0
        }

    cells_ss5x5contiguouscluster = set()
    cells_ss5x5contiguous = set()
    cells_ss5x5cluster = set()
    cells_ss5x5all = set()

    # shower shape that requires contiguous cells in the 5x5 in the cluster:
    # recursively, starting with the max energy cell (cellMaxId):
    # 1. look for the crossing cells
    # 2. for each cell, check that it is not already in the set of cells to be used in the calculation
    # 3. if so, determine if it is in the 5x5
    # 4. if so, determine if it is in the cluster
    # 5. if so, add to the set of cells to be used in the calculation
    # 6. recurse
    cells_ss5x5contiguouscluster.add(cellMaxId)
    checkedCells = set()
    checkedCells.add(cellMaxId)

    def getContiguousClusterCells(cellId):
        crossCells = getCrossCells(cellId)
        for crossCellId in crossCells:
            if crossCellId in checkedCells:
                continue
            checkedCells.add(crossCellId)
            if crossCellId > 17663 or crossCellId < 0:
                continue
            if cell_e[crossCellId] < 0.1 or np.isnan(cell_e[crossCellId]):
                continue
            if crossCellId not in cells5x5:
                continue
            if cell_cluster_index[crossCellId] != icluster:
                continue
            cells_ss5x5contiguouscluster.add(crossCellId)
            getContiguousClusterCells(crossCellId)

    getContiguousClusterCells(cellMaxId)

    # shower shape that requires contiguous cells in the 5x5, not necessarily in the cluster:
    # recursively, starting with the max energy cell (cellMaxId):
    # 1. look for the crossing cells
    # 2. for each cell, check that it is not already in the set of cells to be used in the calculation
    # 3. if so, determine if it is in the 5x5
    # 4. if so, add it to the set of cells to be used in the calculation
    # 5. recurse
    # calculate
    cells_ss5x5contiguous.add(cellMaxId)
    checkedCells = set()
    checkedCells.add(cellMaxId)

    def getContiguousCells(cellId):
        crossCells = getCrossCells(cellId)
        for crossCellId in crossCells:
            if crossCellId in checkedCells:
                continue
            checkedCells.add(crossCellId)
            if crossCellId > 17663 or crossCellId < 0:
                continue
            if cell_e[crossCellId] < 0.1 or np.isnan(cell_e[crossCellId]):
                continue
            if crossCellId not in cells5x5:
                continue
            cells_ss5x5contiguous.add(crossCellId)
            getContiguousCells(crossCellId)

    getContiguousCells(cellMaxId)

    # shower shape that takes all cells in the 5x5 in the cluster, not necessarily contiguous
    # for each cell in the 5x5, add it to the set of cells to be used in the calculation if it is in the cluster

    # shower shape that takes all cells in the 5x5, not necessarily contiguous or in the cluster
    # for each cell in the 5x5, add it to the set of cells to be used in the calculation

    for cellId in cells5x5:
        if cellId > 17663 or cellId < 0:
            continue
        if cell_e[cellId] < 0.1 or np.isnan(cell_e[cellId]):
            continue
        cells_ss5x5all.add(cellId)
        if cell_cluster_index[cellId] == icluster:
            cells_ss5x5cluster.add(cellId)

    ss5x5contiguouscluster = calculateShowerShapeFromCells(cells_ss5x5contiguouscluster, cell_e, cluster_e)
    ss5x5contiguous = calculateShowerShapeFromCells(cells_ss5x5contiguous, cell_e, cluster_e)
    ss5x5cluster = calculateShowerShapeFromCells(cells_ss5x5cluster, cell_e, cluster_e)
    ss5x5all = calculateShowerShapeFromCells(cells_ss5x5all, cell_e, cluster_e)

    return {
        '5x5contiguouscluster': ss5x5contiguouscluster,
        '5x5contiguous': ss5x5contiguous,
        '5x5cluster': ss5x5cluster,
        '5x5all': ss5x5all
    }


def getNMCPhotons(cluster_nmc_truth, cluster_mc_truth_index, mc_truth_pdg_code):
    photonCount = 0
    for i in range(min(cluster_nmc_truth, len(cluster_mc_truth_index))):
        mcindex = cluster_mc_truth_index[i]
        if mcindex == 65535:
            continue
        if mc_truth_pdg_code[mcindex] == 22:
            photonCount = photonCount + 1

    return photonCount


# see the ntuplizer for how mc_truth_is_prompt_photon is determined
# it is only filled for photons and electrons (and positrons), so that needs to be checked first
# current definition of prompt photon cluster: any cluster that contains
# at least one photon/electron/positron that came from a prompt photon
def getIsPrompt(cluster_nmc_truth, cluster_mc_truth_index, mc_truth_pdg_code, mc_truth_is_prompt_photon):
    for i in range(min(cluster_nmc_truth, len(cluster_mc_truth_index))):
        mcindex = cluster_mc_truth_index[i]
        # if this is the dummy index, ignore it
        if mcindex == 65535:
            continue
        # if this is not a photon or electron, ignore it
        if mc_truth_pdg_code[mcindex] not in (11, -11, 22):
            continue
        # if this is from a prompt photon, then the cluster is considered a prompt photon cluster
        if mc_truth_is_prompt_photon[mcindex]:
            return True
    return False


def getTruthPtAndComponents(cluster_nmc_truth, cluster_mc_truth_index, mc_truth_pdg_code, mc_truth_is_prompt_photon, mc_truth_pt):
    truthpt = 0
    truthcomponents = 0
    for i in range(min(cluster_nmc_truth, len(cluster_mc_truth_index))):
        mcindex = cluster_mc_truth_index[i]
        # if this is the dummy index, ignore it
        if mcindex == 65535:
            continue
        # if this is not a photon or electron, ignore it
        if mc_truth_pdg_code[mcindex] not in (11, -11, 22):
            continue
        # if it is not from a prompt photon, ignore it
        if not mc_truth_is_prompt_photon[mcindex]:
            continue

        # account for the truth particles in this cluster: 1 = photon, 100 = electron, 10000 = positron
        # surely I could do 1/10/100, but juuuuuust in case
        if mc_truth_pdg_code[mcindex] == 22:
            truthcomponents += 1
        elif mc_truth_pdg_code[mcindex] == 11:
            truthcomponents += 100
        elif mc_truth_pdg_code[mcindex] == -11:
            truthcomponents += 10000

        # and then add the truth pt to the total
        truthpt += mc_truth_pt[mcindex]

    return truthpt, truthcomponents


def getParentPi0Pt(cluster_nmc_truth, cluster_mc_truth_index, mc_truth_first_parent_pdg_code, mc_truth_first_parent_pt):
    parentpi0pt = -1
    for i in range(min(cluster_nmc_truth, len(cluster_mc_truth_index))):
        mcindex = cluster_mc_truth_index[i]
        # if this is the dummy index, ignore it
        if mcindex == 65535:
            continue
        # if the parent is not a pi0, ignore it
        # not sure if pi0 can have negative pdg code, but check just in case
        if abs(mc_truth_first_parent_pdg_code[mcindex]) != 111:
            continue
        # for now, take the min pt of the pi0s we see
        # unclear how often we end up with more than one, but just in case,
        # we'll be consistent about which one we take
        if parentpi0pt == -1:
            parentpi0pt = mc_truth_first_parent_pt[mcindex]
        else:
            parentpi0pt = min(parentpi0pt, mc_truth_first_parent_pt[mcindex])

    return parentpi0pt


def main(ntuplefilenames, csvfilename):
    warnings.simplefilter('ignore', RuntimeWarning)
    start = time.time()

    with open(csvfilename, 'wb') as csvfile:
        fieldnames = ['cluster_pt', 'cluster_eta', 'cluster_phi', 'cluster_e_cross', 'cluster_e', 'cluster_e_max', 'cluster_ncell',
                      'cluster_tof', 'cluster_distance_to_bad_channel', 'cluster_Lambda', 'cluster_NN1',
                      'cluster_ecross_e', 'cluster_ecross_emax', 'cluster_nlocal_maxima', 'cluster_nmc_photon',
                      'cluster_mc_is_prompt_photon', 'cluster_mc_parentpi0pt', 'cluster_mc_truth_pt', 'cluster_mc_truth_components',
                      'cluster_iso_tpc_01', 'cluster_iso_tpc_02', 'cluster_iso_tpc_03', 'cluster_iso_tpc_04',
                      'cluster_iso_tpc_01_sub', 'cluster_iso_tpc_02_sub', 'cluster_iso_tpc_03_sub', 'cluster_iso_tpc_04_sub',
                      'cluster_frixione_tpc_04_02', 'cluster_frixione_tpc_04_05', 'cluster_frixione_tpc_04_10',
                      'cluster_5x5contiguouscluster', 'cluster_5x5contiguous', 'cluster_5x5cluster', 'cluster_5x5all',
                      'cluster_weight_with_raa',
                      'ue_estimate_tpc_const', 'weights', 'centrality_v0m', 'isINT7', 'isCentral', 'isSemiCentral', 'isEMCEGA']
        csvwriter = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        csvwriter.writeheader()
        clustercount = 0
        totalevents = 0

        ntuplefilenames.sort()
        # initial loop to count nevents per pthat bin, because that's what we actually need for weights
        pthat_nevents = {}
        for ntuplefilename in ntuplefilenames:
            if 'pthat' not in ntuplefilename:
                continue
            pthat = ntuplefilename.split('pthat')[1][0]
            if pthat not in pthat_nevents:
                pthat_nevents[pthat] = 0
            rootfile = ROOT.TFile.Open(ntuplefilename, 'READ')
            tree = rootfile.Get('AliAnalysisTaskNTGJ/_tree_event')
            nevents = tree.GetEntries()
            pthat_nevents[pthat] = pthat_nevents[pthat] + nevents
            rootfile.Close()

        for (ifile, ntuplefilename) in enumerate(ntuplefilenames):
            print('{0} Processing {1} ({2}/{3})'.format(datetime.datetime.now(), os.path.basename(ntuplefilename), ifile + 1, len(ntuplefilenames)))

            rootfile = ROOT.TFile.Open(ntuplefilename, 'READ')
            tree = rootfile.Get('AliAnalysisTaskNTGJ/_tree_event')
            nevents = tree.GetEntries()
            totalevents = totalevents + nevents

            # can't loop through these in the normal way
            acluster_lambda_square = rnp.tree2array(tree, branches='cluster_lambda_square')
            acluster_s_nphoton = rnp.tree2array(tree, branches='cluster_s_nphoton')
            acluster_mc_truth_index = rnp.tree2array(tree, branches='cluster_mc_truth_index')
            atrigger_mask = rnp.tree2array(tree, branches='trigger_mask')

            # compute MC weights or set to 1
            aeg_cross_section = rnp.tree2array(tree, branches='eg_cross_section')
            aeg_ntrial = rnp.tree2array(tree, branches='eg_ntrial')
            if 'pthat' in ntuplefilename:
                for production in fixedWeightProductions:
                    if production in ntuplefilename:
                        fixedWeights = True
                if fixedWeights:
                    aweights = np.full(nevents, getFixedWeight(ntuplefilename))
                else:
                    avg_eg_ntrial = np.mean(aeg_ntrial)
                    pthatnevents = pthat_nevents[ntuplefilename.split('pthat')[1][0]]
                    aweights = np.divide(aeg_cross_section, avg_eg_ntrial * pthatnevents)
            else:
                aweights = np.ones(nevents)

            # track previous run number so that we don't have to recalculate the trigger IDs if it didn't change
            previous_run_number = -1

            # loop through events
            for ievent in range(nevents):
                tree.GetEntry(ievent)

                ncluster = getattr(tree, 'ncluster')
                cluster_pt = getattr(tree, 'cluster_pt')
                cluster_eta = getattr(tree, 'cluster_eta')
                cluster_phi = getattr(tree, 'cluster_phi')
                cluster_e_cross = getattr(tree, 'cluster_e_cross')
                cluster_e = getattr(tree, 'cluster_e')
                cluster_e_max = getattr(tree, 'cluster_e_max')
                cluster_ncell = getattr(tree, 'cluster_ncell')
                cluster_tof = getattr(tree, 'cluster_tof')
                cluster_nlocal_maxima = getattr(tree, 'cluster_nlocal_maxima')
                cluster_distance_to_bad_channel = getattr(tree, 'cluster_distance_to_bad_channel')
                cluster_cell_id_max = getattr(tree, 'cluster_cell_id_max')
                cluster_iso_tpc_01 = getattr(tree, 'cluster_iso_tpc_01')
                cluster_iso_tpc_02 = getattr(tree, 'cluster_iso_tpc_02')
                cluster_iso_tpc_03 = getattr(tree, 'cluster_iso_tpc_03')
                cluster_iso_tpc_04 = getattr(tree, 'cluster_iso_tpc_04')
                cluster_iso_tpc_01_ue = getattr(tree, 'cluster_iso_tpc_01_ue')
                cluster_iso_tpc_02_ue = getattr(tree, 'cluster_iso_tpc_02_ue')
                cluster_iso_tpc_03_ue = getattr(tree, 'cluster_iso_tpc_03_ue')
                cluster_iso_tpc_04_ue = getattr(tree, 'cluster_iso_tpc_04_ue')
                cluster_frixione_tpc_04_02 = getattr(tree, 'cluster_frixione_tpc_04_02')
                cluster_frixione_tpc_04_05 = getattr(tree, 'cluster_frixione_tpc_04_05')
                cluster_frixione_tpc_04_10 = getattr(tree, 'cluster_frixione_tpc_04_10')

                cluster_nmc_truth = getattr(tree, 'cluster_nmc_truth')
                mc_truth_pdg_code = getattr(tree, 'mc_truth_pdg_code')
                mc_truth_pt = getattr(tree, 'mc_truth_pt')
                mc_truth_first_parent_pdg_code = getattr(tree, 'mc_truth_first_parent_pdg_code')
                mc_truth_first_parent_pt = getattr(tree, 'mc_truth_first_parent_pt')
                try:
                    mc_truth_is_prompt_photon = getattr(tree, 'mc_truth_is_prompt_photon')
                except AttributeError:
                    mc_truth_is_prompt_photon = None

                centrality_v0m = getattr(tree, 'centrality_v0m')
                ue_estimate_tpc_const = getattr(tree, 'ue_estimate_tpc_const')
                run_number = getattr(tree, "run_number")

                if run_number != previous_run_number:
                    kINT7TriggerIds = getINT7TriggerIds(run_number)
                    kCentralTriggerIds = getCentralTriggerIds(run_number)
                    kSemiCentralTriggerIds = getSemiCentralTriggerIds(run_number)
                    kEMCEGATriggerIds = getEMCEGATriggerIds(run_number)
                    previous_run_number = run_number

                cell_e = getattr(tree, 'cell_e')
                cell_cluster_index = getattr(tree, 'cell_cluster_index')

                # event-based values
                weight = aweights[ievent]
                # combine the trigger masks into one number
                triggerMask = atrigger_mask[ievent][0] + (int(atrigger_mask[ievent][1]) << 50)

                isINT7 = isEventSelected(kINT7TriggerIds, triggerMask)
                isCentral = isEventSelected(kCentralTriggerIds, triggerMask)
                isSemiCentral = isEventSelected(kSemiCentralTriggerIds, triggerMask)
                isEMCEGA = isEventSelected(kEMCEGATriggerIds, triggerMask)

                for icluster in range(ncluster):
                    # cluster-based values
                    pt = cluster_pt[icluster]
                    eta = cluster_eta[icluster]
                    phi = cluster_phi[icluster]
                    e_cross = cluster_e_cross[icluster]
                    e = cluster_e[icluster]
                    e_max = cluster_e_max[icluster]
                    ncell = cluster_ncell[icluster]
                    tof = cluster_tof[icluster]
                    nlm = struct.unpack('1b', cluster_nlocal_maxima[icluster])[0]
                    dbc = cluster_distance_to_bad_channel[icluster]
                    cell_id_max = cluster_cell_id_max[icluster]
                    iso_tpc_01 = cluster_iso_tpc_01[icluster]
                    iso_tpc_02 = cluster_iso_tpc_02[icluster]
                    iso_tpc_03 = cluster_iso_tpc_03[icluster]
                    iso_tpc_04 = cluster_iso_tpc_04[icluster]
                    iso_tpc_01_ue = cluster_iso_tpc_01_ue[icluster]
                    iso_tpc_02_ue = cluster_iso_tpc_02_ue[icluster]
                    iso_tpc_03_ue = cluster_iso_tpc_03_ue[icluster]
                    iso_tpc_04_ue = cluster_iso_tpc_04_ue[icluster]
                    iso_tpc_01_sub = iso_tpc_01 + iso_tpc_01_ue - (ue_estimate_tpc_const * 0.1 * 0.1 * np.pi)
                    iso_tpc_02_sub = iso_tpc_02 + iso_tpc_02_ue - (ue_estimate_tpc_const * 0.2 * 0.2 * np.pi)
                    iso_tpc_03_sub = iso_tpc_03 + iso_tpc_03_ue - (ue_estimate_tpc_const * 0.3 * 0.3 * np.pi)
                    iso_tpc_04_sub = iso_tpc_04 + iso_tpc_04_ue - (ue_estimate_tpc_const * 0.4 * 0.4 * np.pi)
                    frixione_tpc_04_02 = cluster_frixione_tpc_04_02[icluster]
                    frixione_tpc_04_05 = cluster_frixione_tpc_04_05[icluster]
                    frixione_tpc_04_10 = cluster_frixione_tpc_04_10[icluster]
                    lambda0 = acluster_lambda_square[ievent][icluster][0]
                    nn1 = acluster_s_nphoton[ievent][icluster][1]
                    nmc_photon = getNMCPhotons(cluster_nmc_truth[icluster], acluster_mc_truth_index[ievent][icluster], mc_truth_pdg_code)
                    if mc_truth_is_prompt_photon is None:
                        isprompt = False
                    else:
                        isprompt = getIsPrompt(cluster_nmc_truth[icluster], acluster_mc_truth_index[ievent][icluster],
                                               mc_truth_pdg_code, mc_truth_is_prompt_photon)

                    if isprompt:
                        truth_pt, mc_truth_components = getTruthPtAndComponents(cluster_nmc_truth[icluster],
                                                                                acluster_mc_truth_index[ievent][icluster],
                                                                                mc_truth_pdg_code,
                                                                                mc_truth_is_prompt_photon,
                                                                                mc_truth_pt)
                        parentpi0pt = -1
                    else:
                        truth_pt = -1
                        mc_truth_components = -1
                        parentpi0pt = getParentPi0Pt(cluster_nmc_truth[icluster], acluster_mc_truth_index[ievent][icluster],
                                                     mc_truth_first_parent_pdg_code, mc_truth_first_parent_pt)

                    if pt < 15:
                        continue

                    # not sure how this happens with the pt cut, but avoiding ZeroDivisionErrors is good
                    if e == 0 or e_max == 0:
                        continue

                    # this takes time, so put it after the pT cut
                    weightwithraa = getWeightWithRaa502charged(weight, parentpi0pt, centrality_v0m)
                    showershapes5x5 = calculateShowerShapes5x5(icluster, e, cell_id_max, cell_e, cell_cluster_index)

                    row = {}
                    row['cluster_pt'] = pt
                    row['cluster_eta'] = eta
                    row['cluster_phi'] = phi
                    row['cluster_e_cross'] = e_cross
                    row['cluster_e'] = e
                    row['cluster_e_max'] = e_max
                    row['cluster_ncell'] = ncell
                    row['cluster_tof'] = tof
                    row['cluster_nlocal_maxima'] = nlm
                    row['cluster_distance_to_bad_channel'] = dbc
                    row['ue_estimate_tpc_const'] = ue_estimate_tpc_const
                    row['cluster_iso_tpc_01'] = iso_tpc_01
                    row['cluster_iso_tpc_02'] = iso_tpc_02
                    row['cluster_iso_tpc_03'] = iso_tpc_03
                    row['cluster_iso_tpc_04'] = iso_tpc_04
                    row['cluster_iso_tpc_01_sub'] = iso_tpc_01_sub
                    row['cluster_iso_tpc_02_sub'] = iso_tpc_02_sub
                    row['cluster_iso_tpc_03_sub'] = iso_tpc_03_sub
                    row['cluster_iso_tpc_04_sub'] = iso_tpc_04_sub
                    row['cluster_frixione_tpc_04_02'] = frixione_tpc_04_02
                    row['cluster_frixione_tpc_04_05'] = frixione_tpc_04_05
                    row['cluster_frixione_tpc_04_10'] = frixione_tpc_04_10
                    row['cluster_Lambda'] = lambda0
                    row['cluster_NN1'] = nn1
                    row['cluster_nmc_photon'] = nmc_photon
                    row['cluster_mc_is_prompt_photon'] = isprompt
                    row['cluster_mc_parentpi0pt'] = parentpi0pt
                    row['cluster_mc_truth_pt'] = truth_pt
                    row['cluster_mc_truth_components'] = mc_truth_components
                    row['cluster_ecross_e'] = e_cross / e
                    row['cluster_ecross_emax'] = e_cross / e_max
                    row['cluster_weight_with_raa'] = weightwithraa
                    row['weights'] = weight
                    row['centrality_v0m'] = centrality_v0m
                    row['isINT7'] = isINT7
                    row['isCentral'] = isCentral
                    row['isSemiCentral'] = isSemiCentral
                    row['isEMCEGA'] = isEMCEGA

                    row['cluster_5x5contiguouscluster'] = showershapes5x5['5x5contiguouscluster']
                    row['cluster_5x5contiguous'] = showershapes5x5['5x5contiguous']
                    row['cluster_5x5cluster'] = showershapes5x5['5x5cluster']
                    row['cluster_5x5all'] = showershapes5x5['5x5all']

                    csvwriter.writerow(row)
                    clustercount = clustercount + 1

            rootfile.Close()

    end = time.time()
    duration = end - start
    if duration < 120:
        timetext = '{0:0.0f} seconds'.format(duration)
    elif duration < 3600:
        timetext = '{0:0.0f} minutes'.format(duration / 60.0)
    else:
        timetext = '{0:0.0f} hours {1:0.0f} minutes'.format(duration / 3600, (duration % 3600) / 60.0)

    csvsize = float(os.path.getsize(csvfilename))
    if csvsize < 1024:
        sizetext = '{0} B'.format(csvsize)
    elif csvsize < 1024 * 1024:
        sizetext = '{0} kB'.format(csvsize / 1024)
    elif csvsize < 1024 * 1024 * 1024:
        sizetext = '{0:0.1f} MB'.format(csvsize / (1024 * 1024))
    else:
        sizetext = '{0:0.1f} GB'.format(csvsize / (1024 * 1024 * 1024))

    print('Took {0} to produce {1} ({2}, {3:0.0f} clusters, {4:0.0f} events)'.format(timetext, os.path.basename(csvfilename), sizetext, clustercount, totalevents))


if __name__ == '__main__':
    configfilename = sys.argv[1]
    filetype = sys.argv[2]

    if filetype not in ('data', 'gjmc', 'jjmc'):
        print('File type {0} not recognized; must be one of data, gjmc, jjmc'.format(filetype))
        exit()

    with open(configfilename) as configfile:
        config = yaml.safe_load(configfile)

    ntuplefilenames = config['filelists']['ntuples'][filetype]
    csvfilename = config['filelists']['clustercsvs'][filetype]
    main(ntuplefilenames, csvfilename)
