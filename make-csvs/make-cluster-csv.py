from __future__ import print_function

import csv
import datetime
import math
import numpy as np
import os
import ROOT
import root_numpy as rnp
import sys
import time
import warnings


def getSuperModule(cellId):
    if cellId < 11520:
        return math.floor(cellId / 1152)
    elif cellId < 12288:
        return 10 + math.floor((cellId - 11520) / 384)
    elif cellId < 16896:
        return 12 + math.floor((cellId - 12288) / 768)
    else:
        return 18 + (cellId - 16896) / 384


def getNphi(sm):
    if sm < 10:
        return 24
    elif sm < 12:
        return 8
    elif sm < 18:
        return 24
    else:
        return 8


def getIetaIphi(cellId, sm, nphi):
    if sm < 10:
        n0 = sm * 1152
    elif sm < 12:
        n0 = 11520 + (sm - 10) * 384
    elif sm < 18:
        n0 = 12288 + (sm - 12) * 768
    else:
        n0 = 16896 + (sm - 18) * 384

    n1 = cellId - n0

    ieta = 2 * math.floor(n1 / (2 * nphi)) + 1 - (n1 % 2)
    iphi = math.floor(n1 / 2) % nphi

    return ieta, iphi


def get5x5Cells(cellId):
    sm = getSuperModule(cellId)
    nphi = getNphi(sm)
    cells = []

    for dn in range(0, 10, 2):
        cells.append(cellId - 2 * nphi - 4 + dn)

        if cellId % 2 == 0:
            cells.append(cellId - 3 + dn)
        else:
            cells.append(cellId - 2 * nphi - 5 + dn)

        cells.append(cellId - 4 + dn)

        if cellId % 2 == 0:
            cells.append(cellId + 2 * nphi - 3 + dn)
        else:
            cells.append(cellId - 5 + dn)

        cells.append(cellId + 2 * nphi - 4 + dn)

    return cells


def getCrossCells(cellId):
    sm = getSuperModule(cellId)
    nphi = getNphi(sm)
    cells = []

    if cellId % 2 == 0:
        cells.append(cellId + 1)
    else:
        cells.append(cellId - 2 * nphi - 1)

    cells.append(cellId - 2)
    cells.append(cellId + 2)

    if cellId % 2 == 0:
        cells.append(cellId + 2 * nphi + 1)
    else:
        cells.append(cellId - 1)

    return cells


# port from http://alidoc.cern.ch/AliRoot/v5-09-11/_ali_e_m_c_a_l_rec_point_8cxx_source.html#l01002
def calculateShowerShapeFromCells(cellIds, cell_e, cluster_e):
    wtot = 0.0
    x = 0.0
    z = 0.0
    dxx = 0.0
    dzz = 0.0
    dxz = 0.0

    for cellId in cellIds:
        sm = getSuperModule(cellId)
        nphi = getNphi(sm)
        ieta, iphi = getIetaIphi(cellId, sm, nphi)

        # this handles shared clusters but also doesn't do anything for
        # clusters in a single supermodule, so don't bother checking
        # specifically for shared clusters
        if (sm % 2):
            ieta = ieta + 48

        w = max(0, 4.5 + math.log(cell_e[cellId] / cluster_e))
        dxx = dxx + w * ieta * ieta
        x = x + w * ieta
        dzz = dzz + w * iphi * iphi
        z = z + w * iphi
        dxz = dxz + w * ieta * iphi
        wtot = wtot + w

    if wtot > 0:
        x = x / wtot
        z = z / wtot
        dxx = dxx / wtot - x * x
        dzz = dzz / wtot - z * z
        dxz = dxz / wtot - x * z

    return 0.5 * (dxx + dzz) + math.sqrt(0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz)


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


def main(ntuplefilenames, csvfilename):
    warnings.simplefilter('ignore', RuntimeWarning)
    start = time.time()

    with open(csvfilename, 'wb') as csvfile:
        fieldnames = ['cluster_pt', 'cluster_eta', 'cluster_phi', 'cluster_e_cross', 'cluster_e', 'cluster_e_max', 'cluster_ncell',
                      'cluster_tof', 'cluster_distance_to_bad_channel', 'cluster_Lambda', 'cluster_NN1',
                      'cluster_ecross_e', 'cluster_ecross_emax', 'cluster_nlocal_maxima', 'cluster_nmc_photon',
                      'cluster_mc_is_prompt_photon',
                      'cluster_iso_tpc_01', 'cluster_iso_tpc_02', 'cluster_iso_tpc_03', 'cluster_iso_tpc_04',
                      'cluster_iso_tpc_01_sub', 'cluster_iso_tpc_02_sub', 'cluster_iso_tpc_03_sub', 'cluster_iso_tpc_04_sub',
                      'cluster_frixione_tpc_04_02', 'cluster_frixione_tpc_04_05', 'cluster_frixione_tpc_04_10',
                      'ue_estimate_tpc_const', 'weights', 'centrality_v0m']
        csvwriter = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        csvwriter.writeheader()
        clustercount = 0
        totalevents = 0

        ntuplefilenames.sort()
        for (ifile, ntuplefilename) in enumerate(ntuplefilenames):
            print('{0} Processing {1} ({2}/{3})'.format(datetime.datetime.now(), os.path.basename(ntuplefilename), ifile + 1, len(ntuplefilenames)))

            rootfile = ROOT.TFile.Open(ntuplefilename, 'READ')
            tree = rootfile.Get('AliAnalysisTaskNTGJ/_tree_event')
            nevents = tree.GetEntries()
            totalevents = totalevents + nevents

            ancluster = rnp.tree2array(tree, branches='ncluster')
            acluster_pt = rnp.tree2array(tree, branches='cluster_pt')
            acluster_eta = rnp.tree2array(tree, branches='cluster_eta')
            acluster_phi = rnp.tree2array(tree, branches='cluster_phi')
            acluster_e_cross = rnp.tree2array(tree, branches='cluster_e_cross')
            acluster_e = rnp.tree2array(tree, branches='cluster_e')
            acluster_e_max = rnp.tree2array(tree, branches='cluster_e_max')
            acluster_ncell = rnp.tree2array(tree, branches='cluster_ncell')
            acluster_tof = rnp.tree2array(tree, branches='cluster_tof')
            acluster_nlocal_maxima = rnp.tree2array(tree, branches='cluster_nlocal_maxima')
            acluster_distance_to_bad_channel = rnp.tree2array(tree, branches='cluster_distance_to_bad_channel')
            acluster_iso_tpc_01 = rnp.tree2array(tree, branches='cluster_iso_tpc_01')
            acluster_iso_tpc_02 = rnp.tree2array(tree, branches='cluster_iso_tpc_02')
            acluster_iso_tpc_03 = rnp.tree2array(tree, branches='cluster_iso_tpc_03')
            acluster_iso_tpc_04 = rnp.tree2array(tree, branches='cluster_iso_tpc_04')
            acluster_iso_tpc_01_ue = rnp.tree2array(tree, branches='cluster_iso_tpc_01_ue')
            acluster_iso_tpc_02_ue = rnp.tree2array(tree, branches='cluster_iso_tpc_02_ue')
            acluster_iso_tpc_03_ue = rnp.tree2array(tree, branches='cluster_iso_tpc_03_ue')
            acluster_iso_tpc_04_ue = rnp.tree2array(tree, branches='cluster_iso_tpc_04_ue')
            acluster_frixione_tpc_04_02 = rnp.tree2array(tree, branches='cluster_frixione_tpc_04_02')
            acluster_frixione_tpc_04_05 = rnp.tree2array(tree, branches='cluster_frixione_tpc_04_05')
            acluster_frixione_tpc_04_10 = rnp.tree2array(tree, branches='cluster_frixione_tpc_04_10')
            acluster_lambda_square = rnp.tree2array(tree, branches='cluster_lambda_square')
            acluster_s_nphoton = rnp.tree2array(tree, branches='cluster_s_nphoton')

            acluster_nmc_truth = rnp.tree2array(tree, branches='cluster_nmc_truth')
            acluster_mc_truth_index = rnp.tree2array(tree, branches='cluster_mc_truth_index')
            amc_truth_pdg_code = rnp.tree2array(tree, branches='mc_truth_pdg_code')
            amc_truth_is_prompt_photon = rnp.tree2array(tree, branches='mc_truth_is_prompt_photon')

            acentrality_v0m = rnp.tree2array(tree, branches='centrality_v0m')
            aeg_cross_section = rnp.tree2array(tree, branches='eg_cross_section')
            aeg_ntrial = rnp.tree2array(tree, branches='eg_ntrial')
            aue_estimate_tpc_const = rnp.tree2array(tree, branches='ue_estimate_tpc_const')

            # compute MC weights or set to 1
            if 'pthat' in ntuplefilename:
                avg_eg_ntrial = np.mean(aeg_ntrial)
                nevents = len(aeg_ntrial)
                aweights = np.divide(aeg_cross_section, avg_eg_ntrial * nevents)
            else:
                aweights = np.ones_like(acentrality_v0m)

            # loop through events
            for ievent in range(nevents):
                # event-based values
                weight = aweights[ievent]
                centrality = acentrality_v0m[ievent]
                ue_estimate_tpc_const = aue_estimate_tpc_const[ievent]

                for icluster in range(ancluster[ievent]):
                    # cluster-based values
                    pt = acluster_pt[ievent][icluster]
                    eta = acluster_eta[ievent][icluster]
                    phi = acluster_phi[ievent][icluster]
                    e_cross = acluster_e_cross[ievent][icluster]
                    e = acluster_e[ievent][icluster]
                    e_max = acluster_e_max[ievent][icluster]
                    ncell = acluster_ncell[ievent][icluster]
                    tof = acluster_tof[ievent][icluster]
                    nlm = acluster_nlocal_maxima[ievent][icluster]
                    dbc = acluster_distance_to_bad_channel[ievent][icluster]
                    iso_tpc_01 = acluster_iso_tpc_01[ievent][icluster]
                    iso_tpc_02 = acluster_iso_tpc_02[ievent][icluster]
                    iso_tpc_03 = acluster_iso_tpc_03[ievent][icluster]
                    iso_tpc_04 = acluster_iso_tpc_04[ievent][icluster]
                    iso_tpc_01_ue = acluster_iso_tpc_01_ue[ievent][icluster]
                    iso_tpc_02_ue = acluster_iso_tpc_02_ue[ievent][icluster]
                    iso_tpc_03_ue = acluster_iso_tpc_03_ue[ievent][icluster]
                    iso_tpc_04_ue = acluster_iso_tpc_04_ue[ievent][icluster]
                    iso_tpc_01_sub = iso_tpc_01 + iso_tpc_01_ue - (ue_estimate_tpc_const * 0.1 * 0.1 * np.pi)
                    iso_tpc_02_sub = iso_tpc_02 + iso_tpc_02_ue - (ue_estimate_tpc_const * 0.2 * 0.2 * np.pi)
                    iso_tpc_03_sub = iso_tpc_03 + iso_tpc_03_ue - (ue_estimate_tpc_const * 0.3 * 0.3 * np.pi)
                    iso_tpc_04_sub = iso_tpc_04 + iso_tpc_04_ue - (ue_estimate_tpc_const * 0.4 * 0.4 * np.pi)
                    frixione_tpc_04_02 = acluster_frixione_tpc_04_02[ievent][icluster]
                    frixione_tpc_04_05 = acluster_frixione_tpc_04_05[ievent][icluster]
                    frixione_tpc_04_10 = acluster_frixione_tpc_04_10[ievent][icluster]
                    lambda0 = acluster_lambda_square[ievent][icluster][0]
                    nn1 = acluster_s_nphoton[ievent][icluster][1]
                    nmc_photon = getNMCPhotons(acluster_nmc_truth[ievent][icluster], acluster_mc_truth_index[ievent][icluster], amc_truth_pdg_code[ievent])
                    isprompt = getIsPrompt(acluster_nmc_truth[ievent][icluster], acluster_mc_truth_index[ievent][icluster],
                                           amc_truth_pdg_code[ievent], amc_truth_is_prompt_photon[ievent])

                    if pt < 15:
                        continue

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
                    row['cluster_ecross_e'] = e_cross / e
                    row['cluster_ecross_emax'] = e_cross / e_max
                    row['weights'] = weight
                    row['centrality_v0m'] = centrality

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
    ntuplefilenames = sys.argv[1:-1]
    csvfilename = sys.argv[-1]
    main(ntuplefilenames, csvfilename)
