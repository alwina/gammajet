from __future__ import print_function

from array import array
import datetime
import fastjet
import os
import ROOT
import struct
import sys
import yaml

from alice_emcal import calculateShowerShapes5x5
from alice_triggers import getINT7TriggerIds, getCentralTriggerIds, getSemiCentralTriggerIds, getEMCEGATriggerIds, isEventSelected


def createAuxFile(ntuplefilename):
    rootfile = ROOT.TFile.Open(ntuplefilename, 'READ')
    tree = rootfile.Get('AliAnalysisTaskNTGJ/_tree_event')
    nevents = tree.GetEntries()

    outncluster = array('i', [0])
    cluster_5x5all = array('f', 17664 * [0])

    isINT7 = array('B', [0])
    isCentral = array('B', [0])
    isSemiCentral = array('B', [0])
    isEMCEGA = array('B', [0])

    njet_ak02tpc = array('i', [0])
    jet_ak02tpc_pt_raw = array('f', 10000 * [0])
    jet_ak02tpc_eta = array('f', 10000 * [0])
    jet_ak02tpc_phi = array('f', 10000 * [0])
    jet_ak02tpc_area = array('f', 10000 * [0])
    jet_ak02tpc_multiplicity_raw = array('H', 10000 * [0])

    outfile = ROOT.TFile.Open(ntuplefilename.replace('.root', '_AUX.root'), 'RECREATE')
    outtree = ROOT.TTree('ntupleaux', 'ntupleaux')
    outtree.Branch('ncluster', outncluster, 'ncluster/i')
    outtree.Branch('cluster_5x5all', cluster_5x5all, 'cluster5x5all[ncluster]/F')
    outtree.Branch('isINT7', isINT7, 'isINT7/O')
    outtree.Branch('isCentral', isCentral, 'isCentral/O')
    outtree.Branch('isSemiCentral', isSemiCentral, 'isSemiCentral/O')
    outtree.Branch('isEMCEGA', isEMCEGA, 'isEMCEGA/O')
    outtree.Branch('njet_ak02tpc', njet_ak02tpc, 'njet_ak02tpc/i')
    # this naming is in line with the ntuplizer, but it is with the UE subtracted
    outtree.Branch('jet_ak02tpc_pt_raw', jet_ak02tpc_pt_raw, 'jet_ak02tpc_pt_raw[njet_ak02tpc]/F')
    outtree.Branch('jet_ak02tpc_eta', jet_ak02tpc_eta, 'jet_ak02tpc_eta[njet_ak02tpc]/F')
    outtree.Branch('jet_ak02tpc_phi', jet_ak02tpc_phi, 'jet_ak02tpc_phi[njet_ak02tpc]/F')
    outtree.Branch('jet_ak02tpc_area', jet_ak02tpc_area, 'jet_ak02tpc_area[njet_ak02tpc]/F')
    outtree.Branch('jet_ak02tpc_multiplicity_raw', jet_ak02tpc_multiplicity_raw, 'jet_ak02tpc_multiplicity_raw[njet_ak02tpc]/s')

    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.2)
    areadef = fastjet.AreaDefinition(fastjet.VoronoiAreaSpec())

    # track previous run number so that we don't have to recalculate the trigger IDs if it didn't change
    previous_run_number = -1

    for ievent in range(nevents):
        if ievent % 25000 == 0:
            print('{0} Processing event {1}/{2}'.format(datetime.datetime.now(), ievent, nevents))
        tree.GetEntry(ievent)

        run_number = getattr(tree, 'run_number')
        trigger_mask = getattr(tree, 'trigger_mask')

        ncluster = getattr(tree, 'ncluster')
        cluster_e = getattr(tree, 'cluster_e')
        cluster_cell_id_max = getattr(tree, 'cluster_cell_id_max')
        cell_e = getattr(tree, 'cell_e')
        cell_cluster_index = getattr(tree, 'cell_cluster_index')

        ntrack = getattr(tree, 'ntrack')
        track_pt = getattr(tree, 'track_pt')
        track_eta = getattr(tree, 'track_eta')
        track_phi = getattr(tree, 'track_phi')
        track_e = getattr(tree, 'track_e')
        ue_estimate_tpc_const = getattr(tree, 'ue_estimate_tpc_const')

        # trigger info
        if run_number != previous_run_number:
            kINT7TriggerIds = getINT7TriggerIds(run_number)
            kCentralTriggerIds = getCentralTriggerIds(run_number)
            kSemiCentralTriggerIds = getSemiCentralTriggerIds(run_number)
            kEMCEGATriggerIds = getEMCEGATriggerIds(run_number)
            previous_run_number = run_number

        # combine the trigger masks into one number
        if len(trigger_mask) > 1:
            triggerMask = struct.unpack('2L', trigger_mask)[0] + (struct.unpack('2L', trigger_mask)[1] << 50)
        else:
            triggerMask = struct.unpack('1L', trigger_mask)[0]

        isINT7[0] = isEventSelected(kINT7TriggerIds, triggerMask)
        isCentral[0] = isEventSelected(kCentralTriggerIds, triggerMask)
        isSemiCentral[0] = isEventSelected(kSemiCentralTriggerIds, triggerMask)
        isEMCEGA[0] = isEventSelected(kEMCEGATriggerIds, triggerMask)

        # cluster info
        for icluster in range(ncluster):
            e = cluster_e[icluster]
            cell_id_max = cluster_cell_id_max[icluster]
            showershapes5x5 = calculateShowerShapes5x5(icluster, e, cell_id_max, cell_e, cell_cluster_index)
            cluster_5x5all[icluster] = showershapes5x5['5x5all']

        outncluster[0] = ncluster

        # jet info
        pjs = []
        for itrack in range(ntrack):
            pt = track_pt[itrack]
            eta = track_eta[itrack]
            phi = track_phi[itrack]
            e = track_e[itrack]

            # cut obviously nonsense tracks
            if abs(eta) > 0.9:
                continue

            if pt > 100000:
                continue

            # this cut should already exist, but check just in case
            if pt < 0.15:
                continue

            lv = ROOT.Math.PtEtaPhiEVector(pt, eta, phi, e)
            pjs.append(fastjet.PseudoJet(lv.Px(), lv.Py(), lv.Pz(), lv.E()))

        # need to keep this in scope in order for other stuff to work
        csa = fastjet.ClusterSequenceArea(pjs, jetdef, areadef)
        jets = csa.inclusive_jets()
        for ijet, jet in enumerate(jets):
            jet_ak02tpc_pt_raw[ijet] = jet.perp() - (ue_estimate_tpc_const * jet.area())
            jet_ak02tpc_eta[ijet] = jet.eta()
            jet_ak02tpc_phi[ijet] = jet.phi_std()
            jet_ak02tpc_area[ijet] = jet.area()
            jet_ak02tpc_multiplicity_raw[ijet] = len(jet.constituents())

        njet_ak02tpc[0] = len(jets)

        outtree.Fill()

    outfile.Write()
    outfile.Close()
    rootfile.Close()


def main(ntuplefilenames):
    for (ifile, ntuplefilename) in enumerate(ntuplefilenames):
        print('{0} Creating aux file for {1} ({2}/{3})'.format(datetime.datetime.now(), os.path.basename(ntuplefilename), ifile + 1, len(ntuplefilenames)))
        createAuxFile(ntuplefilename)
    print('Done')


if __name__ == '__main__':
    configfilename = sys.argv[1]
    with open(configfilename) as configfile:
        config = yaml.safe_load(configfile)
    ntuplefilenames = []
    ntuplefilenames.extend(config['filelists']['ntuples']['data'])
    ntuplefilenames.extend(config['filelists']['ntuples']['gjmc'])
    ntuplefilenames.extend(config['filelists']['ntuples']['jjmc'])
    main(ntuplefilenames)
