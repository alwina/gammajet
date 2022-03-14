from array import array
import datetime
import fastjet
import os
import ROOT
import sys
import yaml

from alice_emcal import calculateShowerShapes5x5
from alice_mc import getIsPrompt, pdgCodeToCharge
from alice_triggers import getINT7TriggerIds, getCentralTriggerIds, getSemiCentralTriggerIds, getEMCEGATriggerIds, isEventSelected
from utils import tBranchToArray


def createAuxFile(ntuplefilename, maxnevents=-1):
    # open ntuple and get tree
    rootfile = ROOT.TFile.Open(ntuplefilename, 'READ')
    tree = rootfile.Get('AliAnalysisTaskNTGJ/_tree_event')
    nevents = tree.GetEntries()
    if maxnevents > 0:
        nevents = min(nevents, maxnevents)

    # set up variables for new branches
    outncluster = array('i', [0])
    cluster_5x5all = array('f', 17664 * [0])
    cluster_is_prompt = array('B', 17664 * [0])

    isINT7 = array('B', [0])
    isCentral = array('B', [0])
    isSemiCentral = array('B', [0])
    isEMCEGA = array('B', [0])

    njet_ak02tpc = array('I', [0])
    jet_ak02tpc_pt_raw = array('f', 10000 * [0])
    jet_ak02tpc_eta = array('f', 10000 * [0])
    jet_ak02tpc_phi = array('f', 10000 * [0])
    jet_ak02tpc_area = array('f', 10000 * [0])
    jet_ak02tpc_multiplicity_raw = array('H', 10000 * [0])
    
    njet_ak02its = array('I', [0])
    jet_ak02its_pt_raw = array('f', 10000 * [0])
    jet_ak02its_eta = array('f', 10000 * [0])
    jet_ak02its_phi = array('f', 10000 * [0])
    jet_ak02its_area = array('f', 10000 * [0])
    jet_ak02its_multiplicity_raw = array('H', 10000 * [0])

    njet_charged_truth_ak02 = array('I', [0])
    jet_charged_truth_ak02_pt = array('f', 10000 * [0])
    jet_charged_truth_ak02_eta = array('f', 10000 * [0])
    jet_charged_truth_ak02_phi = array('f', 10000 * [0])
    jet_charged_truth_ak02_area = array('f', 10000 * [0])
    jet_charged_truth_ak02_multiplicity = array('H', 10000 * [0])

    # create aux file and tree
    outfile = ROOT.TFile.Open(ntuplefilename.replace('.root', '_AUX.root'), 'RECREATE')
    outtree = ROOT.TTree('ntupleaux', 'ntupleaux')

    # set up new branches
    outtree.Branch('ncluster', outncluster, 'ncluster/i')
    outtree.Branch('cluster_5x5all', cluster_5x5all, 'cluster5x5all[ncluster]/F')
    outtree.Branch('cluster_is_prompt', cluster_is_prompt, 'cluster_is_prompt[ncluster]/O')

    outtree.Branch('isINT7', isINT7, 'isINT7/O')
    outtree.Branch('isCentral', isCentral, 'isCentral/O')
    outtree.Branch('isSemiCentral', isSemiCentral, 'isSemiCentral/O')
    outtree.Branch('isEMCEGA', isEMCEGA, 'isEMCEGA/O')

    outtree.Branch('njet_ak02tpc', njet_ak02tpc, 'njet_ak02tpc/i')
    outtree.Branch('jet_ak02tpc_pt_raw', jet_ak02tpc_pt_raw, 'jet_ak02tpc_pt_raw[njet_ak02tpc]/F')
    outtree.Branch('jet_ak02tpc_eta', jet_ak02tpc_eta, 'jet_ak02tpc_eta[njet_ak02tpc]/F')
    outtree.Branch('jet_ak02tpc_phi', jet_ak02tpc_phi, 'jet_ak02tpc_phi[njet_ak02tpc]/F')
    outtree.Branch('jet_ak02tpc_area', jet_ak02tpc_area, 'jet_ak02tpc_area[njet_ak02tpc]/F')
    outtree.Branch('jet_ak02tpc_multiplicity_raw', jet_ak02tpc_multiplicity_raw, 'jet_ak02tpc_multiplicity_raw[njet_ak02tpc]/s')
    
    outtree.Branch('njet_ak02its', njet_ak02its, 'njet_ak02its/i')
    outtree.Branch('jet_ak02its_pt_raw', jet_ak02its_pt_raw, 'jet_ak02its_pt_raw[njet_ak02its]/F')
    outtree.Branch('jet_ak02its_eta', jet_ak02its_eta, 'jet_ak02its_eta[njet_ak02its]/F')
    outtree.Branch('jet_ak02its_phi', jet_ak02its_phi, 'jet_ak02its_phi[njet_ak02its]/F')
    outtree.Branch('jet_ak02its_area', jet_ak02its_area, 'jet_ak02its_area[njet_ak02its]/F')
    outtree.Branch('jet_ak02its_multiplicity_raw', jet_ak02its_multiplicity_raw, 'jet_ak02its_multiplicity_raw[njet_ak02its]/s')

    outtree.Branch('njet_charged_truth_ak02', njet_charged_truth_ak02, 'njet_charged_truth_ak02/i')
    outtree.Branch('jet_charged_truth_ak02_pt', jet_charged_truth_ak02_pt, 'jet_charged_truth_ak02_pt[njet_charged_truth_ak02]/F')
    outtree.Branch('jet_charged_truth_ak02_eta', jet_charged_truth_ak02_eta, 'jet_charged_truth_ak02_eta[njet_charged_truth_ak02]/F')
    outtree.Branch('jet_charged_truth_ak02_phi', jet_charged_truth_ak02_phi, 'jet_charged_truth_ak02_phi[njet_charged_truth_ak02]/F')
    outtree.Branch('jet_charged_truth_ak02_area', jet_charged_truth_ak02_area, 'jet_charged_truth_ak02_area[njet_charged_truth_ak02]/F')
    outtree.Branch('jet_charged_truth_ak02_multiplicity', jet_charged_truth_ak02_multiplicity, 'jet_charged_truth_ak02_multiplicity[njet_charged_truth_ak02]/s')

    # set jet reconstruction parameters
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.2)
    areadef = fastjet.AreaDefinition(fastjet.VoronoiAreaSpec())

    # track previous run number so that we don't have to recalculate the trigger IDs if it didn't change
    previous_run_number = -1

    for ievent in range(nevents):
        if ievent % 25000 == 0:
            print('{0} Processing event {1}/{2}'.format(datetime.datetime.now(), ievent, nevents))
        tree.GetEntry(ievent)

        # retrieve event info
        run_number = getattr(tree, 'run_number')
        ue_estimate_tpc_const = getattr(tree, 'ue_estimate_tpc_const')
        ue_estimate_its_const = getattr(tree, 'ue_estimate_its_const')

        # branch_trigger_mask needs to be parsed further,
        # but how it's handled depends on its size, which is variable
        # therefore for now, we just grab the TBranch
        branch_trigger_mask = getattr(tree, 'trigger_mask')

        # retrieve cluster info
        ncluster = getattr(tree, 'ncluster')
        cluster_e = getattr(tree, 'cluster_e')
        cluster_cell_id_max = getattr(tree, 'cluster_cell_id_max')
        cluster_nmc_truth = tBranchToArray(getattr(tree, 'cluster_nmc_truth'), 'i', ncluster)
        cluster_mc_truth_index = tBranchToArray(getattr(tree, 'cluster_mc_truth_index'), 's', (ncluster, 32))

        # retrieve cell info
        cell_e = getattr(tree, 'cell_e')
        cell_cluster_index = getattr(tree, 'cell_cluster_index')

        # retrieve MC info
        nmc_truth = getattr(tree, 'nmc_truth')
        mc_truth_pt = getattr(tree, 'mc_truth_pt')
        mc_truth_eta = getattr(tree, 'mc_truth_eta')
        mc_truth_phi = getattr(tree, 'mc_truth_phi')
        mc_truth_e = getattr(tree, 'mc_truth_e')
        mc_truth_pdg_code = tBranchToArray(getattr(tree, 'mc_truth_pdg_code'), 'S', nmc_truth)
        try:
            mc_truth_is_prompt_photon = tBranchToArray(getattr(tree, 'mc_truth_is_prompt_photon'), 'O', nmc_truth)
        except AttributeError:
            mc_truth_is_prompt_photon = None

        # retrieve track info
        ntrack = getattr(tree, 'ntrack')
        track_pt = getattr(tree, 'track_pt')
        track_eta = getattr(tree, 'track_eta')
        track_phi = getattr(tree, 'track_phi')
        track_e = getattr(tree, 'track_e')

        # trigger info
#         if run_number != previous_run_number:
#             kINT7TriggerIds = getINT7TriggerIds(run_number)
#             kCentralTriggerIds = getCentralTriggerIds(run_number)
#             kSemiCentralTriggerIds = getSemiCentralTriggerIds(run_number)
#             kEMCEGATriggerIds = getEMCEGATriggerIds(run_number)
#             previous_run_number = run_number

#         # combine the trigger masks into one number
#         if len(branch_trigger_mask) > 1:
#             fullTriggerMask = tBranchToArray(branch_trigger_mask, 'l', 2)
#             triggerMask = fullTriggerMask[0] + (fullTriggerMask[1] << 50)
#         else:
#             triggerMask = tBranchToArray(branch_trigger_mask, 'l', 1)

#         isINT7[0] = isEventSelected(kINT7TriggerIds, triggerMask)
#         isCentral[0] = isEventSelected(kCentralTriggerIds, triggerMask)
#         isSemiCentral[0] = isEventSelected(kSemiCentralTriggerIds, triggerMask)
#         isEMCEGA[0] = isEventSelected(kEMCEGATriggerIds, triggerMask)
        isINT7[0] = False
        isCentral[0] = False
        isSemiCentral[0] = False
        isEMCEGA[0] = False

        # cluster info
        for icluster in range(ncluster):
            e = cluster_e[icluster]
            cell_id_max = cluster_cell_id_max[icluster]
            showershapes5x5 = calculateShowerShapes5x5(icluster, e, cell_id_max, cell_e, cell_cluster_index)
            cluster_5x5all[icluster] = showershapes5x5['5x5all']
            if mc_truth_is_prompt_photon is None:
                cluster_is_prompt[icluster] = False
            else:
                cluster_is_prompt[icluster] = getIsPrompt(cluster_nmc_truth[icluster], cluster_mc_truth_index[icluster],
                                                          mc_truth_pdg_code, mc_truth_is_prompt_photon)

        outncluster[0] = ncluster

        # jet R=0.2 info
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

            # this cut should already exist, but do it just in case
            if pt < 0.15:
                continue

            lv = ROOT.Math.PtEtaPhiEVector(pt, eta, phi, e)
            pjs.append(fastjet.PseudoJet(lv.Px(), lv.Py(), lv.Pz(), lv.E()))

        # need to keep this (csa) in scope in order for other stuff to work
        csa = fastjet.ClusterSequenceArea(pjs, jetdef, areadef)
        jets = csa.inclusive_jets()
        
        # for pp (17q, 18b10[ab], 18l2[ab]) and eventually also for p-Pb, these are ITS-only track jets. for everything else, they are TPC jets.
        if '17q' in ntuplefilename or '18b10' in ntuplefilename or '18l2' in ntuplefilename:
            njet_ak02its[0] = len(jets)
            for ijet, jet in enumerate(jets):
                jet_ak02its_pt_raw[ijet] = jet.perp() - (ue_estimate_its_const * jet.area())
                jet_ak02its_eta[ijet] = jet.eta()
                jet_ak02its_phi[ijet] = jet.phi_std()
                jet_ak02its_area[ijet] = jet.area()
                jet_ak02its_multiplicity_raw[ijet] = len(jet.constituents())
        else:
            njet_ak02tpc[0] = len(jets)
            for ijet, jet in enumerate(jets):
                jet_ak02tpc_pt_raw[ijet] = jet.perp() - (ue_estimate_tpc_const * jet.area())
                jet_ak02tpc_eta[ijet] = jet.eta()
                jet_ak02tpc_phi[ijet] = jet.phi_std()
                jet_ak02tpc_area[ijet] = jet.area()
                jet_ak02tpc_multiplicity_raw[ijet] = len(jet.constituents())

        # truth charged jet R=0.2 info
        pjs = []
        for imc in range(nmc_truth):
            # see note in alice_mc.py for why we use the pdg_code and not mc_truth_charge
            pdg_code = mc_truth_pdg_code[imc]
            charge = pdgCodeToCharge[pdg_code]
            if charge == 0:
                continue

            pt = mc_truth_pt[imc]
            eta = mc_truth_eta[imc]
            phi = mc_truth_phi[imc]
            e = mc_truth_e[imc]

            lv = ROOT.Math.PtEtaPhiEVector(pt, eta, phi, e)
            pjs.append(fastjet.PseudoJet(lv.Px(), lv.Py(), lv.Pz(), lv.E()))

        csa = fastjet.ClusterSequenceArea(pjs, jetdef, areadef)
        jets = csa.inclusive_jets()
        njet_charged_truth_ak02[0] = len(jets)
        for ijet, jet in enumerate(jets):
            jet_charged_truth_ak02_pt[ijet] = jet.perp()
            jet_charged_truth_ak02_eta[ijet] = jet.eta()
            jet_charged_truth_ak02_phi[ijet] = jet.phi_std()
            jet_charged_truth_ak02_area[ijet] = jet.area()
            jet_charged_truth_ak02_multiplicity[ijet] = len(jet.constituents())

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
    filetype = sys.argv[2]

    if filetype not in ('data', 'gjmc', 'jjmc'):
        print('File type {0} not recognized; must be one of data, gjmc, jjmc'.format(filetype))
        exit()

    with open(configfilename) as configfile:
        config = yaml.safe_load(configfile)

    ntuplefilenames = config['filelists']['ntuples'][filetype]
    main(ntuplefilenames)
