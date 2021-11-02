from __future__ import print_function

from array import array
import datetime
import os
import ROOT
import root_numpy as rnp
import sys
import yaml

from alice_emcal import calculateShowerShapes5x5
from alice_triggers import getINT7TriggerIds, getCentralTriggerIds, getSemiCentralTriggerIds, getEMCEGATriggerIds, isEventSelected


def createAuxFile(ntuplefilename):
    rootfile = ROOT.TFile.Open(ntuplefilename, 'READ')
    tree = rootfile.Get('AliAnalysisTaskNTGJ/_tree_event')
    atrigger_mask = rnp.tree2array(tree, branches='trigger_mask')
    nevents = tree.GetEntries()

    outncluster = array('i', [0])
    cluster_5x5all = array('f', 17664 * [0])
    isINT7 = array('i', [0])
    isCentral = array('i', [0])
    isSemiCentral = array('i', [0])
    isEMCEGA = array('i', [0])

    outfile = ROOT.TFile.Open(ntuplefilename.replace('.root', '_AUX.root'), 'RECREATE')
    outtree = ROOT.TTree('ntupleaux', 'ntupleaux')
    outtree.Branch('ncluster', outncluster, 'ncluster/i')
    outtree.Branch('cluster_5x5all', cluster_5x5all, 'cluster5x5all[ncluster]/F')
    outtree.Branch('isINT7', isINT7, 'isINT7/O')
    outtree.Branch('isCentral', isCentral, 'isCentral/O')
    outtree.Branch('isSemiCentral', isSemiCentral, 'isSemiCentral/O')
    outtree.Branch('isEMCEGA', isEMCEGA, 'isEMCEGA/O')

    # track previous run number so that we don't have to recalculate the trigger IDs if it didn't change
    previous_run_number = -1

    for ievent in range(nevents):
        tree.GetEntry(ievent)
        ncluster = getattr(tree, 'ncluster')
        cluster_e = getattr(tree, 'cluster_e')
        cluster_cell_id_max = getattr(tree, 'cluster_cell_id_max')
        run_number = getattr(tree, 'run_number')
        cell_e = getattr(tree, 'cell_e')
        cell_cluster_index = getattr(tree, 'cell_cluster_index')

        if run_number != previous_run_number:
            kINT7TriggerIds = getINT7TriggerIds(run_number)
            kCentralTriggerIds = getCentralTriggerIds(run_number)
            kSemiCentralTriggerIds = getSemiCentralTriggerIds(run_number)
            kEMCEGATriggerIds = getEMCEGATriggerIds(run_number)
            previous_run_number = run_number

        # combine the trigger masks into one number
        triggerMask = atrigger_mask[ievent][0] + (int(atrigger_mask[ievent][1]) << 50)

        isINT7[0] = isEventSelected(kINT7TriggerIds, triggerMask)
        isCentral[0] = isEventSelected(kCentralTriggerIds, triggerMask)
        isSemiCentral[0] = isEventSelected(kSemiCentralTriggerIds, triggerMask)
        isEMCEGA[0] = isEventSelected(kEMCEGATriggerIds, triggerMask)

        for icluster in range(ncluster):
            e = cluster_e[icluster]
            cell_id_max = cluster_cell_id_max[icluster]
            showershapes5x5 = calculateShowerShapes5x5(icluster, e, cell_id_max, cell_e, cell_cluster_index)
            cluster_5x5all[icluster] = showershapes5x5['5x5all']

        outncluster[0] = ncluster
        outtree.Fill()

    outfile.Write()
    outfile.Close()
    rootfile.Close()


def main(ntuplefilenames):
    for (ifile, ntuplefilename) in enumerate(ntuplefilenames):
        print('{0} Creating aux file for {2} ({3}/{4})'.format(datetime.datetime.now(), os.path.basename(ntuplefilename), ifile + 1, len(ntuplefilenames)))
        createAuxFile(ntuplefilename)


if __name__ == '__main__':
    configfilename = sys.argv[1]
    with open(configfilename) as configfile:
        config = yaml.safe_load(configfile)
    ntuplefilenames = []
    ntuplefilenames.extend(config['filelists']['ntuples']['data'])
    ntuplefilenames.extend(config['filelists']['ntuples']['gjmc'])
    ntuplefilenames.extend(config['filelists']['ntuples']['jjmc'])
    main(ntuplefilenames)
