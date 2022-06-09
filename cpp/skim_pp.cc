#include <TLorentzVector.h>
#include <iostream>
#include <math.h>

#include "skim_pp.h"
#include "config_parser.h"
#include "shared_defs.h"

int main(int argc, char *argv[])
{
	if (argc < 3) {
		std::cout <<  "Format: [command] [config file] [fileset]" << std::endl;
		exit(EXIT_FAILURE);
	}

	/*--------------------------------------------------------------
	Initial setup
	--------------------------------------------------------------*/
	std::string fileset = argv[2];
    if (fileset != "data" && fileset != "gjmc" && fileset != "jjmc" && fileset != "gjmcwithneutrals") {
        std::cout << "fileset " << fileset << " not recognized; must be one of data, gjmc, jjmc, gjmcwithneutrals" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    YAML::Node configrunperiod = YAML::LoadFile(argv[1]);
	allconfigs.push_back(configrunperiod);
	YAML::Node configsystem = YAML::LoadFile(configrunperiod["systemconfig"].as<std::string>());
	allconfigs.push_back(configsystem);
	YAML::Node configglobal = YAML::LoadFile(configsystem["globalconfig"].as<std::string>());
	allconfigs.push_back(configglobal);
	parseConfig();
	printCutSummary();

	/*--------------------------------------------------------------
	Loop through files
	--------------------------------------------------------------*/
	YAML::Node filenames = configrunperiod["filelists"]["ntuples"][fileset];

	// create output file - need to do it this way to not get the memory-resident TTree thing
	openFilesAndGetTTrees(filenames[0].as<std::string>());
	setBranchAddresses();
	TFile* fout = new TFile((TString) configrunperiod["filelists"]["skimmedntuples"][fileset].as<std::string>(), "RECREATE");
	TTree* outtree = _tree_event->CloneTree(0);
	setupNewBranches(outtree);

	for (YAML::const_iterator fileit = filenames.begin(); fileit != filenames.end(); fileit++) {
		std::string root_filename = fileit->as<std::string>();
		openFilesAndGetTTrees(root_filename);
		setBranchAddresses();

		/*--------------------------------------------------------------
		Loop through events
		--------------------------------------------------------------*/
		nevents = _tree_event->GetEntries();
		for (long ievent = 0; ievent < nevents; ievent++) {
			// load this event
			_tree_event->GetEntry(ievent);
			auxtree->GetEntry(ievent);

			fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nevents);

			copyBranchValues();

			// always fill the first event
			if (ievent == 0) {
				outtree->Fill();
			} else {
				// check to see if there is a sufficiently high-pT cluster to keep the event
				float maxClusterPtInEvent = -1;
				for (long icluster = 0; icluster < ncluster; icluster++) {
					maxClusterPtInEvent = std::max(maxClusterPtInEvent, cluster_pt[icluster]);
				}

				if (maxClusterPtInEvent < cluster_pt_min) continue;

				outtree->Fill();
			}

		} // end event loop

		// close files
		file->Close();
		if (auxfile != NULL) {
			auxfile->Close();
		}
		std::cout << std::endl;
	} // end file loop

	/*--------------------------------------------------------------
	Write outputs to file
	--------------------------------------------------------------*/
	outtree->AutoSave();

	fout->Close();
	std::cout << "Ending" << std::endl;
	return EXIT_SUCCESS;
}


// Print cut summary
void printCutSummary()
{
	std::cout << "Cluster pT min: " << cluster_pt_min << std::endl;
}

void openFilesAndGetTTrees(std::string root_filename)
{
	// open main ROOT file
	std::cout << "Opening " << root_filename << std::endl;
	file = TFile::Open((TString) root_filename);

	if (file == NULL) {
		std::cout << "Failed to open file" << std::endl;
		exit(EXIT_FAILURE);
	}

	_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

	if (_tree_event == NULL) {
		_tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
		if (_tree_event == NULL) {
			std::cout << " fail " << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// open aux file
	std::string aux_filename = root_filename.replace(root_filename.find(".root"), 5, "_AUX.root");
	std::cout << "Opening " << aux_filename << std::endl;
	auxfile = TFile::Open((TString) aux_filename);
	auxtree = dynamic_cast<TTree*>(auxfile->Get("ntupleaux"));
	// the aux file is only needed in certain situations, so check for those situations
	// things should still work even without the aux file otherwise
	if (auxfile == NULL) {
		std::cout << "ERROR: no aux file. Exiting" << std::endl;
		exit(EXIT_FAILURE);
	} else {
		auxtree = dynamic_cast<TTree*>(auxfile->Get("ntupleaux"));
		if (auxtree == NULL) {
			std::cout << "ERROR: no aux tree. Exiting" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

void setBranchAddresses()
{
	_tree_event->SetBranchAddress("ncluster", &ncluster);
	_tree_event->SetBranchAddress("cluster_pt", cluster_pt);

	auxtree->SetBranchAddress("cluster_5x5all", IN_cluster_5x5all);

	auxtree->SetBranchAddress("njet_ak02tpc", &IN_njet_ak02tpc);
	auxtree->SetBranchAddress("jet_ak02tpc_pt_raw", IN_jet_ak02tpc_pt_raw);
	auxtree->SetBranchAddress("jet_ak02tpc_eta", IN_jet_ak02tpc_eta);
	auxtree->SetBranchAddress("jet_ak02tpc_phi", IN_jet_ak02tpc_phi);
	auxtree->SetBranchAddress("jet_ak02tpc_area", IN_jet_ak02tpc_area);
	auxtree->SetBranchAddress("jet_ak02tpc_multiplicity_raw", IN_jet_ak02tpc_multiplicity_raw);

	auxtree->SetBranchAddress("njet_ak02its", &IN_njet_ak02its);
	auxtree->SetBranchAddress("jet_ak02its_pt_raw", IN_jet_ak02its_pt_raw);
	auxtree->SetBranchAddress("jet_ak02its_eta", IN_jet_ak02its_eta);
	auxtree->SetBranchAddress("jet_ak02its_phi", IN_jet_ak02its_phi);
	auxtree->SetBranchAddress("jet_ak02its_area", IN_jet_ak02its_area);
	auxtree->SetBranchAddress("jet_ak02its_multiplicity_raw", IN_jet_ak02its_multiplicity_raw);

	auxtree->SetBranchAddress("njet_charged_truth_ak02", &IN_njet_charged_truth_ak02);
	auxtree->SetBranchAddress("jet_charged_truth_ak02_pt", IN_jet_charged_truth_ak02_pt);
	auxtree->SetBranchAddress("jet_charged_truth_ak02_eta", IN_jet_charged_truth_ak02_eta);
	auxtree->SetBranchAddress("jet_charged_truth_ak02_phi", IN_jet_charged_truth_ak02_phi);
	auxtree->SetBranchAddress("jet_charged_truth_ak02_area", IN_jet_charged_truth_ak02_area);
	auxtree->SetBranchAddress("jet_charged_truth_ak02_multiplicity", IN_jet_charged_truth_ak02_multiplicity);

}

void setupNewBranches(TTree* outtree)
{
	outtree->Branch("ncluster", &outncluster, "ncluster/i");
	outtree->Branch("cluster_5x5all", cluster_5x5all, "cluster5x5all[ncluster]/F");

	outtree->Branch("njet_ak02tpc", &njet_ak02tpc, "njet_ak02tpc/i");
	outtree->Branch("jet_ak02tpc_pt_raw", jet_ak02tpc_pt_raw, "jet_ak02tpc_pt_raw[njet_ak02tpc]/F");
	outtree->Branch("jet_ak02tpc_eta", jet_ak02tpc_eta, "jet_ak02tpc_eta[njet_ak02tpc]/F");
	outtree->Branch("jet_ak02tpc_phi", jet_ak02tpc_phi, "jet_ak02tpc_phi[njet_ak02tpc]/F");
	outtree->Branch("jet_ak02tpc_area", jet_ak02tpc_area, "jet_ak02tpc_area[njet_ak02tpc]/F");
	outtree->Branch("jet_ak02tpc_multiplicity_raw", jet_ak02tpc_multiplicity_raw, "jet_ak02tpc_multiplicity_raw[njet_ak02tpc]/s");

	outtree->Branch("njet_ak02its", &njet_ak02its, "njet_ak02its/i");
	outtree->Branch("jet_ak02its_pt_raw", jet_ak02its_pt_raw, "jet_ak02its_pt_raw[njet_ak02its]/F");
	outtree->Branch("jet_ak02its_eta", jet_ak02its_eta, "jet_ak02its_eta[njet_ak02its]/F");
	outtree->Branch("jet_ak02its_phi", jet_ak02its_phi, "jet_ak02its_phi[njet_ak02its]/F");
	outtree->Branch("jet_ak02its_area", jet_ak02its_area, "jet_ak02its_area[njet_ak02its]/F");
	outtree->Branch("jet_ak02its_multiplicity_raw", jet_ak02its_multiplicity_raw, "jet_ak02its_multiplicity_raw[njet_ak02its]/s");

	outtree->Branch("njet_charged_truth_ak02", &njet_charged_truth_ak02, "njet_charged_truth_ak02/i");
	outtree->Branch("jet_charged_truth_ak02_pt", jet_charged_truth_ak02_pt, "jet_charged_truth_ak02_pt[njet_charged_truth_ak02]/F");
	outtree->Branch("jet_charged_truth_ak02_eta", jet_charged_truth_ak02_eta, "jet_charged_truth_ak02_eta[njet_charged_truth_ak02]/F");
	outtree->Branch("jet_charged_truth_ak02_phi", jet_charged_truth_ak02_phi, "jet_charged_truth_ak02_phi[njet_charged_truth_ak02]/F");
	outtree->Branch("jet_charged_truth_ak02_area", jet_charged_truth_ak02_area, "jet_charged_truth_ak02_area[njet_charged_truth_ak02]/F");
	outtree->Branch("jet_charged_truth_ak02_multiplicity", jet_charged_truth_ak02_multiplicity, "jet_charged_truth_ak02_multiplicity[njet_charged_truth_ak02]/s");
}

void copyBranchValues()
{
	outncluster = ncluster;
	for (int i = 0; i < ncluster; i++) {
		cluster_5x5all[i] = IN_cluster_5x5all[i];
	}

	njet_ak02tpc = IN_njet_ak02tpc;
	for (int i = 0; i < njet_ak02tpc; i++) {
		jet_ak02tpc_pt_raw[i] = IN_jet_ak02tpc_pt_raw[i];
		jet_ak02tpc_eta[i] = IN_jet_ak02tpc_eta[i];
		jet_ak02tpc_phi[i] = IN_jet_ak02tpc_phi[i];
		jet_ak02tpc_area[i] = IN_jet_ak02tpc_area[i];
		jet_ak02tpc_multiplicity_raw[i] = IN_jet_ak02tpc_multiplicity_raw[i];
	}

	njet_ak02its = IN_njet_ak02its;
	for (int i = 0; i < njet_ak02its; i++) {
		jet_ak02its_pt_raw[i] = IN_jet_ak02its_pt_raw[i];
		jet_ak02its_eta[i] = IN_jet_ak02its_eta[i];
		jet_ak02its_phi[i] = IN_jet_ak02its_phi[i];
		jet_ak02its_area[i] = IN_jet_ak02its_area[i];
		jet_ak02its_multiplicity_raw[i] = IN_jet_ak02its_multiplicity_raw[i];
	}

	njet_charged_truth_ak02 = IN_njet_charged_truth_ak02;
	for (int i = 0; i < njet_charged_truth_ak02; i++) {
		jet_charged_truth_ak02_pt[i] = IN_jet_charged_truth_ak02_pt[i];
		jet_charged_truth_ak02_eta[i] = IN_jet_charged_truth_ak02_eta[i];
		jet_charged_truth_ak02_phi[i] = IN_jet_charged_truth_ak02_phi[i];
		jet_charged_truth_ak02_area[i] = IN_jet_charged_truth_ak02_area[i];
		jet_charged_truth_ak02_multiplicity[i] = IN_jet_charged_truth_ak02_multiplicity[i];
	}
}