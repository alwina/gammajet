#include <TLorentzVector.h>
#include <iostream>
#include <math.h>

#include "same_event.h"
#include "config_parser.h"
#include "shared_defs.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cout <<  "Format: [command] [config file] [nevents (optional)]" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (argc > 2) {
		nevents_max = std::stol(argv[2]);
	} else {
		nevents_max = 999999999999999;
	}

	/*--------------------------------------------------------------
	Initial setup
	--------------------------------------------------------------*/
	YAML::Node configrunperiod = YAML::LoadFile(argv[1]);
	allconfigs.push_back(configrunperiod);
	YAML::Node configsystem = YAML::LoadFile(configrunperiod["systemconfig"].as<std::string>());
	allconfigs.push_back(configsystem);
	YAML::Node configglobal = YAML::LoadFile(configsystem["globalconfig"].as<std::string>());
	allconfigs.push_back(configglobal);
	parseConfig();
	printCutSummary();

	// set up THnSparses
	initializeTHnSparses();
	Double_t trig[ndimTrig];
	Double_t corr[ndimCorr];

	/*--------------------------------------------------------------
	Loop through files
	--------------------------------------------------------------*/
	YAML::Node filenames = configrunperiod["filelists"]["ntuples"]["data"];
	for (YAML::const_iterator fileit = filenames.begin(); fileit != filenames.end(); fileit++) {
		std::string root_filename = fileit->as<std::string>();
		openFilesAndGetTTrees(root_filename);
        setBranchAddresses();

		/*--------------------------------------------------------------
		Loop through events
		--------------------------------------------------------------*/
		nevents = std::min(_tree_event->GetEntries(), nevents_max);
		for (long ievent = 0; ievent < nevents; ievent++) {
			// for some reason, loading these 2 events from this ntuple causes a segfault
			if (root_filename == "/global/project/projectdirs/alice/NTuples/PbPb/15o_pass2_cluster15.root") {
				if (ievent == 7894 || ievent == 7895) {
					std::cout << std::endl << "skipping event " << ievent;
					if (ievent == 7895) {
						std::cout << std::endl;
					}
					continue;
				}
			}

			// load this event
			_tree_event->GetEntry(ievent);
			if (auxfile != NULL) {
				auxtree->GetEntry(ievent);
			}
			fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nevents);

			// event selection
			if (abs(primary_vertex[2]) > 10) continue;
			if (primary_vertex[2] == 0.00) continue;
			if (do_pile && is_pileup_from_spd_5_08) continue;

			/*--------------------------------------------------------------
			Loop through clusters
			--------------------------------------------------------------*/
			for (long icluster = 0; icluster < ncluster; icluster++) {
				// apply cluster cuts
				if (rejectCluster(icluster)) continue;

				// determine whether the cluster is isolated
				float isolation = getIsolation(icluster);
				isIsolated = GetIsIsolated(isolation, centrality_v0m, isoconfig);
				if (not(isIsolated)) continue;

				// determine whether it is SR or BR (or neither), calculate purity, and fill trigger THnSparse
				float shower = getShower(icluster);
				isSignal = (shower > srmin) and (shower < srmax);
				isBackground = (shower > brmin) and (shower < brmax);
				if (not(isSignal or isBackground)) continue;

				float purity = getPurity(cluster_pt[icluster], centrality_v0m, purityconfig);

				trig[0] = centrality_v0m;
				trig[1] = cluster_pt[icluster];

				if (isSignal) {
					purity_weight = 1.0 / purity;
					hTrigSR->Fill(trig);
				}

				if (isBackground) {
					purity_weight = 1.0 / purity - 1;
					hTrigBR->Fill(trig);
				}

				/*--------------------------------------------------------------
				Loop through jets
				--------------------------------------------------------------*/
				for (long ijet = 0; ijet < njet; ijet++) {
					// apply jet cuts
					if (not(jet_pt_raw[ijet] > jet_pt_min and jet_pt_raw[ijet] < jet_pt_max)) continue;
					if (not(abs(jet_eta[ijet]) < jet_eta_max)) continue;

					// calculate observables
					float deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[icluster] - jet_phi[ijet]));
					float jetpt = jet_pt_raw[ijet];
					float ptratio = jetpt / cluster_pt[icluster];

					corr[0] = centrality_v0m;
					corr[1] = cluster_pt[icluster];
					corr[2] = deltaphi;
					corr[3] = jetpt;
					corr[4] = ptratio;

					// fill correlation THnSparse
					if (isSignal) {
						hCorrSR->Fill(corr, purity_weight);
					}

					if (isBackground) {
						hCorrBR->Fill(corr, purity_weight);
					}
				} // end jet loop
			} // end cluster loop
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
	// Write to fout
	TFile* fout;
	fout = new TFile((TString) configrunperiod["filelists"]["correlations"]["sameevent"].as<std::string>(), "RECREATE");
	std::cout << "Writing to file" << std::endl;

	hTrigSR->Write();
	hCorrSR->Write();
	hTrigBR->Write();
	hCorrBR->Write();

	fout->Close();
	std::cout << "Ending" << std::endl;
	return EXIT_SUCCESS;
}

/*--------------------------------------------------------------
Set up THnSparses
hCorrSR: cluster-jet correlations for the signal region
hTrigSR: counting the number of clusters in each bin in the signal region
hCorrBR: cluster-jet correlations for the bkg region
hTrigBR: counting the number of clusters in each bin in the bkg region
--------------------------------------------------------------*/
void initializeTHnSparses()
{
	int nbinsClusterPt = 2 * (cluster_pt_max - cluster_pt_min);

	// dimensions: centrality, cluster pT
	ndimTrig = 2;
	Int_t nbinsTrig[ndimTrig];
	Double_t minbinsTrig[ndimTrig];
	Double_t maxbinsTrig[ndimTrig];

	// centrality
	nbinsTrig[0] = 10;
	minbinsTrig[0] = 0;
	maxbinsTrig[0] = 100;

	// cluster pT
	nbinsTrig[1] = nbinsClusterPt;
	minbinsTrig[1] = cluster_pt_min;
	maxbinsTrig[1] = cluster_pt_max;

	hTrigSR = new THnSparseF("hTrigSR", "Number of clusters (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
	hTrigBR = new THnSparseF("hTrigBR", "Number of clusters (BR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

	// dimensions: centrality, cluster pT, delta phi, jet pT, pT ratio
	ndimCorr = 5;
	Int_t nbinsCorr[ndimCorr];
	Double_t minbinsCorr[ndimCorr];
	Double_t maxbinsCorr[ndimCorr];

	// centrality
	nbinsCorr[0] = 10;
	minbinsCorr[0] = 0;
	maxbinsCorr[0] = 100;

	// cluster pT
	nbinsCorr[1] = nbinsClusterPt;
	minbinsCorr[1] = cluster_pt_min;
	maxbinsCorr[1] = cluster_pt_max;

	// deltaphi
	nbinsCorr[2] = 120;
	minbinsCorr[2] = deltaphi_min;
	maxbinsCorr[2] = deltaphi_max;

	// jetpt
	nbinsCorr[3] = 120;
	minbinsCorr[3] = jetpt_min;
	maxbinsCorr[3] = jetpt_max;

	// ptratio
	nbinsCorr[4] = 120;
	minbinsCorr[4] = ptratio_min;
	maxbinsCorr[4] = ptratio_max;

	hCorrSR = new THnSparseF("hCorrSR", "Correlations (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrBR = new THnSparseF("hCorrBR", "Correlations (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrSR->Sumw2();
	hCorrBR->Sumw2();
}

// Apply cluster cuts
bool rejectCluster(int icluster)
{
	if (not(cluster_pt[icluster] > cluster_pt_min and cluster_pt[icluster] < cluster_pt_max)) return true;
	if (not(abs(cluster_eta[icluster]) < cluster_eta_max)) return true;
	if (not(cluster_ncell[icluster] >= cluster_ncell_min)) return true;
	if (not(cluster_e_cross[icluster] / cluster_e_max[icluster] > cluster_ecross_emax_min)) return true;
	if (not(cluster_distance_to_bad_channel[icluster] >= cluster_dbc_min)) return true;
	if (not(cluster_nlocal_maxima[icluster] < cluster_nlm_max)) return true;
	if (not(abs(cluster_tof[icluster]) < cluster_tof_max)) return true;
	return false;
}

// Calculate isolation
float getIsolation(int icluster)
{
    float isolation;
	if (isovar == "cluster_iso_tpc_04") {
		isolation = cluster_iso_tpc_04[icluster];
	} else if (isovar == "cluster_iso_its_04") {
		isolation = cluster_iso_its_04[icluster];
	} else if (isovar == "cluster_iso_its_04_sub") {
		isolation = cluster_iso_its_04[icluster] + cluster_iso_its_04_ue[icluster] - ue_estimate_its_const * M_PI * 0.4 * 0.4;
	} else if (isovar == "cluster_iso_tpc_02_sub") {
		isolation = cluster_iso_tpc_02[icluster] + cluster_iso_tpc_02_ue[icluster] - ue_estimate_tpc_const * M_PI * 0.2 * 0.2;
	} else if (isovar == "cluster_iso_tpc_04_sub") {
		isolation = cluster_iso_tpc_04[icluster] + cluster_iso_tpc_04_ue[icluster] - ue_estimate_tpc_const * M_PI * 0.4 * 0.4;
	} else if (isovar == "cluster_frixione_tpc_04_02") {
		isolation = cluster_frixione_tpc_04_02[icluster];
	} else if (isovar == "cluster_frixione_its_04_02") {
		isolation = cluster_frixione_its_04_02[icluster];
	} else {
		std::cout << "ERROR: Isolation variable " << isovar << " not recognized. Aborting" << std::endl;
		exit(EXIT_FAILURE);
	}
	return isolation;
}

// Get shower shape value
float getShower(int icluster)
{
    float shower;
	if (shower_shape == "cluster_Lambda") {
		shower = cluster_lambda_square[icluster][0];
	} else if (shower_shape == "cluster_NN1") {
		shower = cluster_s_nphoton[icluster][1];
	} else if (shower_shape == "cluster_emax_over_e") {
		shower = cluster_e_max[icluster] / cluster_e[icluster];
	} else if (shower_shape == "cluster_5x5all") {
		shower = cluster_5x5all[icluster];
	} else {
		std::cout << "ERROR: Shower shape variable " << shower_shape << " not recognized. Aborting" << std::endl;
		exit(EXIT_FAILURE);
	}
	return shower;
}

// Print cut summary
void printCutSummary()
{
	std::cout << "Cluster pT range: " << cluster_pt_min << "-" << cluster_pt_max << std::endl;
	std::cout << "Cluster eta max: " << cluster_eta_max << std::endl;
	std::cout << "Cluster ncell min: " << cluster_ncell_min << std::endl;
	std::cout << "Cluster Ecross/Emax min: " << cluster_ecross_emax_min << std::endl;
	std::cout << "Cluster dist to bad channel min: " << cluster_dbc_min << std::endl;
	std::cout << "Cluster nlocal maxima max: " << cluster_nlm_max << std::endl;
	std::cout << "Cluster TOF max: " << cluster_tof_max << std::endl;
	std::cout << "Shower shape SR range: " << srmin << "-" << srmax << std::endl;
	std::cout << "Shower shape BR range: " << brmin << "-" << brmax << std::endl;
	std::cout << "Jet type: " << jettype << std::endl;
	std::cout << "Jet pT range: " << jet_pt_min << "-" << jet_pt_max << std::endl;
	std::cout << "Jet eta max: " << jet_eta_max << std::endl;
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
	// the aux file is only needed in certain situations, so check for those situations
	// things should still work even without the aux file otherwise
	if (auxfile == NULL) {
		if (jettype == "ak02tpc") {
			std::cout << "ERROR: No aux file; cannot get jet type " << jettype << std::endl;
			exit(EXIT_FAILURE);
		}
		if (shower_shape == "cluster_5x5all") {
			std::cout << "ERROR: No aux file; cannot get shower shape " << shower_shape << std::endl;
			exit(EXIT_FAILURE);
		}
	} else {
		auxtree = dynamic_cast<TTree*>(auxfile->Get("ntupleaux"));
	}
}

// Set the branch addresses of the branches in the TTrees
void setBranchAddresses()
{
	_tree_event->SetBranchStatus("*mc*", 0);

	// event Addresses
	_tree_event->SetBranchAddress("primary_vertex", primary_vertex);
	_tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
	_tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
	_tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
	_tree_event->SetBranchAddress("centrality_v0m", &centrality_v0m);

	//track Addresses
	_tree_event->SetBranchAddress("primary_vertex", primary_vertex);
	_tree_event->SetBranchAddress("ntrack", &ntrack);
	_tree_event->SetBranchAddress("track_e", track_e);
	_tree_event->SetBranchAddress("track_pt", track_pt);
	_tree_event->SetBranchAddress("track_eta", track_eta);
	_tree_event->SetBranchAddress("track_phi", track_phi);
	_tree_event->SetBranchAddress("track_eta_emcal", track_eta_emcal);
	_tree_event->SetBranchAddress("track_phi_emcal", track_phi_emcal);
	_tree_event->SetBranchAddress("track_quality", track_quality);
	_tree_event->SetBranchAddress("track_its_ncluster", &track_its_ncluster);
	_tree_event->SetBranchAddress("track_its_chi_square", &track_its_chi_square);
	_tree_event->SetBranchAddress("track_dca_xy", &track_dca_xy);
	_tree_event->SetBranchAddress("track_dca_z", &track_dca_z);

	//Cluster Addresses
	_tree_event->SetBranchAddress("ncluster", &ncluster);
	_tree_event->SetBranchAddress("cluster_e", cluster_e);
	_tree_event->SetBranchAddress("cluster_e_max", cluster_e_max);
	_tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
	_tree_event->SetBranchAddress("cluster_pt", cluster_pt);
	_tree_event->SetBranchAddress("cluster_eta", cluster_eta);
	_tree_event->SetBranchAddress("cluster_phi", cluster_phi);
	_tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
	_tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
	_tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
	_tree_event->SetBranchAddress("cluster_iso_tpc_02", cluster_iso_tpc_02);
	_tree_event->SetBranchAddress("cluster_iso_tpc_04", cluster_iso_tpc_04);
	_tree_event->SetBranchAddress("cluster_iso_its_04", cluster_iso_its_04);
	_tree_event->SetBranchAddress("cluster_frixione_tpc_04_02", cluster_frixione_tpc_04_02);
	_tree_event->SetBranchAddress("cluster_frixione_its_04_02", cluster_frixione_its_04_02);
	_tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
	_tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);

	_tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
	_tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
	_tree_event->SetBranchAddress("cell_e", cell_e);

	_tree_event->SetBranchAddress("cluster_tof", cluster_tof);
	_tree_event->SetBranchAddress("cluster_iso_its_04_ue", cluster_iso_its_04_ue);
	_tree_event->SetBranchAddress("cluster_iso_tpc_02_ue", cluster_iso_tpc_02_ue);
	_tree_event->SetBranchAddress("cluster_iso_tpc_04_ue", cluster_iso_tpc_04_ue);

	if (auxfile != NULL) {
		auxtree->SetBranchAddress("cluster_5x5all", cluster_5x5all);
	}

	// Jet addresses
	// switch based on jet type
	if (jettype == "ak04tpc") {
		_tree_event->SetBranchAddress("njet_ak04tpc", &njet);
		_tree_event->SetBranchAddress("jet_ak04tpc_pt_raw", jet_pt_raw);
		_tree_event->SetBranchAddress("jet_ak04tpc_eta", jet_eta);
		_tree_event->SetBranchAddress("jet_ak04tpc_phi", jet_phi);
	} else if (jettype == "ak02tpc") {
		auxtree->SetBranchAddress("njet_ak02tpc", &njet);
		auxtree->SetBranchAddress("jet_ak02tpc_pt_raw", jet_pt_raw);
		auxtree->SetBranchAddress("jet_ak02tpc_eta", jet_eta);
		auxtree->SetBranchAddress("jet_ak02tpc_phi", jet_phi);
	} else if (jettype == "ak04its") {
		_tree_event->SetBranchAddress("njet_ak04its", &njet);
		_tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_pt_raw);
		_tree_event->SetBranchAddress("jet_ak04its_eta", jet_eta);
		_tree_event->SetBranchAddress("jet_ak04its_phi", jet_phi);
	} else {
		std::cout << "ERROR: Jet type " << jettype << " not recognized. Aborting" << std::endl;
		exit(EXIT_FAILURE);
	}
}