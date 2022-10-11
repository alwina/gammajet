#include <TLorentzVector.h>
#include <iostream>
#include <math.h>

#include "same_event_mc.h"
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
	float avg_eg_ntrial = 0;
	float weight = 0;

	/*--------------------------------------------------------------
	Loop through files
	--------------------------------------------------------------*/
    std::vector<std::string> filenames = configrunperiod["filelists"]["ntuples"]["gjmc"].as<std::vector<std::string>>();
    std::vector<std::string> jjmcfilenames = configrunperiod["filelists"]["ntuples"]["jjmc"].as<std::vector<std::string>>();
    filenames.insert(filenames.end(), jjmcfilenames.begin(), jjmcfilenames.end());
	for (auto & root_filename : filenames) {
		openFilesAndGetTTrees(root_filename);
		setBranchAddresses();
		avg_eg_ntrial = getAvgEgNtrial(root_filename);

		/*--------------------------------------------------------------
		Loop through events
		--------------------------------------------------------------*/
		nevents = std::min(_tree_event->GetEntries(), nevents_max);
		for (long ievent = 0; ievent < nevents; ievent++) {
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
			// if (!isINT7) continue;

			weight = eg_cross_section / avg_eg_ntrial;
            
			matchJetsInEvent();
            // if there are any reco jets matched with multiple truth jets, throw this away
            std::set<int> setMatchedReco(matchedReco.begin(), matchedReco.end());
            if (setMatchedReco.size() != matchedReco.size()) {
                // clear the vectors first
                matchedJetIndices.clear();
                unmatchedTruth.clear();
                unmatchedReco.clear();
                matchedReco.clear();
                continue;
            }

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
				isSignal = GetIsSignal(shower, centrality_v0m, ssconfig);
				isBackground = GetIsBackground(shower, centrality_v0m, ssconfig);
				if (not(isSignal or isBackground)) continue;
                
                // determine if this is a non-decay photon
                bool isNonDecay = getIsClusterNonDecay(icluster);

				float purity = getPurity(cluster_pt[icluster], centrality_v0m, purityconfig);

				trig[0] = centrality_v0m;
				trig[1] = cluster_pt[icluster];

				if (isSignal) {
					purity_weight = 1.0 / purity;
					hTrigSR->Fill(trig, weight);
                    
                    if (cluster_is_prompt[icluster]) {
                        hTrigSRPromptTruth->Fill(trig, weight);
                    }
                    
                    if (isNonDecay) {
                        hTrigSRNonDecay->Fill(trig, weight);
                    }
				}

				if (isBackground) {
					purity_weight = 1.0 / purity - 1;
					hTrigBR->Fill(trig, weight);

                    if (isNonDecay) {
                        hTrigBRNonDecay->Fill(trig, weight);
                    }
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
					corr[5] = cluster_pt[icluster] * sin(deltaphi);

					// fill correlation THnSparse
					if (isSignal) {
						hCorrSR->Fill(corr, weight * purity_weight);
						hCorr1ptSR->Fill(corr, weight * purity_weight / jetpt);
					}

					if (isBackground) {
						hCorrBR->Fill(corr, weight * purity_weight);
						hCorr1ptBR->Fill(corr, weight * purity_weight / jetpt);
					}
				} // end jet loop
                
				/*--------------------------------------------------------------
				Loop through matched jets, but only if it's a non-decay photon
				--------------------------------------------------------------*/
                if (not(isNonDecay)) continue;
                
				for (auto matchedIndex : matchedJetIndices) {
					int ireco = matchedIndex.second;
                    
					float j_reco_pt = jet_pt_raw[ireco];
					float j_reco_eta = jet_eta[ireco];
					float j_reco_phi = jet_phi[ireco];

					// apply jet cuts on the reco jets
					if (not(j_reco_pt > jet_pt_min and j_reco_pt < jet_pt_max)) continue;
					if (not(abs(j_reco_eta) < jet_eta_max)) continue;

					float deltaphi_reco = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[icluster] - j_reco_phi));
                    
					corr[0] = centrality_v0m;
					corr[1] = cluster_pt[icluster];
					corr[2] = deltaphi_reco;
					corr[3] = j_reco_pt;
					corr[4] = j_reco_pt / cluster_pt[icluster];
					corr[5] = cluster_pt[icluster] * sin(deltaphi_reco);
                    
					if (isSignal) {
						hCorrSRNonDecayMatchedJet->Fill(corr, weight * purity_weight);
						hCorr1ptSRNonDecayMatchedJet->Fill(corr, weight * purity_weight / j_reco_pt);
					}

					if (isBackground) {
						hCorrBRNonDecayMatchedJet->Fill(corr, weight * purity_weight);
						hCorr1ptBRNonDecayMatchedJet->Fill(corr, weight * purity_weight / j_reco_pt);
					}
                } // end matched jet loop
			} // end cluster loop
            
			matchedJetIndices.clear();
			unmatchedTruth.clear();
			unmatchedReco.clear();
            matchedReco.clear();
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
	fout = new TFile("root/sameEvent20g3.root", "RECREATE");
	std::cout << "Writing to file" << std::endl;

    hTrigSRPromptTruth->Write();
	hTrigSR->Write();
	hCorrSR->Write();
	hCorr1ptSR->Write();
	hTrigBR->Write();
	hCorrBR->Write();
	hCorr1ptBR->Write();
	hTrigSRNonDecay->Write();
	hCorrSRNonDecayMatchedJet->Write();
	hCorr1ptSRNonDecayMatchedJet->Write();
	hTrigBRNonDecay->Write();
	hCorrBRNonDecayMatchedJet->Write();
	hCorr1ptBRNonDecayMatchedJet->Write();

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

	hTrigSRPromptTruth = new THnSparseF("hTrigSRPromptTruth", "Number of true prompt clusters (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
	hTrigSR = new THnSparseF("hTrigSR", "Number of clusters (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
	hTrigBR = new THnSparseF("hTrigBR", "Number of clusters (BR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
	hTrigSRNonDecay = new THnSparseF("hTrigSRNonDecay", "Number of true non-decay clusters (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
	hTrigBRNonDecay = new THnSparseF("hTrigBRNonDecay", "Number of true non-decay clusters (BR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

	// dimensions: centrality, cluster pT, delta phi, jet pT, pT ratio, jet kT
	ndimCorr = 6;
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

	// jetkt
	nbinsCorr[5] = 120;
	minbinsCorr[5] = 0;
	maxbinsCorr[5] = 60;

	hCorrSR = new THnSparseF("hCorrSR", "Correlations (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrBR = new THnSparseF("hCorrBR", "Correlations (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorr1ptSR = new THnSparseF("hCorr1ptSR", "Correlations with 1/jetpt weight (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorr1ptBR = new THnSparseF("hCorr1ptBR", "Correlations with 1/jetpt weight (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrSR->Sumw2();
	hCorrBR->Sumw2();
	hCorr1ptSR->Sumw2();
	hCorr1ptBR->Sumw2();

	hCorrSRNonDecayMatchedJet = new THnSparseF("hCorrSRNonDecayMatchedJet", "Truth-like correlations (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrBRNonDecayMatchedJet = new THnSparseF("hCorrBRNonDecayMatchedJet", "Truth-like correlations (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorr1ptSRNonDecayMatchedJet = new THnSparseF("hCorr1ptSRNonDecayMatchedJet", "Truth-like correlations with 1/jetpt weight (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorr1ptBRNonDecayMatchedJet = new THnSparseF("hCorr1ptBRNonDecayMatchedJet", "Truth-like correlations with 1/jetpt weight (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrSRNonDecayMatchedJet->Sumw2();
	hCorrBRNonDecayMatchedJet->Sumw2();
	hCorr1ptSRNonDecayMatchedJet->Sumw2();
	hCorr1ptBRNonDecayMatchedJet->Sumw2();

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
	return false;
}

// figure out if a cluster has a non-decay photon
bool getIsClusterNonDecay(int icluster)
{
    if (cluster_is_prompt[icluster]) return true;
    
    bool foundNonDecayPhoton = false;
    for (UInt_t itruth = 0; itruth < 32; itruth++) {
        unsigned short mcindex = cluster_mc_truth_index[icluster][itruth];
        if (mcindex == 65535 || mcindex == 0) continue;
        std::cout << mcindex << ": " << mc_truth_pdg_code[mcindex] << std::endl;
        if (not(mc_truth_pdg_code[mcindex] == 22)) continue;
        int parent_pdg_code = mc_truth_first_parent_pdg_code[mcindex];
        std::cout << std::endl << "Found photon with parent " << parent_pdg_code << std::endl;
        switch (abs(parent_pdg_code) % 100000) {
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 21:
            case 22:
                foundNonDecayPhoton = true;
                break;
            default:
                return false;
        }
    }
    
    return foundNonDecayPhoton;
}

// Match truth and reco jets
void matchJetsInEvent()
{
    // start with all reco jets unmatched
	for (int ireco = 0; ireco < njet; ireco++) {
		unmatchedReco.insert(ireco);
	}

	// look at each truth jet
	for (int itruth = 0; itruth < njet_charged_truth; itruth++) {
		int matchedRecoIndex = -1;
		bool skipTruthJet = false;

		float eta_truth = jet_charged_truth_eta[itruth];
		float phi_truth = jet_charged_truth_phi[itruth];

		// look at each reco jet and find the one with the smallest distance to this truth jet
		for (int ireco = 0; ireco < njet; ireco++) {
			float eta_reco = jet_eta[ireco];
			float phi_reco = jet_phi[ireco];

			float distance2 = ((eta_truth - eta_reco) * (eta_truth - eta_reco)) + (TVector2::Phi_mpi_pi(phi_truth - phi_reco) * TVector2::Phi_mpi_pi(phi_truth - phi_reco));

			// distance threshold is 0.6 * R
			if (distance2 < (0.6 * 0.2) * (0.6 * 0.2)) {
				// if this truth jet has already found a match, skip this truth jet entirely
				if (matchedRecoIndex > -1) {
					skipTruthJet = true;
				}
				matchedRecoIndex = ireco;
			}
		}

		if (skipTruthJet) continue;

		// if a reco jet was found to match this truth jet
		if (matchedRecoIndex > -1) {
			// if that matched jet is outside the acceptance,
			// consider the truth jet to be unmatched
			if (abs(jet_eta[matchedRecoIndex]) > jet_eta_max) {
				unmatchedTruth.insert(itruth);
			} else {
				// save off this matched pair
				matchedIndex = std::make_pair(itruth, matchedRecoIndex);
				matchedJetIndices.push_back(matchedIndex);
				// remove the reco jet from the set of unmatched reco jets
				unmatchedReco.erase(matchedRecoIndex);
                matchedReco.push_back(matchedRecoIndex);
			}
		} else {
			unmatchedTruth.insert(itruth);
		}
	}
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
	} else if (isovar == "cluster_iso_its_02_sub") {
		isolation = cluster_iso_its_02[icluster] + cluster_iso_its_02_ue[icluster] - ue_estimate_its_const * M_PI * 0.2 * 0.2; 
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

float getAvgEgNtrial(std::string filename)
{
	// can't seem to figure out how to do this dynamically, so we'll just calculate it once and hard-code it here
	// it's faster anyway
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat1_18q_295913_celltrack.root") return 609.831558824;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat2_18q_296381_celltrack.root") return 39.2319029752;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat3_18q_296623_celltrack.root") return 17.1794083568;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat4_18q_296312_celltrack.root") return 14.1022386412;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat5_18q_296197_celltrack.root") return 11.9054275273;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat6_18q_296419_celltrack.root") return 9.4368308794;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat1_18q_cent1030_int7.root") return 994.049271389;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat2_18q_cent1030_int7.root") return 40.2123641696;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat3_18q_cent1030_int7.root") return 17.2359678303;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat4_18q_cent1030_int7.root") return 14.0834208953;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat5_18q_cent1030_int7.root") return 11.9398200279;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat6_18q_cent1030_int7.root") return 9.47514747474;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat1_18q_295586_celltrack.root") return 1060733.62722;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat2_18q_296509_celltrack.root") return 80066.0888592;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat3_18q_296280_celltrack.root") return 14245.3726383;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat4_18q_296065_celltrack.root") return 2724.94568333;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat5_18q_296415_celltrack.root") return 694.716856253;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat6_18q_296135_celltrack.root") return 221.5977468;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat7_18q_296511_celltrack.root") return 78.8354605923;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat8_18q_295819_celltrack.root") return 23.6664469397;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat1_18r_297193_celltrack.root") return 446021.122817;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat2_18r_297590_celltrack.root") return 83209.2395006;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat3_18r_297415_celltrack.root") return 12543.0506085;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat4_18r_297588_celltrack.root") return 2669.61427539;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat5_18r_297441_celltrack.root") return 700.724918225;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat6_18r_297218_celltrack.root") return 221.092738388;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat7_18r_297317_celltrack.root") return 78.9736794019;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat8_18r_297541_celltrack.root") return 23.8052742002;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat1_18q_cent1030_int7.root") return 1068767.82634;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat2_18q_cent1030_int7.root") return 113655.662778;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat3_18q_cent1030_int7.root") return 15863.7317306;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat4_18q_cent1030_int7.root") return 3180.62166183;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat5_18q_cent1030_int7.root") return 762.109017096;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat6_18q_cent1030_int7.root") return 235.597678758;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat7_18q_cent1030_int7.root") return 82.3022021864;
    if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3c_pthat8_18q_cent1030_int7.root") return 24.400433146;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3a/20g3a_pthat1_297588.root") return 714.190082645;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3a/20g3a_pthat2_297588.root") return 41.5663471778;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3a/20g3a_pthat3_297588.root") return 17.5137374618;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3a/20g3a_pthat4_297588.root") return 14.8231126313;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3a/20g3a_pthat5_297588.root") return 12.3522151054;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3a/20g3a_pthat6_297588.root") return 9.58740424463;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3c/20g3c_pthat1_297588.root") return 1610992.77193;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3c/20g3c_pthat2_297588.root") return 166551.245902;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3c/20g3c_pthat3_297588.root") return 21003.3947368;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3c/20g3c_pthat4_297588.root") return 3561.51378446;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3c/20g3c_pthat5_297588.root") return 813.13949382;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3c/20g3c_pthat6_297588.root") return 255.870021723;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3c/20g3c_pthat7_297588.root") return 89.7408734602;
    if (filename == "/global/project/projectdirs/alice/NTuples/MC/20g3c/20g3c_pthat8_297588.root") return 26.4355370716;

	return 0;
}

void setBranchAddresses()
{
	_tree_event->SetBranchStatus("*mc*", 0);

	// event addresses
	_tree_event->SetBranchAddress("primary_vertex", primary_vertex);
	_tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
	_tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
	_tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
	_tree_event->SetBranchAddress("centrality_v0m", &centrality_v0m);
	_tree_event->SetBranchAddress("eg_cross_section", &eg_cross_section);
    
	// MC addresses
	_tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
	_tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
	_tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
	_tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
	_tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
	_tree_event->SetBranchAddress("mc_truth_is_prompt_photon", mc_truth_is_prompt_photon);

	// track addresses
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

	// cluster addresses
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
	_tree_event->SetBranchAddress("cluster_iso_its_02", cluster_iso_its_02);
	_tree_event->SetBranchAddress("cluster_iso_its_04", cluster_iso_its_04);
	_tree_event->SetBranchAddress("cluster_frixione_tpc_04_02", cluster_frixione_tpc_04_02);
	_tree_event->SetBranchAddress("cluster_frixione_its_04_02", cluster_frixione_its_04_02);
	_tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
	_tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);

	_tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
	_tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
	_tree_event->SetBranchAddress("cell_e", cell_e);

	_tree_event->SetBranchAddress("cluster_tof", cluster_tof);
	_tree_event->SetBranchAddress("cluster_iso_its_02_ue", cluster_iso_its_02_ue);
	_tree_event->SetBranchAddress("cluster_iso_its_04_ue", cluster_iso_its_04_ue);
	_tree_event->SetBranchAddress("cluster_iso_tpc_02_ue", cluster_iso_tpc_02_ue);
	_tree_event->SetBranchAddress("cluster_iso_tpc_04_ue", cluster_iso_tpc_04_ue);

	if (auxfile != NULL) {
		auxtree->SetBranchAddress("cluster_5x5all", cluster_5x5all);
        auxtree->SetBranchAddress("cluster_is_prompt", cluster_is_prompt);
		auxtree->SetBranchAddress("isINT7", &isINT7);
	}

	// jet addresses
	// switch based on jet type
	if (jettype == "ak04tpc") {
		_tree_event->SetBranchAddress("njet_ak04tpc", &njet);
		_tree_event->SetBranchAddress("jet_ak04tpc_pt_raw", jet_pt_raw);
		_tree_event->SetBranchAddress("jet_ak04tpc_eta", jet_eta);
		_tree_event->SetBranchAddress("jet_ak04tpc_phi", jet_phi);
		_tree_event->SetBranchAddress("jet_ak04tpc_area", jet_area);
		_tree_event->SetBranchAddress("jet_ak04tpc_multiplicity_raw", jet_multiplicity_raw);

		_tree_event->SetBranchAddress("njet_charged_truth_ak04", &njet_charged_truth);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_pt", jet_charged_truth_pt);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_eta", jet_charged_truth_eta);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_phi", jet_charged_truth_phi);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_area", jet_charged_truth_area);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_multiplicity", jet_charged_truth_multiplicity);
	} else if (jettype == "ak02tpc") {
		auxtree->SetBranchAddress("njet_ak02tpc", &njet);
		auxtree->SetBranchAddress("jet_ak02tpc_pt_raw", jet_pt_raw);
		auxtree->SetBranchAddress("jet_ak02tpc_eta", jet_eta);
		auxtree->SetBranchAddress("jet_ak02tpc_phi", jet_phi);
		auxtree->SetBranchAddress("jet_ak02tpc_area", jet_area);
		auxtree->SetBranchAddress("jet_ak02tpc_multiplicity_raw", jet_multiplicity_raw);

		auxtree->SetBranchAddress("njet_charged_truth_ak02", &njet_charged_truth);
		auxtree->SetBranchAddress("jet_charged_truth_ak02_pt", jet_charged_truth_pt);
		auxtree->SetBranchAddress("jet_charged_truth_ak02_eta", jet_charged_truth_eta);
		auxtree->SetBranchAddress("jet_charged_truth_ak02_phi", jet_charged_truth_phi);
		auxtree->SetBranchAddress("jet_charged_truth_ak02_area", jet_charged_truth_area);
		auxtree->SetBranchAddress("jet_charged_truth_ak02_multiplicity", jet_charged_truth_multiplicity);
	} else if (jettype == "ak04its") {
		_tree_event->SetBranchAddress("njet_ak04its", &njet);
		_tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_pt_raw);
		_tree_event->SetBranchAddress("jet_ak04its_eta", jet_eta);
		_tree_event->SetBranchAddress("jet_ak04its_phi", jet_phi);
		_tree_event->SetBranchAddress("jet_ak04its_area", jet_area);
		_tree_event->SetBranchAddress("jet_ak04its_multiplicity_raw", jet_multiplicity_raw);

		_tree_event->SetBranchAddress("njet_charged_truth_ak04", &njet_charged_truth);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_pt", jet_charged_truth_pt);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_eta", jet_charged_truth_eta);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_phi", jet_charged_truth_phi);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_area", jet_charged_truth_area);
		_tree_event->SetBranchAddress("jet_charged_truth_ak04_multiplicity", jet_charged_truth_multiplicity);
	} else {
		std::cout << "ERROR: Jet type " << jettype << " not recognized. Aborting" << std::endl;
		exit(EXIT_FAILURE);
	}
}