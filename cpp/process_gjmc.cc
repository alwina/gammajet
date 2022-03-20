#include <TLorentzVector.h>
#include <TH2.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include "yaml-cpp/yaml.h"

#include "process_gjmc.h"
#include "config_parser.h"
#include "shared_defs.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cout << "Format: [command] [config file] [nevents (optional)]" << std::endl;
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

	initializeRooUnfoldResponses();
	initializeTHnSparses();
	Double_t trig[ndimTrig];
	Double_t corr[ndimCorr];
	Double_t photonresolution[ndimPhotonRes];
	Double_t jetresolution[ndimJetRes];
	Double_t h4[4];
	float avg_eg_ntrial = 0;
	float weight = 0;
    int nVeryHighPtJets = 0;

	/*--------------------------------------------------------------
	Loop through files
	--------------------------------------------------------------*/
	YAML::Node filenames = configrunperiod["filelists"]["ntuples"]["gjmc"];
	for (YAML::const_iterator fileit = filenames.begin(); fileit != filenames.end(); fileit++) {
		std::string root_filename = fileit->as<std::string>();
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
			if (!isINT7) continue;

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
			for (ULong64_t icluster = 0; icluster < ncluster; icluster++) {
				// apply cluster cuts
				if (rejectCluster(icluster)) continue;

				float photon_truth_pt;
				float photon_truth_eta;
				float photon_truth_phi;
				int npromptphotons = 0;
				int promptpdgsum = 0;

				// look for the prompt photon contained in this cluster
				for (UInt_t itruth = 0; itruth < 32; itruth++) {
					unsigned short mcindex = cluster_mc_truth_index[icluster][itruth];
					if (mcindex == 65535) continue;
					// std::cout << std::endl << "PDG: " << mc_truth_pdg_code[mcindex] << ", is prompt photon: " << mc_truth_is_prompt_photon[mcindex];
					if (abs(mc_truth_pdg_code[mcindex]) == 11 || mc_truth_pdg_code[mcindex] == 22) {
						if (mc_truth_is_prompt_photon[mcindex]) {
							photon_truth_pt = mc_truth_pt[mcindex];
							photon_truth_eta = mc_truth_eta[mcindex];
							photon_truth_phi = mc_truth_phi[mcindex];
							npromptphotons++;
							promptpdgsum += mc_truth_pdg_code[mcindex];
						}
					}
				}

				if (npromptphotons > 1) {
					std::cout << std::endl << "Event " << ievent << " cluster " << icluster << " has " << npromptphotons << " prompt photons (PDG sum " << promptpdgsum << ")" << std::endl;
				}

				// determine whether the cluster is isolated
				float isolation = getIsolation(icluster);
				isIsolated = GetIsIsolated(isolation, centrality_v0m, isoconfig);
				if (not(isIsolated)) continue;

				// determine whether it is SR
				float shower = getShower(icluster);
				isSignal = (shower > srmin) and (shower < srmax);
				if (not(isSignal)) continue;

				// get the centrality bin number
				int centbin = getCentBinNumber(centrality_v0m);
				int ptbin = getPtBinNumber(cluster_pt[icluster]);

				// fill photon response matrices
				photonptResponses[centbin][ptbin].Fill(cluster_pt[icluster], photon_truth_pt, weight);
				photonetaResponses[centbin][ptbin].Fill(cluster_eta[icluster], photon_truth_eta, weight);
				photonphiResponses[centbin][ptbin].Fill(cluster_phi[icluster], photon_truth_phi, weight);

				trig[0] = centrality_v0m;
				trig[1] = cluster_pt[icluster];
				hTrigSR->Fill(trig, weight);

				photonresolution[0] = centrality_v0m;
				photonresolution[1] = cluster_pt[icluster];

				jetresolution[0] = centrality_v0m;
				jetresolution[1] = cluster_pt[icluster];

				// fill photon pt and phi resolution
				photonresolution[2] = (cluster_pt[icluster] - photon_truth_pt) / photon_truth_pt;
				hPhotonPtResolution->Fill(photonresolution, weight);
				photonresolution[2] = TVector2::Phi_mpi_pi(cluster_phi[icluster] - photon_truth_phi);
				hPhotonPhiResolution->Fill(photonresolution, weight);

				/*--------------------------------------------------------------
				Loop through all reco jets to fill THnSparses
				--------------------------------------------------------------*/
				for (ULong64_t ijet = 0; ijet < njet; ijet++) {
					if (jet_pt_raw[ijet] < jet_pt_min) continue;
					if (jet_pt_raw[ijet] > jet_pt_max) continue;
					if (abs(jet_eta[ijet]) > jet_eta_max) continue;

					// Observables: delta phi, jet pT, pT ratio
					Float_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[icluster] - jet_phi[ijet]));
					Float_t jetpt = jet_pt_raw[ijet];
					Float_t ptratio = jetpt / cluster_pt[icluster];

					corr[0] = centrality_v0m;
					corr[1] = cluster_pt[icluster];
					corr[2] = deltaphi;
					corr[3] = jetpt;
					corr[4] = ptratio;
					corr[5] = cluster_pt[icluster] * sin(deltaphi);
					hCorrSRAll->Fill(corr, weight);
					hCorr1ptSRAll->Fill(corr, weight / jetpt);
				}

				/*--------------------------------------------------------------
				Loop through all truth jets to fill THnSparses
				--------------------------------------------------------------*/
				for (ULong64_t ijet = 0; ijet < njet_charged_truth; ijet++) {
                    // keep jets regardless of pT in order to estimate the kinematic efficiency
					// if (jet_charged_truth_pt[ijet] < jet_pt_min) continue;
					// if (jet_charged_truth_pt[ijet] > jet_pt_max) continue;
                    if (jet_charged_truth_pt[ijet] < 0) continue;
					if (abs(jet_charged_truth_eta[ijet]) > jet_eta_max) continue;

					// Observables: delta phi, jet pT, pT ratio
					Float_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[icluster] - jet_charged_truth_phi[ijet]));
					Float_t jetpt = jet_charged_truth_pt[ijet];
					Float_t ptratio = jetpt / cluster_pt[icluster];

					corr[0] = centrality_v0m;
					corr[1] = cluster_pt[icluster];
					corr[2] = deltaphi;
					corr[3] = jetpt;
					corr[4] = ptratio;
					corr[5] = cluster_pt[icluster] * sin(deltaphi);
					hCorrSRTruth->Fill(corr, weight);
					hCorr1ptSRTruth->Fill(corr, weight / jetpt);
				}

				/*--------------------------------------------------------------
				Loop through matched jets to fill jet response matrices
				--------------------------------------------------------------*/
				for (auto matchedIndex : matchedJetIndices) {
					int itruth = matchedIndex.first;
					int ireco = matchedIndex.second;

					float j_truth_pt = jet_charged_truth_pt[itruth];
					float j_truth_eta = jet_charged_truth_eta[itruth];
					float j_truth_phi = jet_charged_truth_phi[itruth];
					float j_truth_area = jet_charged_truth_area[itruth];
					float j_truth_mult = jet_charged_truth_multiplicity[itruth];
					float j_reco_pt = jet_pt_raw[ireco];
					float j_reco_eta = jet_eta[ireco];
					float j_reco_phi = jet_phi[ireco];
					float j_reco_area = jet_area[ireco];
					float j_reco_mult = jet_multiplicity_raw[ireco];

					// apply jet cuts on the reco jets
					if (not(j_reco_pt > jet_pt_min and j_reco_pt < jet_pt_max)) continue;
					if (not(abs(j_reco_eta) < jet_eta_max)) continue;

					float deltaphi_truth = TMath::Abs(TVector2::Phi_mpi_pi(photon_truth_phi - j_truth_phi));
					float deltaphi_reco = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[icluster] - j_reco_phi));
					// effective multiplicity: jet multiplicity scaled by pt kept / pt before UE subtraction
					float j_reco_eff_mult = j_reco_mult * (j_reco_pt) / (j_reco_pt + j_reco_area * ue_estimate_tpc_const);
                    
					corr[0] = centrality_v0m;
					corr[1] = cluster_pt[icluster];
					corr[2] = deltaphi_truth;
					corr[3] = j_truth_pt;
					corr[4] = j_truth_pt / photon_truth_pt;
					corr[5] = cluster_pt[icluster] * sin(deltaphi_truth);
					hCorrSRMatchedTruth->Fill(corr, weight);
					hCorr1ptSRMatchedTruth->Fill(corr, weight / j_truth_pt);

					corr[0] = centrality_v0m;
					corr[1] = cluster_pt[icluster];
					corr[2] = deltaphi_reco;
					corr[3] = j_reco_pt;
					corr[4] = j_reco_pt / cluster_pt[icluster];
					corr[5] = cluster_pt[icluster] * sin(deltaphi_reco);
					hCorrSRMatchedReco->Fill(corr, weight);
					hCorr1ptSRMatchedReco->Fill(corr, weight / j_reco_pt);

					if (deltaphi_truth > 7 * M_PI / 8 && deltaphi_reco > 7 * M_PI / 8) {
						ptratiojetptResponses[centbin][ptbin].Fill(j_reco_pt / cluster_pt[icluster], j_reco_pt, j_truth_pt / photon_truth_pt, j_truth_pt, weight);
						ptratioResponses[centbin][ptbin].Fill(j_reco_pt / cluster_pt[icluster], j_truth_pt / photon_truth_pt, weight);
						jetptResponses[centbin][ptbin].Fill(j_reco_pt, j_truth_pt, weight);
						jetetaResponses[centbin][ptbin].Fill(j_reco_eta, j_truth_eta, weight);
						jetphiResponses[centbin][ptbin].Fill(j_reco_phi, j_truth_phi, weight);
						jetareaResponses[centbin][ptbin].Fill(j_reco_area, j_truth_area, weight);
						jetmultResponses[centbin][ptbin].Fill(j_reco_mult, j_truth_mult, weight);
						jeteffmultResponses[centbin][ptbin].Fill(j_reco_eff_mult, j_truth_mult, weight);

						h4[0] = j_reco_pt / cluster_pt[icluster];
						h4[1] = j_truth_pt / photon_truth_pt;
						h4[2] = j_reco_pt;
						h4[3] = j_truth_pt;
						ptratiojetptHists[centbin][ptbin]->Fill(h4, weight);

						// fill jet pt and phi resolution
						jetresolution[2] = j_reco_pt;
						jetresolution[3] = (j_reco_pt - j_truth_pt) / j_truth_pt;
						hJetB2bPtResolution->Fill(jetresolution, weight);

						jetresolution[2] = j_reco_pt;
						jetresolution[3] = TVector2::Phi_mpi_pi(j_reco_phi - j_truth_phi);
						hJetB2bPhiResolution->Fill(jetresolution, weight);
					}
                    
                    // deltaphi requires ptratio < 1.2, which we do at the reco level
                    // it shouldn't really matter because there are very few jets that get cut, but just in case
                    // unclear whether this should also happen at the truth level though
                    if (j_reco_pt / cluster_pt[icluster] < 1.2) {
                        deltaphijetptResponses[centbin][ptbin].Fill(deltaphi_reco, j_reco_pt, deltaphi_truth, j_truth_pt, weight);
                        deltaphiResponses[centbin][ptbin].Fill(deltaphi_reco, deltaphi_truth, weight);

                        h4[0] = deltaphi_reco;
                        h4[1] = deltaphi_truth;
                        h4[2] = j_reco_pt;
                        h4[3] = j_truth_pt;
                        deltaphijetptHists[centbin][ptbin]->Fill(h4, weight);

                        if (j_truth_pt / photon_truth_pt > 1.2) {
                            nVeryHighPtJets++;
                        }
                    }
                    
					// fill jet pt and phi resolution
					jetresolution[2] = j_reco_pt;
					jetresolution[3] = (j_reco_pt - j_truth_pt) / j_truth_pt;
					hJetPtResolution->Fill(jetresolution, weight);

					jetresolution[2] = j_reco_pt;
					jetresolution[3] = j_reco_phi - j_truth_phi;
					hJetPhiResolution->Fill(jetresolution, weight);
				}

				if (keepMisses) {
					for (int itruth : unmatchedTruth) {
						deltaphiResponses[centbin][ptbin].Miss(TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[icluster] - jet_charged_truth_phi[itruth]))), weight;
						jetptResponses[centbin][ptbin].Miss(jet_charged_truth_pt[itruth], weight);
						ptratioResponses[centbin][ptbin].Miss(jet_charged_truth_pt[itruth] / cluster_pt[icluster], weight);
					}
				}

				if (keepFakes) {
					for (int ireco : unmatchedReco) {
						deltaphiResponses[centbin][ptbin].Fake(TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[icluster] - jet_phi[ireco])), weight);
						jetptResponses[centbin][ptbin].Fake(jet_pt_raw[ireco], weight);
						ptratioResponses[centbin][ptbin].Fake(jet_pt_raw[ireco] / cluster_pt[icluster], weight);
					}
				}
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
	TFile* foutrm;
	foutrm = new TFile((TString) configrunperiod["filelists"]["correlations"]["responsematrix"].as<std::string>(), "RECREATE");
    std::cout << "Found " << nVeryHighPtJets << " truth jets with pT > 1.2 * photon pT" << std::endl;
	std::cout << "Writing response matrices to file" << std::endl;

	for (int i = 0; i < ncentralityranges; i++) {
		for (int j = 0; j < nphotonptranges; j++) {
			std::string deltaphijetptname = "deltaphijetptResponse" + std::to_string(i) + std::to_string(j);
			std::string ptratiojetptname = "ptratiojetptResponse" + std::to_string(i) + std::to_string(j);
			std::string deltaphijetpthistname = "deltaphijetptHist" + std::to_string(i) + std::to_string(j);
			std::string ptratiojetpthistname = "ptratiojetptHist" + std::to_string(i) + std::to_string(j);
			std::string deltaphiname = "deltaphiResponse" + std::to_string(i) + std::to_string(j);
			std::string ptrationame = "ptratioResponse" + std::to_string(i) + std::to_string(j);
			std::string jetptname = "jetptResponse" + std::to_string(i) + std::to_string(j);
			std::string jetetaname = "jetetaResponse" + std::to_string(i) + std::to_string(j);
			std::string jetphiname = "jetphiResponse" + std::to_string(i) + std::to_string(j);
			std::string jetareaname = "jetareaResponse" + std::to_string(i) + std::to_string(j);
			std::string jetmultname = "jetmultResponse" + std::to_string(i) + std::to_string(j);
			std::string jeteffmultname = "jeteffmultResponse" + std::to_string(i) + std::to_string(j);
			std::string photonptname = "photonptResponse" + std::to_string(i) + std::to_string(j);
			std::string photonetaname = "photonetaResponse" + std::to_string(i) + std::to_string(j);
			std::string photonphiname = "photonphiResponse" + std::to_string(i) + std::to_string(j);

			deltaphijetptResponses[i][j].Write(deltaphijetptname.c_str());
			ptratiojetptResponses[i][j].Write(ptratiojetptname.c_str());
			deltaphijetptHists[i][j]->Write(deltaphijetpthistname.c_str());
			ptratiojetptHists[i][j]->Write(ptratiojetpthistname.c_str());
			deltaphiResponses[i][j].Write(deltaphiname.c_str());
			ptratioResponses[i][j].Write(ptrationame.c_str());
			jetptResponses[i][j].Write(jetptname.c_str());
			jetetaResponses[i][j].Write(jetetaname.c_str());
			jetphiResponses[i][j].Write(jetphiname.c_str());
			jetareaResponses[i][j].Write(jetareaname.c_str());
			jetmultResponses[i][j].Write(jetmultname.c_str());
			jeteffmultResponses[i][j].Write(jeteffmultname.c_str());
			photonptResponses[i][j].Write(photonptname.c_str());
			photonetaResponses[i][j].Write(photonetaname.c_str());
			photonphiResponses[i][j].Write(photonphiname.c_str());
		}
	}

	foutrm->Close();

	TFile* fouthists;
	fouthists = new TFile((TString) configrunperiod["filelists"]["correlations"]["gjmcchecks"].as<std::string>(), "RECREATE");
	std::cout << "Writing THnSparses to file" << std::endl;

	hTrigSR->Write();
	hCorrSRTruth->Write();
	hCorrSRAll->Write();
    hCorrSRMatchedTruth->Write();
    hCorrSRMatchedReco->Write();
	hCorr1ptSRTruth->Write();
	hCorr1ptSRAll->Write();
    hCorr1ptSRMatchedTruth->Write();
    hCorr1ptSRMatchedReco->Write();
    hPhotonPtResolution->Write();
	hPhotonPhiResolution->Write();
	hJetPtResolution->Write();
	hJetPhiResolution->Write();
	hJetB2bPtResolution->Write();
	hJetB2bPhiResolution->Write();

	fouthists->Close();

	std::cout << "Ending" << std::endl;
	return EXIT_SUCCESS;
}

void printCutSummary()
{
	std::cout << "Cluster pT range: " << cluster_pt_min << "-" << cluster_pt_max << std::endl;
	std::cout << "Cluster eta max: " << cluster_eta_max << std::endl;
	std::cout << "Cluster ncell min: " << cluster_ncell_min << std::endl;
	std::cout << "Cluster Ecross/Emax min: " << cluster_ecross_emax_min << std::endl;
	std::cout << "Cluster dist to bad channel min: " << cluster_dbc_min << std::endl;
	std::cout << "Cluster nlocal maxima max: " << cluster_nlm_max << std::endl;
	std::cout << "Shower shape SR range: " << srmin << "-" << srmax << std::endl;
	std::cout << "Shower shape BR range: " << brmin << "-" << brmax << std::endl;
	std::cout << "Jet type: " << jettype << std::endl;
	std::cout << "Jet pT range: " << jet_pt_min << "-" << jet_pt_max << std::endl;
	std::cout << "Jet eta max: " << jet_eta_max << std::endl;
}

void initializeRooUnfoldResponses()
{
	TH2F* h2DeltaphiJetpt = new TH2F("", "", deltaphi_nbins, deltaphi_min, deltaphi_max, jetpt_nbins, jetpt_min, jetpt_max);
	TH2F* h2PtratioJetpt = new TH2F("", "", ptratio_nbins, ptratio_min, ptratio_max, jetpt_nbins, jetpt_min, jetpt_max);

	Int_t nbinsh4dpjp[4] = {120, 120, 120, 120};
	Double_t minbinsh4dpjp[4] = {deltaphi_min, deltaphi_min, jetpt_min, jetpt_min};
	Double_t maxbinsh4dpjp[4] = {deltaphi_max, deltaphi_max, jetpt_max, jetpt_max};

	Int_t nbinsh4prjp[4] = {120, 120, 120, 120};
	Double_t minbinsh4prjp[4] = {ptratio_min, ptratio_min, jetpt_min, jetpt_min};
	Double_t maxbinsh4prjp[4] = {ptratio_max, ptratio_max, jetpt_max, jetpt_max};

	for (int i = 0; i < ncentralityranges; i++) {
		deltaphijetptResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		ptratiojetptResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		deltaphijetptHists.push_back(std::vector<THnSparseF*> (nphotonptranges));
		ptratiojetptHists.push_back(std::vector<THnSparseF*> (nphotonptranges));
		deltaphiResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		ptratioResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		jetptResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		jetetaResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		jetphiResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		jetareaResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		jetmultResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		jeteffmultResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		photonptResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		photonetaResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));
		photonphiResponses.push_back(std::vector<RooUnfoldResponse> (nphotonptranges));

		for (int j = 0; j < nphotonptranges; j++) {
			THnSparseF* h4DeltaphiJetpt = new THnSparseF("", "", 4, nbinsh4dpjp, minbinsh4dpjp, maxbinsh4dpjp);
			h4DeltaphiJetpt->Sumw2();
			deltaphijetptHists[i][j] = h4DeltaphiJetpt;

			THnSparseF* h4PtratioJetpt = new THnSparseF("", "", 4, nbinsh4prjp, minbinsh4prjp, maxbinsh4prjp);
			h4PtratioJetpt->Sumw2();
			ptratiojetptHists[i][j] = h4PtratioJetpt;

			RooUnfoldResponse deltaphijetptResponse(h2DeltaphiJetpt, h2DeltaphiJetpt);
			deltaphijetptResponses[i][j] = deltaphijetptResponse;

			RooUnfoldResponse ptratiojetptResponse(h2PtratioJetpt, h2PtratioJetpt);
			ptratiojetptResponses[i][j] = ptratiojetptResponse;

			RooUnfoldResponse deltaphiResponse(120, deltaphi_min, deltaphi_max);
			deltaphiResponses[i][j] = deltaphiResponse;

			RooUnfoldResponse ptratioResponse(200, ptratio_min, ptratio_max);
			ptratioResponses[i][j] = ptratioResponse;

			RooUnfoldResponse jetptResponse(120, jetpt_min, jetpt_max);
			jetptResponses[i][j] = jetptResponse;

			RooUnfoldResponse jetetaResponse(140, -0.7, 0.7);
			jetetaResponses[i][j] = jetetaResponse;

			RooUnfoldResponse jetphiResponse(120, -M_PI, M_PI);
			jetphiResponses[i][j] = jetphiResponse;

			RooUnfoldResponse jetareaResponse(40, 0.0, 0.4);
			jetareaResponses[i][j] = jetareaResponse;

			RooUnfoldResponse jetmultResponse(200, 0.0, 100.0);
			jetmultResponses[i][j] = jetmultResponse;

			RooUnfoldResponse jeteffmultResponse(200, 0.0, 100.0);
			jeteffmultResponses[i][j] = jeteffmultResponse;

			RooUnfoldResponse photonptResponse(10 * (cluster_pt_max - cluster_pt_min), cluster_pt_min, cluster_pt_max);
			photonptResponses[i][j] = photonptResponse;

			RooUnfoldResponse photonetaResponse(140, -0.7, 0.7);
			photonetaResponses[i][j] = photonetaResponse;

			RooUnfoldResponse photonphiResponse(120, -M_PI, M_PI);
			photonphiResponses[i][j] = photonphiResponse;
		}
	}
}

/*--------------------------------------------------------------
Set up THnSparses
hTrigSR: Number of real photon clusters in SR
hCorrSRTruth: SR photons correlated with truth jets
hCorrSRAll: SR photons correlated with all jets
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

	hTrigSR = new THnSparseF("hTrigSR", "Number of real photon clusters in SR", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

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

	hCorrSRTruth = new THnSparseF("hCorrSRTruth", "Correlations (Truth jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrSRAll = new THnSparseF("hCorrSRAll", "Correlations (All jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
    hCorrSRMatchedTruth = new THnSparseF("hCorrSRMatchedTruth", "Correlations (Matched truth jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
    hCorrSRMatchedReco = new THnSparseF("hCorrSRMatchedReco", "Correlations (Matched reco jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorr1ptSRTruth = new THnSparseF("hCorr1ptSRTruth", "Correlations with 1/jetpt weight (Truth jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorr1ptSRAll = new THnSparseF("hCorr1ptSRAll", "Correlations with 1/jetpt weight (All jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
    hCorr1ptSRMatchedTruth = new THnSparseF("hCorr1ptSRMatchedTruth", "Correlations with 1/jetpt weight (Matched truth jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
    hCorr1ptSRMatchedReco = new THnSparseF("hCorr1ptSRMatchedReco", "Correlations with 1/jetpt weight (Matched reco jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrSRTruth->Sumw2();
	hCorrSRAll->Sumw2();
	hCorr1ptSRTruth->Sumw2();
	hCorr1ptSRAll->Sumw2();

	// resolutions
	ndimPhotonRes = 3;
	ndimJetRes = 4;

	Int_t nbinsPhotonPtResolution[ndimPhotonRes] = {10, nbinsClusterPt, 100};
	Double_t minbinsPhotonPtResolution[ndimPhotonRes] = {0, cluster_pt_min, -0.5};
	Double_t maxbinsPhotonPtResolution[ndimPhotonRes] = {100, cluster_pt_max, 0.5};
	hPhotonPtResolution = new THnSparseF("hPhotonPtResolution", "(reco - truth) / truth", ndimPhotonRes, nbinsPhotonPtResolution, minbinsPhotonPtResolution, maxbinsPhotonPtResolution);
	hPhotonPtResolution->Sumw2();

	Int_t nbinsJetPtResolution[ndimJetRes] = {10, nbinsClusterPt, 120, 400};
	Double_t minbinsJetPtResolution[ndimJetRes] = {0, cluster_pt_min, 0, -2.0};
	Double_t maxbinsJetPtResolution[ndimJetRes] = {100, cluster_pt_max, 60, 2.0};
	hJetPtResolution = new THnSparseF("hJetPtResolution", "(reco - truth) / truth", ndimJetRes, nbinsJetPtResolution, minbinsJetPtResolution, maxbinsJetPtResolution);
	hJetPtResolution->Sumw2();
	hJetB2bPtResolution = new THnSparseF("hJetB2bPtResolution", "(reco - truth) / truth", ndimJetRes, nbinsJetPtResolution, minbinsJetPtResolution, maxbinsJetPtResolution);
	hJetB2bPtResolution->Sumw2();

	Int_t nbinsPhotonPhiResolution[ndimPhotonRes] = {10, nbinsClusterPt, 120};
	Double_t minbinsPhotonPhiResolution[ndimPhotonRes] = {0, cluster_pt_min, -M_PI};
	Double_t maxbinsPhotonPhiResolution[ndimPhotonRes] = {100, cluster_pt_max, M_PI};
	hPhotonPhiResolution = new THnSparseF("hPhotonPhiResolution", "reco - truth", ndimPhotonRes, nbinsPhotonPhiResolution, minbinsPhotonPhiResolution, maxbinsPhotonPhiResolution);
	hPhotonPhiResolution->Sumw2();

	Int_t nbinsJetPhiResolution[ndimJetRes] = {10, nbinsClusterPt, 120, 120};
	Double_t minbinsJetPhiResolution[ndimJetRes] = {0, cluster_pt_min, 0, -M_PI};
	Double_t maxbinsJetPhiResolution[ndimJetRes] = {100, cluster_pt_max, 60, M_PI};
	hJetPhiResolution = new THnSparseF("hJetPhiResolution", "reco - truth", ndimJetRes, nbinsJetPhiResolution, minbinsJetPhiResolution, maxbinsJetPhiResolution);
	hJetPhiResolution->Sumw2();
	hJetB2bPhiResolution = new THnSparseF("hJetB2bPhiResolution", "reco - truth", ndimJetRes, nbinsJetPhiResolution, minbinsJetPhiResolution, maxbinsJetPhiResolution);
	hJetB2bPhiResolution->Sumw2();
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
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat1_18q_295913_celltrack.root") return 609.83156;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat2_18q_296381_celltrack.root") return 39.231903;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat3_18q_296623_celltrack.root") return 17.179408;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat4_18q_296312_celltrack.root") return 14.102239;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat5_18q_296197_celltrack.root") return 11.905428;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat6_18q_296419_celltrack.root") return 9.4368309;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat1_18q_cent1030_int7.root") return 994.04927;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat2_18q_cent1030_int7.root") return 40.212364;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat3_18q_cent1030_int7.root") return 17.235968;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat4_18q_cent1030_int7.root") return 14.083421;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat5_18q_cent1030_int7.root") return 11.939820;
	if (filename == "/global/project/projectdirs/alice/NTuples/embed/embed_20g3a_pthat6_18q_cent1030_int7.root") return 9.4751475;
	return 0;
}

void setBranchAddresses()
{
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
	_tree_event->SetBranchAddress("cluster_nmc_truth", cluster_nmc_truth);
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

// Apply cluster cuts
bool rejectCluster(int icluster)
{
	if (not(cluster_pt[icluster] > cluster_pt_min and cluster_pt[icluster] < cluster_pt_max)) return true;
	if (not(abs(cluster_eta[icluster]) < cluster_eta_max)) return true;
	if (not(cluster_ncell[icluster] >= cluster_ncell_min)) return true;
	if (not(cluster_e_cross[icluster] / cluster_e_max[icluster] > cluster_ecross_emax_min)) return true;
	if (not(cluster_distance_to_bad_channel[icluster] >= cluster_dbc_min)) return true;
	if (not(cluster_nlocal_maxima[icluster] < cluster_nlm_max)) return true;
	if (not(cluster_is_prompt[icluster])) return true;

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

// need to make a response matrix per centrality
int getCentBinNumber(float centrality)
{
	int binNumber = -1;

	for (int icent = 0; icent < ncentralityranges; icent++) {
		std::pair<float, float> centrange = centralityranges[icent].as<std::pair<float, float>>();
		if (centrality >= centrange.first && centrality < centrange.second) {
			binNumber = icent;
			break;
		}
	}

	return binNumber;
}

// also a response matrix per photon pt
int getPtBinNumber(float pt)
{
	int binNumber = -1;

	for (int ipt = 0; ipt < nphotonptranges; ipt++) {
		std::pair<float, float> ptrange = photonptranges[ipt].as<std::pair<float, float>>();
		if (pt >= ptrange.first && pt < ptrange.second) {
			binNumber = ipt;
			break;
		}
	}

	return binNumber;
}