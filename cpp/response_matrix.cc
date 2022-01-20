#include <iostream>
#include <math.h>

#include "response_matrix.h"
#include "config_parser.h"
#include "shared_defs.h"

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout << "Format: [command] [config file] [nevents (optional)]" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (argc > 2) {
    nevents_max = std::stol(argv[2])
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

  /*--------------------------------------------------------------
  Loop through files
  --------------------------------------------------------------*/
  YAML::Node filenames = configrunperiod["filelists"]["ntuples"]["gjmc"];
  for (YAML::const_iterator it = filenames.begin(); it != filenames.end(); it++) {
    std::string root_filename = fileit->as<std::string>();
    openFilesAndGetTTrees(root_filename);
    setBranchAddresses();

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

      matchJetsInEvent();

      /*--------------------------------------------------------------
      Loop through clusters
      --------------------------------------------------------------*/
      for (ULong64_t icluster = 0; icluster < ncluster; icluster++) {
        // apply cluster cuts
        if (rejectCluster(icluster)) continue;

        // determine whether the cluster is isolated
        float isolation = getIsolation(icluster);
        isIsolated = GetIsIsolated(isolation, centrality_v0m, isoconfig);
        if (not(isIsolated)) continue;

        // determine whether it is SR
        float shower = getShower(icluster);
        isSignal = (shower > srmin) and (shower < srmax);
        if (not(isSignal)) continue;

        // get the centrality bin number
        int centbin = getCentBinNumber(centrality_v0m, centralityranges);

        /*--------------------------------------------------------------
        Loop through matched jets to fill response matrices
        --------------------------------------------------------------*/
        for (auto matchedIndex : matchedJetIndices) {
          int itruth = matchedIndex.first;
          int ireco = matchedIndex.second;

          float j_truth_pt = jet_charged_truth_ak04_pt[itruth];
          float j_truth_phi = jet_charged_truth_ak04_phi[itruth];
          float j_reco_pt = jet_ak04tpc_pt_raw[ireco];
          float j_reco_phi = jet_ak04tpc_phi[ireco];

          float deltaphi_truth = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - j_truth_phi));
          float deltaphi_reco = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - j_reco_phi));

          deltaphiResponses[centbin].Fill(deltaphi_reco, deltaphi_truth);
          jetptResponses[centbin].Fill(j_reco_pt, j_truth_pt);
          ptratioResponses[centbin].Fill(j_reco_pt / cluster_pt[n], j_truth_pt / cluster_pt[n]);
        }

        if (keepMisses) {
          for (int itruth : unmatchedTruth) {
            deltaphiResponses[centbin].Miss(TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - jet_charged_truth_ak04_phi[itruth])));
            jetptResponses[centbin].Miss(jet_charged_truth_ak04_pt[itruth]);
            ptratioResponses[centbin].Miss(jet_charged_truth_ak04_pt[itruth] / cluster_pt[n]);
          }
        }

        if (keepFakes) {
          for (int ireco : unmatchedReco) {
            deltaphiResponses[centbin].Fake(TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - jet_ak04tpc_phi[ireco])));
            jetptResponses[centbin].Fake(jet_ak04tpc_pt_raw[ireco]);
            ptratioResponses[centbin].Fake(jet_ak04tpc_pt_raw[ireco] / cluster_pt[n]);
          }
        }
      } // end cluster loop

      matchedJetIndices.clear();
      unmatchedTruth.clear();
      unmatchedReco.clear();
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
  TFile* fout;
  fout = new TFile((TString) configrunperiod["filelists"]["correlations"]["responsematrix"].as<std::string>(), "RECREATE");
  std::cout << "Writing to file" << std::endl;

  for (int i = 0; i < ncentralityranges; i++) {
    std::string dname = "deltaphiResponse" + std::to_string(i);
    std::string jname = "jetptResponse" + std::to_string(i);
    std::string pname = "ptratioResponse" + std::to_string(i);

    deltaphiResponses[i].Write(dname.c_str());
    jetptResponses[i].Write(jname.c_str());
    ptratioResponses[i].Write(pname.c_str());
  }

  fout->Close();
  std::cout << " ending " << std::endl;
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
  std::cout << "Cluster TOF max: " << cluster_tof_max << std::endl;
  std::cout << "Shower shape SR range: " << srmin << "-" << srmax << std::endl;
  std::cout << "Shower shape BR range: " << brmin << "-" << brmax << std::endl;
  std::cout << "Jet type: " << jettype << std::endl;
  std::cout << "Jet pT range: " << jet_pt_min << "-" << jet_pt_max << std::endl;
  std::cout << "Jet eta max: " << jet_eta_max << std::endl;
}

void initializeRooUnfoldResponses()
{
  std::vector<RooUnfoldResponse> deltaphiResponses;
  std::vector<RooUnfoldResponse> jetptResponses;
  std::vector<RooUnfoldResponse> ptratioResponses;

  for (int i = 0; i < ncentralityranges; i++) {
    RooUnfoldResponse deltaphiResponse(10, 0, M_PI, 10, 0, M_PI);
    RooUnfoldResponse jetptResponse(9, 5, 50, 9, 5, 50);
    RooUnfoldResponse ptratioResponse(10, 0, 2, 10, 0, 2);

    deltaphiResponses.push_back(deltaphiResponse);
    jetptResponses.push_back(jetptResponse);
    ptratioResponses.push_back(ptratioResponse);
  }
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

    _tree_event->SetBranchAddress("njet_charged_truth_ak04", &njet_charged_truth);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_pt", jet_charged_truth_pt);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_eta", jet_charged_truth_eta);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_phi", jet_charged_truth_phi);
  } else if (jettype == "ak02tpc") {
    auxtree->SetBranchAddress("njet_ak02tpc", &njet);
    auxtree->SetBranchAddress("jet_ak02tpc_pt_raw", jet_pt_raw);
    auxtree->SetBranchAddress("jet_ak02tpc_eta", jet_eta);
    auxtree->SetBranchAddress("jet_ak02tpc_phi", jet_phi);

    auxtree->SetBranchAddress("njet_charged_truth_ak02", &njet_charged_truth);
    auxtree->SetBranchAddress("jet_charged_truth_ak02_pt", jet_charged_truth_pt);
    auxtree->SetBranchAddress("jet_charged_truth_ak02_eta", jet_charged_truth_eta);
    auxtree->SetBranchAddress("jet_charged_truth_ak02_phi", jet_charged_truth_phi);
  } else if (jettype == "ak04its") {
    _tree_event->SetBranchAddress("njet_ak04its", &njet);
    _tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_pt_raw);
    _tree_event->SetBranchAddress("jet_ak04its_eta", jet_eta);
    _tree_event->SetBranchAddress("jet_ak04its_phi", jet_phi);

    _tree_event->SetBranchAddress("njet_charged_truth_ak04", &njet_charged_truth);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_pt", jet_charged_truth_pt);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_eta", jet_charged_truth_eta);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_phi", jet_charged_truth_phi);
  } else {
    std::cout << "ERROR: Jet type " << jettype << " not recognized. Aborting" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void matchJetsInEvent()
{
  // first, make all reco jets unmatched
  for (int ireco = 0; ireco < njet; ireco++) {
    unmatchedReco.insert(ireco);
  }

  // look at each truth jet
  for (int itruth = 0; itruth < njet_charged_truth; itruth++) {
    float minDistance2 = 99999999999;
    int minDistanceIndex = -1;

    float eta_truth = jet_charged_truth_eta[itruth];
    float phi_truth = jet_charged_truth_phi[itruth];

    // look at each reco jet and find the one with the smallest distance to this truth jet
    for (int ireco = 0; ireco < njet; ireco++) {
      float eta_reco = jet_eta[ireco];
      float phi_reco = jet_phi[ireco];

      float distance2 = ((eta_truth - eta_reco) * (eta_truth - eta_reco)) + ((phi_truth - phi_reco) * (phi_truth - phi_reco));

      if (distance2 < minDistance2) {
        minDistance2 = distance2;
        minDistanceIndex = ireco;
      }
    }

    // if a reco jet was found to match this truth jet
    if (minDistanceIndex > -1) {
      // if that matched jet is outside the acceptance,
      // consider the truth jet to be unmatched
      if (abs(jet_eta[minDistanceIndex]) > jet_eta_max) {
        unmatchedTruth.insert(itruth);
      } else {
        // save off this matched pair
        matchedIndex = make_pair(itruth, minDistanceIndex);
        matchedJetIndices.push_back(matchedIndex);
        // remove the reco jet from the set of unmatched reco jets
        unmatchedReco.erase(minDistanceIndex);
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

  // skip TOF cut because this is MC, not data
  // for the gamma-jet MC, we want to select prompt clusters only

  bool isPrompt = false;
  for (UInt_t itruth = 0; itruth < cluster_nmc_truth[icluster]; itruth++) {
    unsigned short mcindex = cluster_mc_truth_index[icluster][itruth];
    if (mcindex == 65535) continue;
    if (!(abs(mc_truth_pdg_code[itruth]) == 11 || mc_truth_pdg_code[itruth] == 22)) continue;
    if (mc_truth_is_prompt_photon[itruth]) {
      isPrompt = true;
      break;
    }
  }

  if (not(isPrompt)) return true;

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
int getCentBinNumber(float centrality, YAML::Node centralityranges)
{
  int ncentralityranges = centralityranges.size();
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
