#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <THStack.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <iostream>
#include <fstream>

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

#include "RooUnfoldResponse.h"
#include "yaml-cpp/yaml.h"
#include "shared_defs.h"

using namespace std;

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_ISO_ITS_04_SUB, CLUSTER_ISO_TPC_02_SUB, CLUSTER_ISO_TPC_04_SUB, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};


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


int main(int argc, char *argv[])
{
  if (argc < 1) {
    fprintf(stderr, "Format: [command] [config file]\n");
    exit(EXIT_FAILURE);
  }

  // load config files
  // each config points to the next
  std::vector<YAML::Node> allconfigs;
  YAML::Node configrunperiod = YAML::LoadFile(argv[1]);
  allconfigs.push_back(configrunperiod);
  YAML::Node configsystem = YAML::LoadFile(configrunperiod["systemconfig"].as<std::string>());
  allconfigs.push_back(configsystem);
  YAML::Node configglobal = YAML::LoadFile(configsystem["globalconfig"].as<std::string>());
  allconfigs.push_back(configglobal);

  double srmin = 0;
  double srmax = 0;
  double brmin = 0;
  double brmax = 0;
  double pT_min = 0;
  double pT_max = 0;
  double Eta_max = 0;
  double Cluster_min = 0;
  float Cluster_DtoBad = 0;
  UChar_t Cluster_NLocal_Max = 0;
  double EcrossoverE_min = 0;

  double jet_pt_min = 0;
  double jet_pt_max = 500;
  double jet_eta_max = 10;

  bool do_pile = false;

  //double deta_max = 0;
  isolationDet determiner = CLUSTER_ISO_ITS_04;
  std::string shower_shape = "DNN";
  std::string purity_deviation = "None";

  bool TPC_Iso_Flag = false;

  YAML::Node purityconfig;
  YAML::Node isoconfig;
  YAML::Node centralityranges;

  int ncentralityranges;
  bool keepFakes;
  bool keepMisses;

  // go through the configs backwards; that way more specific settings
  // can override those from the more general configs
  for (auto it = allconfigs.rbegin(); it != allconfigs.rend(); ++it)
  {
    YAML::Node config = *it;

    // parse config file
    // check for existence first, then cast as appropriate
    if (config["showershape"]) {
      srmin = config["showershape"]["srmin"].as<double>();
      srmax = config["showershape"]["srmax"].as<double>();
      brmin = config["showershape"]["brmin"].as<double>();
      brmax = config["showershape"]["brmax"].as<double>();

      shower_shape = config["showershape"]["ssvar"].as<std::string>();
      std::cout << "Shower Shape: " << shower_shape << std::endl;
    }

    if (config["clustercuts"]["all"]["cluster_pt"]) {
      pT_min = config["clustercuts"]["all"]["cluster_pt"]["min"].as<double>();
      pT_max = config["clustercuts"]["all"]["cluster_pt"]["max"].as<double>();
    }

    if (config["clustercuts"]["all"]["cluster_eta"]) {
      Eta_max = config["clustercuts"]["all"]["cluster_eta"]["max"].as<double>();
    }

    if (config["clustercuts"]["all"]["cluster_ncell"]) {
      Cluster_min = config["clustercuts"]["all"]["cluster_ncell"]["incmin"].as<double>();
    }

    if (config["clustercuts"]["all"]["cluster_distance_to_bad_channel"]) {
      Cluster_DtoBad = config["clustercuts"]["all"]["cluster_distance_to_bad_channel"]["incmin"].as<double>();
    }

    if (config["clustercuts"]["all"]["cluster_nlocal_maxima"]) {
      Cluster_NLocal_Max = config["clustercuts"]["all"]["cluster_nlocal_maxima"]["max"].as<double>();
    }

    if (config["clustercuts"]["all"]["cluster_ecross_emax"]) {
      EcrossoverE_min = config["clustercuts"]["all"]["cluster_ecross_emax"]["min"].as<double>();
    }

    if (config["jetcuts"]) {
      jet_pt_min = config["jetcuts"]["jet_ak04tpc_pt_raw"]["min"].as<double>();
      jet_pt_max = config["jetcuts"]["jet_ak04tpc_pt_raw"]["max"].as<double>();
      jet_eta_max = config["jetcuts"]["jet_ak04tpc_eta"]["max"].as<double>();
    }

    if (config["do_pileup_cut"]) {
      do_pile = config["do_pileup_cut"].as<bool>();
    }

    if (config["isolation"]) {
      isoconfig = config["isolation"];
      std::string determinant = config["isolation"]["isovar"].as<std::string>();

      if (determinant == "cluster_iso_tpc_04") {
        determiner = CLUSTER_ISO_TPC_04;
        std::cout << "Isolation Variable: cluster_iso_tpc_04" << std::endl;
      }

      else if (determinant == "cluster_iso_its_04") {
        determiner = CLUSTER_ISO_ITS_04;
        std::cout << "Isolation Variable: cluster_iso_its_04" << std::endl;
      }

      else if (determinant == "cluster_iso_its_04_sub") {
        determiner = CLUSTER_ISO_ITS_04_SUB;
        std::cout << "Isolation Variable: cluster_iso_its_04_sub" << std::endl;
      }

      else if (determinant == "cluster_iso_tpc_02_sub") {
        determiner = CLUSTER_ISO_TPC_02_SUB;
        TPC_Iso_Flag = true;
        std::cout << "Isolation Variable: cluster_iso_tpc_02_sub" << std::endl;
      }

      else if (determinant == "cluster_iso_tpc_04_sub") {
        determiner = CLUSTER_ISO_TPC_04_SUB;
        TPC_Iso_Flag = true;
        std::cout << "Isolation Variable: cluster_iso_tpc_04_sub" << std::endl;
      }

      else if (determinant == "cluster_frixione_tpc_04_02") {
        determiner = CLUSTER_FRIXIONE_TPC_04_02;
        std::cout << "Isolation Variable: cluster_frixione_tpc_04_02" << std::endl;
      }

      else if (determinant == "cluster_frixione_its_04_02") {
        determiner = CLUSTER_FRIXIONE_ITS_04_02;
        std::cout << "Isolation Variable: cluster_frixione_its_04_02" << std::endl;
      }

      else {
        std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    if (config["Purity_Dev"]) {
      purity_deviation = config["Purity_Dev"].as<std::string>();
      std::cout << "Purity Deviation Change: " << purity_deviation << std::endl;
    }

    if (config["purity"]) {
      purityconfig = config["purity"];
    }

    if (config["responsematrix"]) {
      keepMisses = config["responsematrix"]["keepmisses"].as<bool>();
      keepFakes = config["responsematrix"]["keepfakes"].as<bool>();
    }

    if (config["centralityranges"]) {
      centralityranges = config["centralityranges"];
      ncentralityranges = centralityranges.size();
    }
  }

  /*--------------------------------------------------------------
  Setting up RooUnfoldResponse objects
  --------------------------------------------------------------*/
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

  /*--------------------------------------------------------------
  Setting up local variables to be linked with ROOT branches
  --------------------------------------------------------------*/
  //Events
  Bool_t is_pileup_from_spd_5_08;
  Double_t primary_vertex[3];
  Float_t ue_estimate_its_const;
  Float_t ue_estimate_tpc_const;
  Float_t centrality_v0m;

  //Tracks
  UInt_t ntrack;
  Float_t track_e[NTRACK_MAX];
  Float_t track_pt[NTRACK_MAX];
  Float_t track_eta[NTRACK_MAX];
  Float_t track_phi[NTRACK_MAX];
  Float_t track_eta_emcal[NTRACK_MAX];
  Float_t track_phi_emcal[NTRACK_MAX];
  UChar_t track_quality[NTRACK_MAX];
  UChar_t track_its_ncluster[NTRACK_MAX];
  Float_t track_its_chi_square[NTRACK_MAX];
  Float_t track_dca_xy[NTRACK_MAX];
  Float_t track_dca_z[NTRACK_MAX];

  //Clusters
  UInt_t ncluster;
  Float_t cluster_e[NTRACK_MAX];
  Float_t cluster_e_max[NTRACK_MAX];
  Float_t cluster_e_cross[NTRACK_MAX];
  Float_t cluster_pt[NTRACK_MAX];
  Float_t cluster_eta[NTRACK_MAX];
  Float_t cluster_phi[NTRACK_MAX];
  Float_t cluster_iso_tpc_02[NTRACK_MAX];
  Float_t cluster_iso_tpc_04[NTRACK_MAX];
  Float_t cluster_iso_its_04[NTRACK_MAX];
  Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
  Float_t cluster_frixione_its_04_02[NTRACK_MAX];
  Float_t cluster_s_nphoton[NTRACK_MAX][4];
  UInt_t cluster_nmc_truth[NTRACK_MAX];
  unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
  Int_t cluster_ncell[NTRACK_MAX];
  UShort_t  cluster_cell_id_max[NTRACK_MAX];
  Float_t cluster_lambda_square[NTRACK_MAX][2];
  Float_t cell_e[17664];
  Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
  UChar_t cluster_nlocal_maxima[NTRACK_MAX];

  Float_t cluster_tof[NTRACK_MAX];
  Float_t cluster_iso_its_04_ue[NTRACK_MAX];
  Float_t cluster_iso_tpc_02_ue[NTRACK_MAX];
  Float_t cluster_iso_tpc_04_ue[NTRACK_MAX];

  // Jets
  UInt_t njet_ak04tpc;
  Float_t jet_ak04tpc_pt_raw[NTRACK_MAX];
  Float_t jet_ak04tpc_eta[NTRACK_MAX];
  Float_t jet_ak04tpc_phi[NTRACK_MAX];

  UInt_t njet_charged_truth_ak04;
  Float_t jet_charged_truth_ak04_pt[NTRACK_MAX];
  Float_t jet_charged_truth_ak04_eta[NTRACK_MAX];
  Float_t jet_charged_truth_ak04_phi[NTRACK_MAX];

  //MC
  unsigned int nmc_truth;
  Float_t mc_truth_pt[NTRACK_MAX];
  Float_t mc_truth_eta[NTRACK_MAX];
  Float_t mc_truth_phi[NTRACK_MAX];
  short mc_truth_pdg_code[NTRACK_MAX];
  short mc_truth_first_parent_pdg_code[NTRACK_MAX];
  char mc_truth_charge[NTRACK_MAX];
  Bool_t mc_truth_is_prompt_photon[NTRACK_MAX];

  Float_t mc_truth_first_parent_e[NTRACK_MAX];
  Float_t mc_truth_first_parent_pt[NTRACK_MAX];
  Float_t mc_truth_first_parent_eta[NTRACK_MAX];
  Float_t mc_truth_first_parent_phi[NTRACK_MAX];
  UChar_t mc_truth_status[NTRACK_MAX];
  //Float_t eg_cross_section;
  //Int_t   eg_ntrial;

  //Cluster Cut Summary
  fprintf(stderr, "%d: CLUSTER CUT SUMMARY \n ", __LINE__);
  fprintf(stderr, "%d: pT_max =  %f \n ", __LINE__, pT_max);
  fprintf(stderr, "%d: eta max = %f \n ", __LINE__, Eta_max);
  fprintf(stderr, "%d: SR Lambda max = %f \n ", __LINE__, srmax);
  fprintf(stderr, "%d: ncell min = %f \n ", __LINE__, Cluster_min);
  fprintf(stderr, "%d: Ecross/Emax = %f \n ", __LINE__, EcrossoverE_min);
  fprintf(stderr, "%d: Dist. bad channel = %f \n ", __LINE__, Cluster_DtoBad);

  YAML::Node filenames = configrunperiod["filelists"]["ntuples"]["gjmc"];
  for (YAML::const_iterator it = filenames.begin(); it != filenames.end(); it++) {
    std::string root_file = it->as<std::string>();
    std::cout << "Opening " << root_file << std::endl;
    TFile *file = TFile::Open((TString)root_file);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();

    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

    if (_tree_event == NULL) {
      _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
      if (_tree_event == NULL) {
        std::cout << " fail " << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    // Set the branch addresses of the branches in the TTrees
    //event Addresses
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

    //_tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    //_tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);

    // MC addresses
    _tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
    _tree_event->SetBranchAddress("mc_truth_is_prompt_photon", mc_truth_is_prompt_photon);

    // Jet addresses
    _tree_event->SetBranchAddress("njet_ak04tpc", &njet_ak04tpc);
    _tree_event->SetBranchAddress("jet_ak04tpc_pt_raw", jet_ak04tpc_pt_raw);
    _tree_event->SetBranchAddress("jet_ak04tpc_eta", jet_ak04tpc_eta);
    _tree_event->SetBranchAddress("jet_ak04tpc_phi", jet_ak04tpc_phi);

    _tree_event->SetBranchAddress("njet_charged_truth_ak04", &njet_charged_truth_ak04);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_pt", jet_charged_truth_ak04_pt);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_eta", jet_charged_truth_ak04_eta);
    _tree_event->SetBranchAddress("jet_charged_truth_ak04_phi", jet_charged_truth_ak04_phi);


    //IMPORTANT BOOLEAN VARIABLES
    Bool_t Signal = false;
    Bool_t Background = false;
    Bool_t Isolated = false;

    Long64_t nentries = _tree_event->GetEntries();

    //MAIN CORRELATION LOOP
    for (Long64_t ievent = 0; ievent < nentries ; ievent++) {
      //for(Long64_t ievent = 0; ievent < 10000 ; ievent++){
      _tree_event->GetEntry(ievent);
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

      Float_t purity_weight = 0;
      Float_t BR_purity_weight = 0;
      bool first_cluster = true;
      //if (not(first_cluster)) continue;

      //Event Selection
      if (TMath::Abs(primary_vertex[2]) > 10) continue;
      if (primary_vertex[2] == 0.00) continue;
      if (do_pile && is_pileup_from_spd_5_08) continue;

      // Jet matching - once per event
      vector<pair<int, int>> matchedJetIndices;
      set<int> unmatchedTruth;
      set<int> unmatchedReco;
      pair<int, int> matchedIndex;

      // first, make all reco jets unmatched
      for (int ireco = 0; ireco < njet_ak04tpc; ireco++) {
        unmatchedReco.insert(ireco);
      }

      for (int itruth = 0; itruth < njet_charged_truth_ak04; itruth++) {
        float minDistance2 = 99999999999;
        int minDistanceIndex = -1;

        float eta_truth = jet_charged_truth_ak04_eta[itruth];
        float phi_truth = jet_charged_truth_ak04_phi[itruth];

        for (int ireco = 0; ireco < njet_ak04tpc; ireco++) {
          float eta_reco = jet_ak04tpc_eta[ireco];
          float phi_reco = jet_ak04tpc_phi[ireco];
          float distance2 = ((eta_truth - eta_reco) * (eta_truth - eta_reco)) + ((phi_truth - phi_reco) * (phi_truth - phi_reco));

          if (distance2 < minDistance2) {
            minDistance2 = distance2;
            minDistanceIndex = ireco;
          }
        }

        if (minDistanceIndex > -1) {
          if (abs(jet_ak04tpc_eta[minDistanceIndex]) > jet_eta_max) {
            unmatchedTruth.insert(itruth);
          } else {
            matchedIndex = make_pair(itruth, minDistanceIndex);
            matchedJetIndices.push_back(matchedIndex);
            // remove the reco jet from the set of unmatched reco jets
            unmatchedReco.erase(minDistanceIndex);
          }
        } else {
          unmatchedTruth.insert(itruth);
        }
      }

      // Cluster loop
      for (ULong64_t n = 0; n < ncluster; n++) {
        if ( not(cluster_pt[n] > pT_min and cluster_pt[n] < pT_max)) continue; //select pt of photons
        if ( not(TMath::Abs(cluster_eta[n]) < Eta_max)) continue;           //cut edges of detector
        if ( not(cluster_ncell[n] >= Cluster_min)) continue;                 //removes clusters with 1 or 2 cells
        if ( not(cluster_e_cross[n] / cluster_e_max[n] > EcrossoverE_min)) continue; //removes "spiky" clusters
        if ( not(cluster_distance_to_bad_channel[n] >= Cluster_DtoBad)) continue; //removes clusters near bad channels
        if ( not(cluster_nlocal_maxima[n] < 3)) continue; //require to have at most 2 local maxima.

        // select on prompt photon
        bool isPrompt = false;
        for (UInt_t itruth = 0; itruth < cluster_nmc_truth[n]; itruth++) {
          unsigned short mcindex = cluster_mc_truth_index[n][itruth];
          if (mcindex == 65535) continue;
          if (!(abs(mc_truth_pdg_code[itruth]) == 11 || mc_truth_pdg_code[itruth] == 22)) continue;
          if (mc_truth_is_prompt_photon[itruth]) {
            isPrompt = true;
            break;
          }
        }
        if ( not(isPrompt)) continue;

        float isolation;
        if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
        else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
        else if (determiner == CLUSTER_ISO_ITS_04_SUB)
          isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n] - ue_estimate_its_const * 3.1416 * 0.4 * 0.4;
        else if (determiner == CLUSTER_ISO_TPC_04_SUB) {
          isolation = cluster_iso_tpc_04[n] + cluster_iso_tpc_04_ue[n] - ue_estimate_tpc_const * 3.1416 * 0.4 * 0.4;
        }
        else if (determiner == CLUSTER_ISO_TPC_02_SUB) {
          isolation = cluster_iso_tpc_02[n] + cluster_iso_tpc_02_ue[n] - ue_estimate_tpc_const * 3.1416 * 0.2 * 0.2;
        }
        else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
        else isolation = cluster_frixione_its_04_02[n];

        Isolated = GetIsIsolated(isolation, centrality_v0m, isoconfig);

        float shower = -1;
        if (shower_shape == "cluster_Lambda") {
          shower = cluster_lambda_square[n][0];
        }
        else if (shower_shape == "cluster_NN1") {
          shower = cluster_s_nphoton[n][1];
        }
        else if (shower_shape == "cluster_emax_over_e") {
          shower = cluster_e_max[n] / cluster_e[n];
        }

        Signal = (shower > srmin) and (shower < srmax);
        Background = (shower > brmin) and (shower < brmax);

        if (not(Signal)) continue;

        // get the centrality bin number
        int centbin = getCentBinNumber(centrality_v0m, centralityranges);

        // Filling the response matrices
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

        first_cluster = false;
      }//for nclusters

      matchedJetIndices.clear();
      unmatchedTruth.clear();
      unmatchedReco.clear();

    } //for nevents
    //}//end loop over samples
    file->Close();
    std::cout << std::endl;
  }

  // Write to fout
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