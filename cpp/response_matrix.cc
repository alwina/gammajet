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

using namespace std;

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_ISO_ITS_04_SUB, CLUSTER_ISO_TPC_04_SUB, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};


Float_t Get_Purity_ErrFunction(Float_t pT_GeV, std::string deviation, bool Is_pp, bool Use_TPC) {

  // fprintf(stderr,"\n USE TPC FLAG = ");
  // fputs(Use_TPC ? "true \n" : "false \n", stdout);

  // fprintf(stderr,"\n Proton-Proton FLAG = ");
  // fputs(Is_pp ? "true \n" : "false \n", stdout);

  Float_t purity_val = 0;

  Float_t par[3] = {0.549684905516,
                    8.44338685256,
                    13.3454091464
                   };

  //purity_val = par[0]*TMath::Erf((pT_GeV-par[1])/par[2]);
  //fprintf(stderr,"\n %d: Purity ITS = %f",__LINE__,purity_val);

  if (Use_TPC) { //pPb
    par[0] = 0.569429959156;
    par[1] = 8.1906528936;
    par[2] = 11.8993694765;

    // purity_val = par[0]*TMath::Erf((pT_GeV-par[1])/par[2]);
    // fprintf(stderr,"\n %d: Purity TPC = %f",__LINE__,purity_val);
    // fprintf(stderr,"\n %d: Purity TPC = %f",__LINE__,purity_val);

  }

  if (Is_pp) {
    // fprintf(stderr,"\n");
    // fprintf(stderr,"\n PP SELECTED \n");
    par[0] = 0.500229283252;
    par[1] = 9.016920902665;
    par[2] = 11.373299838596;
  }
  //order of conditionals ensures pp always overwrights TPC purity

  if (strcmp(deviation.data(), "Plus") == 0) {
    par[0] = 0.60750016509;
    par[1] = 7.05184155403;
    par[2] = 13.6116163603;
  }

  if (strcmp(deviation.data(), "Minus") == 0) {
    par[0] = 0.479958593235;
    par[1] = 9.05392932723;
    par[2] = 10.2061359452;
  }

  purity_val = par[0] * TMath::Erf((pT_GeV - par[1]) / par[2]);
  //fprintf(stderr,"\n %d: Cluster pT = %f, Purity = %f \n",__LINE__,pT_GeV,purity_val);
  return purity_val;
}


int main(int argc, char *argv[])
{
  if (argc < 3) {
    fprintf(stderr, "Format: [command] [root file] [pp or pPb] \n");
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");


  bool Is_pp = false;


  std::string coll_system = argv[2];
  if (strcmp(coll_system.c_str(), "pp") == 0)
    Is_pp = true;

  if (Is_pp)
    fprintf(stderr, "\n PROTON PROTON SELECTED \n \n");

  //Config File
  FILE* config = fopen("Corr_config.yaml", "r");
  double DNN_min = 0;
  double DNN_max = 0;
  double DNN_Bkgd = 0;
  double Lambda0_cut = 0;
  double Emax_min = 0;
  double Emax_max = 0;
  double pT_min = 0;
  double pT_max = 0;
  double Eta_max = 0;
  double Cluster_min = 0;
  float Cluster_DtoBad = 0;
  UChar_t Cluster_NLocal_Max = 0;
  double EcrossoverE_min = 0;
  double cluster_time = 20;

  bool do_pile = false;

  float track_pT_min = 0.0;
  float track_pT_max = 0.0;
  int Track_Cut_Bit = 0;
  int track_chi_max = 0;
  double iso_max = 0;
  double noniso_min = 0;
  double noniso_max = 0;
  //double deta_max = 0;
  isolationDet determiner = CLUSTER_ISO_ITS_04;
  int n_eta_bins = 0;
  int n_phi_bins = 0;
  std::string shower_shape = "DNN";
  std::string purity_deviation = "None";

  bool TPC_Iso_Flag = false;

  // parse config file
  // check for existence first, then cast as appropriate
  if (config["DNN_min"]) {
    DNN_min = config["DNN_min"].as<double>();
  }

  if (config["DNN_max"]) {
    DNN_max = config["DNN_max"].as<double>();
  }

  if (config["DNN_BKGD"]) {
    DNN_Bkgd = config["DNN_BKGD"].as<double>();
  }

  if (config["Lambda0_cut"]) {
    Lambda0_cut = config["Lambda0_cut"].as<double>();
  }

  if (config["EMax_EClus_min"]) {
    Emax_min = config["EMax_EClus_min"].as<double>();
  }

  if (config["EMax_EClus_max"]) {
    Emax_max = config["EMax_EClus_max"].as<double>();
  }

  if (config["pT_min"]) {
    pT_min = config["pT_min"].as<double>();
  }

  if (config["pT_max"]) {
    pT_max = config["pT_max"].as<double>();
  }

  if (config["Eta_max"]) {
    Eta_max = config["Eta_max"].as<double>();
  }

  if (config["Cluster_min"]) {
    Cluster_min = config["Cluster_min"].as<double>();
  }

  if (config["Cluster_dist_to_bad_channel"]) {
    Cluster_DtoBad = config["Cluster_dist_to_bad_channel"].as<double>();
  }

  if (config["Cluster_N_Local_Maxima"]) {
    Cluster_NLocal_Max = config["Cluster_N_Local_Maxima"].as<double>();
  }

  if (config["EcrossoverE_min"]) {
    EcrossoverE_min = config["EcrossoverE_min"].as<double>();
  }

  if (config["Cluster_Time"]) {
    cluster_time = config["Cluster_Time"].as<double>();
  }

  if (config["iso_max"]) {
    iso_max = config["iso_max"].as<double>();
  }

  if (config["noniso_min"]) {
    noniso_min = config["noniso_min"].as<double>();
  }

  if (config["noniso_max"]) {
    noniso_max = config["noniso_max"].as<double>();
  }

  if (config["do_pileup_cut"]) {
    do_pile = config["do_pileup_cut"].as<bool>();
  }

  if (config["N_Phi_Bins"]) {
    n_phi_bins = config["N_Phi_Bins"].as<int>();
  }

  if (config["N_Eta_Bins"]) {
    n_eta_bins = config["N_Eta_Bins"].as<int>();
  }

  if (config["Cluster_isolation_determinant"]) {
    std::string determinant = config["Cluster_isolation_determinant"].as<std::string>();

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

  if (config["Shower_Shape"]) {
    shower_shape = config["Shower_Shape"].as<std::string>();
    std::cout << "Shower Shape: " << shower_shape << std::endl;
  }

  if (config["Purity_Dev"]) {
    purity_deviation = config["Purity_Dev"].as<std::string>();
    std::cout << "Purity Deviation Change: " << purity_deviation << std::endl;
  }

  //Purity Handling

  //Following purities for pT range: 12.5,13.2,14.4,15.8
  int  N_pT_Ranges = 5;
  float pT_Ranges[5] = {12.0, 15.0, 20.0, 25.0, 40.0};
  float purities[4] = {0};
  float purity_Uncertainties[4] = {0};
  float Cluster_Purity = 0;
  float Cluster_Purity_Uncertainty = 0;

  // if (strcmp(shower_shape.data(), "DNN") == 0){
  //   purities[0] = 0.207;
  //   purities[1] = 0.255;
  //   purities[2] = 0.326;
  //   purities[3] = 0.372;
  //   purities[4] = 0.447;
  //   purities[5] = 0.502;
  //   purities[6] = 0.533; //Extrapolating last bin
  //   purities[7] = 0.533;


  //   purity_Uncertainties[0] = 0.030;
  //   purity_Uncertainties[1] = 0.037;
  //   purity_Uncertainties[2] = 0.042;
  //   purity_Uncertainties[3] = 0.050;
  //   purity_Uncertainties[4] = 0.056;
  //   purity_Uncertainties[5] = 0.062;
  //   purity_Uncertainties[6] = 0.058;
  //   purity_Uncertainties[7] = 0.058;

  // }

  if (strcmp(shower_shape.data(), "Lambda") == 0) {
    purities[0] = 0.206;
    purities[1] = 0.341;
    purities[2] = 0.471;
    purities[3] = 0.546;

    purity_Uncertainties[0] = 0.0301;
    purity_Uncertainties[1] = 0.0305;
    purity_Uncertainties[2] = 0.0503;
    purity_Uncertainties[3] = 0.0572;

    if (Is_pp) {

      purities[0] = 0.198;
      purities[1] = 0.318;
      purities[2] = 0.470;
      purities[3] = 0.487;

      purity_Uncertainties[0] = 0.0490;
      purity_Uncertainties[1] = 0.0400;
      purity_Uncertainties[2] = 0.0712;
      purity_Uncertainties[3] = 0.1221;


    }
  }

  /*--------------------------------------------------------------
  Setting up RooUnfoldResponse objects
  --------------------------------------------------------------*/
  bool keepMisses = false;
  bool keepFakes = false;

  RooUnfoldResponse deltaphiResponse(10, 0, M_PI, 10, 0, M_PI);
  RooUnfoldResponse jetptResponse(9, 5, 50, 9, 5, 50);
  RooUnfoldResponse ptratioResponse(10, 0, 2, 10, 0, 2);


  //for (int iarg = 1; iarg < argc; iarg++) {
  int iarg = 1;
  TString root_file = (TString)argv[iarg];
  std::cout << "Opening: " << (TString)argv[iarg] << std::endl;

  TFile *file = TFile::Open(root_file);

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
  std::cout << " Total Number of entries in TTree: " << nentries << std::endl;

  //Cluster Cut Summary
  fprintf(stderr, "%d: CLUSTER CUT SUMMARY \n ", __LINE__);
  fprintf(stderr, "%d: pT_max =  %f \n ", __LINE__, pT_max);
  fprintf(stderr, "%d: eta max = %f \n ", __LINE__, Eta_max);
  fprintf(stderr, "%d: SR Lambda max = %f \n ", __LINE__, Lambda0_cut);
  fprintf(stderr, "%d: ncell min = %f \n ", __LINE__, Cluster_min);
  fprintf(stderr, "%d: Ecross/Emax = %f \n ", __LINE__, EcrossoverE_min);
  fprintf(stderr, "%d: Dist. bad channel = %f \n ", __LINE__, Cluster_DtoBad);
  fprintf(stderr, "%d: cluster tof = %f \n ", __LINE__, cluster_time);

  //MAIN CORRELATION LOOP

  fprintf(stderr, "\n Looping for main correlation functions \n");
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
        if (abs(jet_ak04tpc_eta[minDistanceIndex]) > 0.5) {
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
      if ( not(cluster_ncell[n] > Cluster_min)) continue;                 //removes clusters with 1 or 2 cells
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
      else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
      else isolation = cluster_frixione_its_04_02[n];

      Isolated = (isolation < iso_max);
      if (not(Isolated)) continue;

      if (strcmp(shower_shape.data(), "Lambda") == 0) {
        Signal = ((cluster_lambda_square[n][0] > 0.1) and (cluster_lambda_square[n][0] < Lambda0_cut));
        Background = (cluster_lambda_square[n][0] > 0.6);
      }

      else if (strcmp(shower_shape.data(), "DNN") == 0) {
        Signal = ( (cluster_s_nphoton[n][1] > DNN_min) && (cluster_s_nphoton[n][1] < DNN_max));
        Background = (cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < DNN_Bkgd);
      }

      else if (strcmp(shower_shape.data(), "EMax") == 0) {
        Signal = (cluster_e_max[n] / cluster_e[n] > Emax_max);
        Background = (cluster_e_max[n] / cluster_e[n] < Emax_min);
      }

      if (not(Signal)) continue;

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

        deltaphiResponse.Fill(deltaphi_reco, deltaphi_truth);
        jetptResponse.Fill(j_reco_pt, j_truth_pt);
        ptratioResponse.Fill(j_reco_pt / cluster_pt[n], j_truth_pt / cluster_pt[n]);
      }

      if (keepMisses) {
        for (int itruth : unmatchedTruth) {
          deltaphiResponse.Miss(TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - jet_charged_truth_ak04_phi[itruth])));
          jetptResponse.Miss(jet_charged_truth_ak04_pt[itruth]);
          ptratioResponse.Miss(jet_charged_truth_ak04_pt[itruth] / cluster_pt[n]);
        }
      }

      if (keepFakes) {
        for (int ireco : unmatchedReco) {
          deltaphiResponse.Fake(TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - jet_ak04tpc_phi[ireco])));
          jetptResponse.Fake(jet_ak04tpc_pt_raw[ireco]);
          ptratioResponse.Fake(jet_ak04tpc_pt_raw[ireco] / cluster_pt[n]);
        }
      }

      first_cluster = false;
    }//for nclusters

    matchedJetIndices.clear();
    unmatchedTruth.clear();
    unmatchedReco.clear();
    
  } //for nevents
  //}//end loop over samples

  // Write to fout
  TFile* fout;
  fout = new TFile("responseMatrix.root", "RECREATE");
  std::cout << "Writing to file" << std::endl;

  deltaphiResponse.Write();
  jetptResponse.Write();
  ptratioResponse.Write();

  fout->Close();
  file->Close();
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}