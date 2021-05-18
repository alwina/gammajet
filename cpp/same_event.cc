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

#include "yaml-cpp/yaml.h"

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_ISO_ITS_04_SUB, CLUSTER_ISO_TPC_02_SUB, CLUSTER_ISO_TPC_04_SUB, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};


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
  if (argc < 1) {
    fprintf(stderr, "Format: [command] [config file]\n");
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");
  bool Is_pp = false;

  //Config File
  YAML::Node config = YAML::LoadFile(argv[1]);
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
  double cluster_time = 20;

  double jet_pt_min = 0;
  double jet_pt_max = 500;
  double jet_eta_max = 10;

  bool do_pile = false;

  double iso_max = 0;
  double noniso_min = 0;
  double noniso_max = 0;
  //double deta_max = 0;
  isolationDet determiner = CLUSTER_ISO_ITS_04;
  std::string shower_shape = "DNN";
  std::string purity_deviation = "None";

  bool TPC_Iso_Flag = false;

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

  if (config["clustercuts"]["data"]["cluster_tof"]) {
    cluster_time = config["clustercuts"]["data"]["cluster_tof"]["max"].as<double>();
  }

  if (config["isolation"]) {
    iso_max = config["isolation"]["isocut"].as<double>();
    noniso_min = config["isolation"]["antiisocutlow"].as<double>();
    noniso_max = config["isolation"]["antiisocuthigh"].as<double>();
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
  Setting up THnSparses
  hCorrSR: cluster-jet correlations for the signal region
  hTrigSR: counting the number of clusters in each bin in the signal region
  hCorrBR: cluster-jet correlations for the bkg region
  hTrigBR: counting the number of clusters in each bin in the bkg region
  --------------------------------------------------------------*/

  // dimensions: centrality, cluster pT
  Int_t ndimTrig = 2;
  Int_t nbinsTrig[ndimTrig] = {10, 50};
  Double_t minbinsTrig[ndimTrig] = {0, 15};
  Double_t maxbinsTrig[ndimTrig] = {100, 40};
  THnSparseF* hTrigSR = new THnSparseF("hTrigSR", "Number of clusters (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
  THnSparseF* hTrigBR = new THnSparseF("hTrigBR", "Number of clusters (BR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

  Double_t trigSR[ndimTrig];
  Double_t trigBR[ndimTrig];

  // dimensions: centrality, cluster pT, delta phi, jet pT, pT ratio
  Int_t ndimCorr = 5;
  Int_t nbinsCorr[ndimCorr] = {10, 50, 120, 120, 120};
  Double_t minbinsCorr[ndimCorr] = {0, 15, 0, 0, 0};
  Double_t maxbinsCorr[ndimCorr] = {100, 40, M_PI, 50, 2};
  THnSparseF* hCorrSR = new THnSparseF("hCorrSR", "Correlations (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
  THnSparseF* hCorrBR = new THnSparseF("hCorrBR", "Correlations (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
  hCorrSR->Sumw2();
  hCorrBR->Sumw2();

  Double_t corrSR[ndimCorr];
  Double_t corrBR[ndimCorr];

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

  //MC
  unsigned int nmc_truth;
  Float_t mc_truth_pt[NTRACK_MAX];
  Float_t mc_truth_eta[NTRACK_MAX];
  Float_t mc_truth_phi[NTRACK_MAX];
  short mc_truth_pdg_code[NTRACK_MAX];
  short mc_truth_first_parent_pdg_code[NTRACK_MAX];
  char mc_truth_charge[NTRACK_MAX];

  Float_t mc_truth_first_parent_e[NTRACK_MAX];
  Float_t mc_truth_first_parent_pt[NTRACK_MAX];
  Float_t mc_truth_first_parent_eta[NTRACK_MAX];
  Float_t mc_truth_first_parent_phi[NTRACK_MAX];
  UChar_t mc_truth_status[NTRACK_MAX];
  //Float_t eg_cross_section;
  //Int_t   eg_ntrial;

  YAML::Node filenames = config["filelists"]["data"];
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
    _tree_event->SetBranchStatus("*mc*", 0);

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

    // Jet addresses
    _tree_event->SetBranchAddress("njet_ak04tpc", &njet_ak04tpc);
    _tree_event->SetBranchAddress("jet_ak04tpc_pt_raw", jet_ak04tpc_pt_raw);
    _tree_event->SetBranchAddress("jet_ak04tpc_eta", jet_ak04tpc_eta);
    _tree_event->SetBranchAddress("jet_ak04tpc_phi", jet_ak04tpc_phi);


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
    fprintf(stderr, "%d: SR Lambda max = %f \n ", __LINE__, srmax);
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


      for (ULong64_t n = 0; n < ncluster; n++) {
        if ( not(cluster_pt[n] > pT_min and cluster_pt[n] < pT_max)) continue; //select pt of photons
        if ( not(TMath::Abs(cluster_eta[n]) < Eta_max)) continue;           //cut edges of detector
        if ( not(cluster_ncell[n] >= Cluster_min)) continue;                 //removes clusters with 1 or 2 cells
        if ( not(cluster_e_cross[n] / cluster_e_max[n] > EcrossoverE_min)) continue; //removes "spiky" clusters
        if ( not(cluster_distance_to_bad_channel[n] >= Cluster_DtoBad)) continue; //removes clusters near bad channels
        if ( not(cluster_nlocal_maxima[n] < 3)) continue; //require to have at most 2 local maxima.
        if (not(abs(cluster_tof[n]) < cluster_time)) continue;

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

        Isolated = (isolation < iso_max);

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

        float bkg_weight = 1.0;
        float track_weight = 1.0; //Fake Rate, smearing, efficiency

        if (Background and Isolated) {
          BR_purity_weight = (1.0 / Get_Purity_ErrFunction(cluster_pt[n], purity_deviation, Is_pp, TPC_Iso_Flag) - 1); //(1-p)/p = 1/p - 1

          trigBR[0] = centrality_v0m;
          trigBR[1] = cluster_pt[n];
          hTrigBR->Fill(trigBR);
        }

        if (Signal and Isolated) {
          purity_weight = 1.0 / Get_Purity_ErrFunction(cluster_pt[n], purity_deviation, Is_pp, TPC_Iso_Flag);

          trigSR[0] = centrality_v0m;
          trigSR[1] = cluster_pt[n];
          hTrigSR->Fill(trigSR);
        }

        //Jet Loop

        for (ULong64_t ijet = 0; ijet < njet_ak04tpc; ijet++) {
          if (jet_ak04tpc_pt_raw[ijet] < jet_pt_min) continue;
          if (jet_ak04tpc_pt_raw[ijet] > jet_pt_max) continue;
          if (abs(jet_ak04tpc_eta[ijet]) > jet_eta_max) continue;

          // Observables: delta phi, jet pT, pT ratio
          Float_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - jet_ak04tpc_phi[ijet]));
          Float_t jetpt = jet_ak04tpc_pt_raw[ijet];
          Float_t ptratio = jetpt / cluster_pt[n];

          if (Signal and Isolated) {
            corrSR[0] = centrality_v0m;
            corrSR[1] = cluster_pt[n];
            corrSR[2] = deltaphi;
            corrSR[3] = jetpt;
            corrSR[4] = ptratio;
            hCorrSR->Fill(corrSR, purity_weight);
          }

          if (Background and Isolated) {
            corrBR[0] = centrality_v0m;
            corrBR[1] = cluster_pt[n];
            corrBR[2] = deltaphi;
            corrBR[3] = jetpt;
            corrBR[4] = ptratio;
            hCorrBR->Fill(corrBR, BR_purity_weight);
          }
        }//for ijets
        first_cluster = false;
      }//for nclusters
    } //for nevents
    file->Close();
  }

  // Write to fout
  TFile* fout;
  fout = new TFile("sameEvent.root", "RECREATE");
  std::cout << "Writing to file" << std::endl;

  hTrigSR->Write();
  hCorrSR->Write();
  hTrigBR->Write();
  hCorrBR->Write();

  fout->Close();
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}