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
#include "H5Cpp.h"

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_ISO_ITS_04_SUB, CLUSTER_ISO_TPC_04_SUB, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

using namespace H5;

int main(int argc, char *argv[])
{
  if (argc < 5) {
    fprintf(stderr, "Format: [command] [root file] [MB hdf5 file] [PbPb or pp] [Mix Start] [Mix End]\n");
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");

  //Config File
  FILE* config = fopen("../Corr_config.yaml", "r");
  if (!config){
    fprintf(stderr,"Config YAML not Found. Bailing\n");
    exit(EXIT_FAILURE);
  }
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

  bool TPC_Iso_Flag = false;

  // Zt bins
  //FIXME: Will have to likely set nztbins first, then initialize array
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins + 1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;

  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins + 1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;

  //FIXME: Obviously needs to be put in a header file.
  // Loop through config file
  char line[MAX_INPUT_LENGTH];
  while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
    if (line[0] == '#') continue;

    char key[MAX_INPUT_LENGTH];
    char dummy[MAX_INPUT_LENGTH];
    char value[MAX_INPUT_LENGTH];

    // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
    key[0] = '\0';
    value[0] = '\0';
    sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);

    //Read Config File: Detect Keys
    if (strcmp(key, "DNN_min") == 0) {
      DNN_min = atof(value);
      std::cout << "DNN_min: " << DNN_min << std::endl;
    }

    else if (strcmp(key, "DNN_max") == 0) {
      DNN_max = atof(value);
      std::cout << "DNN_max: " << DNN_max << std::endl;
    }

    else if (strcmp(key, "DNN_BKGD") == 0) {
      DNN_Bkgd = atof(value);
      std::cout << "DNN_BKGD: " << DNN_Bkgd << std::endl;
    }

    else if (strcmp(key, "Lambda0_cut") == 0) {
      Lambda0_cut = atof(value);
      std::cout << "Lambda0_cut: " << Lambda0_cut << std::endl;
    }

    if (strcmp(key, "EMax_EClus_min") == 0) {
      Emax_min = atof(value);
      std::cout << "EMax_EClus_min:" << Emax_min << std::endl;
    }

    else if (strcmp(key, "EMax_EClus_max") == 0) {
      Emax_max = atof(value);
      std::cout << "EMax_EClus_max: " << Emax_max << std::endl;
    }

    else if (strcmp(key, "pT_min") == 0) {
      pT_min = atof(value);
      std::cout << "pT_min: " << pT_min << std::endl;
    }

    else if (strcmp(key, "pT_max") == 0) {
      pT_max = atof(value);
      std::cout << "pT_max: " << pT_max << std::endl;
    }

    else if (strcmp(key, "Eta_max") == 0) {
      Eta_max = atof(value);
      std::cout << "Eta_max: " << Eta_max << std::endl;
    }
    else if (strcmp(key, "Cluster_min") == 0) {
      Cluster_min = atof(value);
      std::cout << "Cluster_min: " << Cluster_min << std::endl;
    }

    else if (strcmp(key, "Cluster_dist_to_bad_channel") == 0) {
      Cluster_DtoBad = atof(value);
      std::cout << "Cluster_DtoBad: " << Cluster_DtoBad << std::endl;
    }

    else if (strcmp(key, "Cluster_N_Local_Maxima") == 0) {
      Cluster_NLocal_Max = atof(value);
      std::cout << "Cluster_NLocal_Max: " << Cluster_NLocal_Max << std::endl;
    }

    else if (strcmp(key, "EcrossoverE_min") == 0) {
      EcrossoverE_min = atof(value);
      std::cout << "EcrossoverE_min; " << EcrossoverE_min << std::endl;
    }

    else if (strcmp(key, "Cluster_Time") == 0) {
      cluster_time = atof(value);
      std::cout << "Cluster_Time: " << cluster_time << std::endl;
    }

    else if (strcmp(key, "iso_max") == 0) {
      iso_max = atof(value);
      std::cout << "iso_max: " << iso_max << std::endl;
    }

    else if (strcmp(key, "noniso_min") == 0) {
      noniso_min = atof(value);
      std::cout << "noniso_min: " << noniso_min << std::endl;
    }

    else if (strcmp(key, "noniso_max") == 0) {
      noniso_max = atof(value);
      std::cout << "noniso_max: " << noniso_max << std::endl;
    }

    else if (strcmp(key, "do_pileup_cut") == 0) {
      if (strcmp(value, "true") == 0)
        do_pile = true;
      std::cout << "do_pileup_cut: " << do_pile << std::endl;
    }

    // else if (strcmp(key, "deta_max") == 0) {
    //     deta_max = atof(value);
    //     std::cout << "deta_max: " << deta_max << std::endl; }

    else if (strcmp(key, "N_Phi_Bins") == 0) {
      n_phi_bins = atoi(value);
      std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl;
    }

    else if (strcmp(key, "N_Eta_Bins") == 0) {
      n_eta_bins = atoi(value);
      std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl;
    }

    else if (strcmp(key, "Track_pT_Min") == 0) {
      track_pT_min = atof(value);
      std::cout << "Track Min pT: " << track_pT_min << std::endl;
    }

    else if (strcmp(key, "Track_pT_Max") == 0) {
      track_pT_max = atof(value);
      std::cout << "Track Max pT: " << track_pT_max << std::endl;
    }

    else if (strcmp(key, "Track_Cut_Bit") == 0) {
      Track_Cut_Bit = atoi(value);
      std::cout << "Track Cut Bit: " << Track_Cut_Bit << std::endl;
    }

    else if (strcmp(key, "Track_Chi_Max") == 0) {
      track_chi_max = atoi(value);
      std::cout << "Track Chi^2 Max: " << track_chi_max << std::endl;
    }


    else if (strcmp(key, "Zt_bins") == 0) {
      nztbins = -1;
      for (const char *v = value; *v != ']';) {
        while (*v != ']' && !isdigit(*v)) v++;
        nztbins++;
        while (*v != ']' && (isdigit(*v) || *v == '.')) v++;
      }

      ztbins = new float[nztbins + 1];
      int i = 0;
      for (const char *v = value; *v != ']' ;) {
        while (*v != ']' && !isdigit(*v)) v++;
        ztbins[i] = atof(v);
        i++;
        while (*v != ']' && (isdigit(*v) || *v == '.')) v++;
      }

      std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
      for (int i = 0; i <= nztbins; i++)
        std::cout << ztbins[i] << ", ";
      std::cout << "}\n";
    }

    else if (strcmp(key, "Pt_bins") == 0) {
      nptbins = -1;
      for (const char *v = value; *v != ']';) {
        while (*v != ']' && !isdigit(*v)) v++;
        nptbins++;
        while (*v != ']' && (isdigit(*v) || *v == '.')) v++;
      }

      ptbins = new float[nptbins + 1];
      int i = 0;
      for (const char *v = value; *v != ']' ;) {
        while (*v != ']' && !isdigit(*v))  v++;
        ptbins[i] = atof(v);
        i++;
        while (*v != ']' && (isdigit(*v) || *v == '.')) v++;
      }

      std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
      for (int i = 0; i <= nptbins; i++)
        std::cout << ptbins[i] << ", ";
      std::cout << "}\n";
    }

    else if (strcmp(key, "Cluster_isolation_determinant") == 0) {
      if (strcmp(value, "cluster_iso_tpc_04") == 0) {
        determiner = CLUSTER_ISO_TPC_04;
        std::cout << "Isolation Variable: cluster_iso_tpc_04" << std::endl;
      }

      else if (strcmp(value, "cluster_iso_its_04") == 0) {
        determiner = CLUSTER_ISO_ITS_04;
        std::cout << "Isolation Variable: cluster_iso_its_04" << std::endl;
      }

      else if (strcmp(value, "cluster_iso_its_04_sub") == 0) {
        determiner = CLUSTER_ISO_ITS_04_SUB;
        std::cout << "Isolation Variable: cluster_iso_its_04_sub" << std::endl;
      }

      else if (strcmp(value, "cluster_iso_tpc_04_sub") == 0) {
        determiner = CLUSTER_ISO_TPC_04_SUB;
        TPC_Iso_Flag = true;
        std::cout << "Isolation Variable: cluster_iso_tpc_04_sub" << std::endl;
      }

      else if (strcmp(value, "cluster_frixione_tpc_04_02") == 0) {
        determiner = CLUSTER_FRIXIONE_TPC_04_02;
        std::cout << "Isolation Variable: cluster_frixione_tpc_04_02" << std::endl;
      }

      else if (strcmp(value, "cluster_frixione_its_04_02") == 0) {
        determiner = CLUSTER_FRIXIONE_ITS_04_02;
        std::cout << "Isolation Variable: cluster_frixione_its_04_02" << std::endl;
      }

      else {
        std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    else if (strcmp(key, "Shower_Shape") == 0) {
      shower_shape = value;
      std::cout << "Shower Shape: " << shower_shape.data() << std::endl;
      //if (strcmp(shower_shape.data(),"Lambda")== 0) std::cout<<"test worked"<<std::endl;
    }

    else std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;

  }
  //end Config Loop

  fclose(config);

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
  Long64_t mix_events[300];

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

  // Set the branch addresses of the branches in the TTrees
  _tree_event->SetBranchStatus("*mc*", 0);

  //event Addresses
  _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
  _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
  _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
  _tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
  _tree_event->SetBranchAddress("centrality_v0m", &centrality_v0m);
  _tree_event->SetBranchAddress("mixed_events", mix_events);

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
  fprintf(stderr, "%d: SR Lambda max = %f \n ", __LINE__, Lambda0_cut);
  fprintf(stderr, "%d: ncell min = %f \n ", __LINE__, Cluster_min);
  fprintf(stderr, "%d: Ecross/Emax = %f \n ", __LINE__, EcrossoverE_min);
  fprintf(stderr, "%d: Dist. bad channel = %f \n ", __LINE__, Cluster_DtoBad);
  fprintf(stderr, "%d: cluster tof = %f \n ", __LINE__, cluster_time);

  /*------------------------------------------------------------------------
    MIXED EVENTS. 
  ------------------------------------------------------------------------*/
  TH1I *SR_mixed_event_counter = new TH1I("SR_mixed_event_counter","Distribution of Mixed Event No.",1E6,0,1E6);
  TH1I *BR_mixed_event_counter = new TH1I("BR_mixed_event_counter","Distribution of Mixed Event No.",1E6,0,1E6);
  //Deicated Histograms, or must be added to THnSparse
  //GetEntries() of above is all that is needed. Distributions serve as additional check to mixing

  const H5std_string hdf5_file_name(argv[2]);
  TString hdf5_file = (TString)argv[2];
  fprintf(stderr,hdf5_file);

  //Track and Purity weights should not be used.
  //But mixing tails should cut on v2 (flow) in PbPb
  bool Is_PbPb = false;
  std::string coll_system = argv[3];
  
  if (strcmp(coll_system.c_str(), "PbPb") == 0)
  {
    Is_PbPb = true;
    fprintf(stderr, "\n LEAD LEAD SELECTED \n \n");
  }
  
  size_t mix_start = atoi(argv[4]);
  size_t mix_end = atoi(argv[5]);

  /* ---------------------------------------------------------------
  Using low level hdf5 API for jets
  ---------------------------------------------------------------*/
  /* See to_hdf5.cc for jet_vars. Useful map below: */
  /* Jet Variariable 0 = jet_ak04tpc_pt_raw */
  /* Jet Variariable 1 = jet_ak04tpc_eta_raw */
  /* Jet Variariable 2 = jet_ak04tpc_phi */

  //open hdf5: Define size of data from file, explicitly allocate memory in hdf5 space and array size
  const H5std_string jet_ds_name( "jet" );
  H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY ); //file name from argv[2]
  DataSet jet_dataset = h5_file.openDataSet( jet_ds_name );
  DataSpace jet_dataspace = jet_dataset.getSpace();

  //Load the dimensions of dataset from file, to be used in array/hyperslab
  const int jet_ndims = jet_dataspace.getSimpleExtentNdims();
  hsize_t jet_maxdims[jet_ndims];
  hsize_t jetdims[jet_ndims];
  jet_dataspace.getSimpleExtentDims(jetdims, jet_maxdims);
  UInt_t njet_max = jetdims[1];
  UInt_t Njet_Vars = jetdims[2];
  fprintf(stderr, "\n%s:%d: n jet variables = %i\n", __FILE__, __LINE__, Njet_Vars);

  //Define array hyperslab will be fed into
  float jet_data_out[1][njet_max][Njet_Vars];

  //Define hyperslab size and offset in  FILE;
  hsize_t jet_offset[3] = {0, 0, 0};
  hsize_t jet_count[3] = {1, njet_max, Njet_Vars};

  /* 
     The Offset is how we iterate over the entire hdf5 file.
     For example, To obtain data for event 68, set the
     offset's to {68, njet_max, Njet_Vars}.
     */


  jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");

  //Define the memory dataspace in which to place hyperslab
  const int RANK_OUT = 3; //# of Dimensions
  DataSpace jet_memspace( RANK_OUT, jetdims );

  //Define memory offset for hypreslab starting at begining:
  hsize_t jet_offset_out[3] = {0};

  //define Dimensions of array, for writing slab to array
  hsize_t jet_count_out[3] = {1, njet_max, Njet_Vars};

  //define space in memory for hyperslab, then write from file to memory
  jet_memspace.selectHyperslab( H5S_SELECT_SET, jet_count_out, jet_offset_out );
  jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "jet dataset read into array: OK");

  //MAIN CORRELATION LOOP

  fprintf(stderr, "\n Looping for main correlation functions \n");
  for (Long64_t ievent = 0; ievent < nentries ; ievent++) {
    /* for(Long64_t ievent = 0; ievent < 10000 ; ievent++){ */
    _tree_event->GetEntry(ievent);
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

    bool first_cluster = true;
    //if (not(first_cluster)) continue;

    //Event Selection
    if (TMath::Abs(primary_vertex[2]) > 10) continue;
    if (primary_vertex[2] == 0.00) continue;
    if (do_pile && is_pileup_from_spd_5_08) continue;


    for (ULong64_t n = 0; n < ncluster; n++) {
      /* fprintf(stderr,"%d: pT = %f, eta = %f, ncell = %i, e_cross_e = %f, d_to_bad = %f, tof = %f\n",__LINE__, */
      /*  cluster_pt[n],cluster_eta[n], cluster_ncell[n], cluster_e_cross[n], cluster_distance_to_bad_channel[n]); */ 
      if ( not(cluster_pt[n] > pT_min and cluster_pt[n] < pT_max)) continue; //select pt of photons
      if ( not(TMath::Abs(cluster_eta[n]) < Eta_max)) continue;           //cut edges of detector
      if ( not(cluster_ncell[n] > Cluster_min)) continue;                 //removes clusters with 1 or 2 cells
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
      else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
      else isolation = cluster_frixione_its_04_02[n];

      Isolated = (isolation < iso_max);
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

      //Not used in Mixing
      /* float bkg_weight = 1.0; */
      /* float track_weight = 1.0; //Fake Rate, smearing, efficiency */

      if (Background and Isolated) {
        trigBR[0] = centrality_v0m;
        trigBR[1] = cluster_pt[n];
        hTrigBR->Fill(trigBR);
      }

      if (Signal and Isolated) {
        trigSR[0] = centrality_v0m;
        trigSR[1] = cluster_pt[n];
        hTrigSR->Fill(trigSR);
      }

      //Jet Loop
      for (Long64_t imix = mix_start; imix < mix_end; imix++){
        Long64_t mix_event = mix_events[imix];
        fprintf(stderr,"\n %s:%d: Mixed event = %lu",__FILE__,__LINE__,mix_event);

        if(mix_event  < 0) continue; //unpaired triggered events set to negative pairings 

        //FIXME: Add Delta Cuts on Event Pairing Here

        //adjust offset for next mixed event
        jet_offset[0]=mix_event;
        jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
        jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );

        for (ULong64_t ijet = 0; ijet < njet_max; ijet++) {
          if (std::isnan(jet_data_out[0][ijet][0])) continue;//hdf5 filled with NaN up to njet_max
          float jetpt = jet_data_out[0][ijet][0];
          float jeteta = jet_data_out[0][ijet][1];   
          float jetphi = jet_data_out[0][ijet][2];   
          if ( jetpt < 5) continue;
          if ( jetpt > 50) continue;
          if ( jeteta> 0.5) continue;

          //FIXME: Move cut values to config file

          //From same_event.cc
          /* for (ULong64_t ijet = 0; ijet < njet_ak04tpc; ijet++) { */
          /*   if (jet_ak04tpc_pt_raw[ijet] < 5) continue; */
          /*   if (jet_ak04tpc_pt_raw[ijet] > 50) continue; */
          /*   if (abs(jet_ak04tpc_eta[ijet]) > 0.5) continue; */

          // Observables: delta phi, jet pT, pT ratio
          Float_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - jetphi));
          Float_t ptratio = jetpt / cluster_pt[n];

          if (Signal and Isolated) {
            corrSR[0] = centrality_v0m;
            corrSR[1] = cluster_pt[n];
            corrSR[2] = deltaphi;
            corrSR[3] = jetpt;
            corrSR[4] = ptratio;
            hCorrSR->Fill(corrSR);
            SR_mixed_event_counter->Fill(mix_event);

          }

          if (Background and Isolated) {
            corrBR[0] = centrality_v0m;
            corrBR[1] = cluster_pt[n];
            corrBR[2] = deltaphi;
            corrBR[3] = jetpt;
            corrBR[4] = ptratio;
            hCorrBR->Fill(corrBR);
            BR_mixed_event_counter->Fill(mix_event);
          }
        }//for ijets
        first_cluster = false;
        }//for mixed
      }//for nclusters
    } //for nevents
    //}//end loop over samples

    // Write to fout
    TFile* fout;
    fout = new TFile("mixedEvent.root", "RECREATE");
    std::cout << "Writing to file" << std::endl;

    hTrigSR->Write();
    hCorrSR->Write();
    SR_mixed_event_counter->Write();
    hTrigBR->Write();
    hCorrBR->Write();
    BR_mixed_event_counter->Write();

    //Seperate zt loops for easier file reading
    fout->Close();
    file->Close();
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
  }
