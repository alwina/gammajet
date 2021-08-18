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
#include "shared_defs.h"

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_ISO_ITS_04_SUB, CLUSTER_ISO_TPC_02_SUB, CLUSTER_ISO_TPC_04_SUB, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};


int main(int argc, char *argv[])
{
  if (argc < 1) {
    fprintf(stderr, "Format: [command] [config file]\n");
    exit(EXIT_FAILURE);
  }

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

  //Cluster Cut Summary
  fprintf(stderr, "%d: CLUSTER CUT SUMMARY \n ", __LINE__);
  fprintf(stderr, "%d: pT_max =  %f \n ", __LINE__, pT_max);
  fprintf(stderr, "%d: eta max = %f \n ", __LINE__, Eta_max);
  fprintf(stderr, "%d: SR Lambda max = %f \n ", __LINE__, srmax);
  fprintf(stderr, "%d: ncell min = %f \n ", __LINE__, Cluster_min);
  fprintf(stderr, "%d: Ecross/Emax = %f \n ", __LINE__, EcrossoverE_min);
  fprintf(stderr, "%d: Dist. bad channel = %f \n ", __LINE__, Cluster_DtoBad);
  fprintf(stderr, "%d: cluster tof = %f \n ", __LINE__, cluster_time);

  /*---------------------------------------------------------------
  Using low level hdf5 API for clusters and jets
  ---------------------------------------------------------------*/
  //Procedure goes something like this:
  //1. Get dataset from the appropriate hdf5 file
  //2. Define dataspace as intermediary between file and array
  //3. Get Dimensinal information from dataspace, used for allocation
  //4. Read hyperslab from disk into hyperslab in memory (annoying)
  //5. Read data in memory dataspace into local arrays
  //6. Repeat steps 2&3 with different offsets.

  //Terms:
  //hyperslap: an n-dimensional "chunk" of the hdf5 data (disk or memory)
  //offset: the n-dimensional starting point of the hyperslab.
  //dataspace: defines size and shape of dataset or attribute
  
  //Documentation:
  //https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html
  //Instead of ranks as strict definitions (RANK_OUT for example), 
  //I use dimensions obtained by reading from the hdf5 file and dataset

  //Define triggered and min-bias H5File's
  H5File triggered_h5_file( triggered_hdf5_file_name, H5F_ACC_RDONLY ); 
  H5File MB_h5_file( MB_hdf5_file_name, H5F_ACC_RDONLY ); //file name from argv[2]

  //get cluster dataset from TRIGGERED file
  const H5std_string cluster_ds_name( "cluster" );
  DataSet cluster_dataset = triggered_h5_file.openDataSet( cluster_ds_name );

  //get event dataset from TRIGGERED file
  const H5std_string event_ds_name( "event" );
  DataSet event_dataset = triggered_h5_file.openDataSet( event_ds_name );

  //Mixed Event pairings are found in the mix dateset, in the TRIGGERED file
  const H5std_string mix_ds_name( "mix" );
  DataSet mix_dataset = triggered_h5_file.openDataSet( mix_ds_name );

  //get Jet Dateset from MIN-BIAS file
  const H5std_string jet_ds_name( "jet" );
  DataSet jet_dataset = MB_h5_file.openDataSet( jet_ds_name );
  /* DataSet mb_event_dataset = MB_h5_file.openDataSet( event_ds_name ); */
  //FIXME:Add MB event information

  //Get DataSpaces from datasets
  DataSpace cluster_dataspace = cluster_dataset.getSpace();
  DataSpace event_dataspace = event_dataset.getSpace();
  DataSpace mix_dataspace = mix_dataset.getSpace();
  DataSpace jet_dataspace = jet_dataset.getSpace();
  /* DataSpace mb_event_dataspace = event_dataset.getSpace(); */

  //Load the dimensions of datasets from file, to be used in dataspace/hyperslab
  //first get the number of dimensions with ExtentNdims
  //Then get the length of each dimension with ExtentDims.
  const int cluster_rank = cluster_dataspace.getSimpleExtentNdims();
  hsize_t clusterdims[cluster_rank];
  cluster_dataspace.getSimpleExtentDims(clusterdims, NULL);
  UInt_t ncluster_max = clusterdims[1];
  UInt_t Ncluster_Vars = clusterdims[2];
  fprintf(stderr, "\n%s:%d: number of cluster variables = %i\n", __FILE__, __LINE__, Ncluster_Vars);

  const int mix_rank = mix_dataspace.getSimpleExtentNdims();
  hsize_t mixdims[mix_rank];
  mix_dataspace.getSimpleExtentDims(mixdims, NULL);
  UInt_t nmix = mixdims[2]
  fprintf(stderr, "\n%s:%d: number of mixed events from file = %i\n", __FILE__, __LINE__, nmix);

  const int event_rank = event_dataspace.getSimpleExtentNdims();
  hsize_t eventdims[event_rank];
  event_dataspace.getSimpleExtentDims(eventdims, NULL);
  hsize_t nevent_max = eventdims[0];
  UInt_t nevent_Vars = eventdims[1];
  fprintf(stderr, "\n%s:%d: number of event variables = %i\n", __FILE__, __LINE__, nevent_Vars);//FIXME: ->Nevent_Vars

  const int jet_rank = jet_dataspace.getSimpleExtentNdims();
  hsize_t jetdims[jet_rank];
  jet_dataspace.getSimpleExtentDims(jetdims, NULL);
  UInt_t njet_max = jetdims[1];
  UInt_t Njet_Vars = jetdims[2];
  fprintf(stderr, "\n%s:%d: number of jet variables = %i\n", __FILE__, __LINE__, Njet_Vars);

  //Local arrays the hyperslabs will be fed into, and ultimatley used in LOOPS
  float cluster_data_out[1][ncluster_max][Ncluster_Vars];
  float mix_data_out[1][1][nmix];
  float event_data_out[1][nevent_Vars];
  float jet_data_out[1][njet_max][Njet_Vars];

  //Define hyperslab size and offset in FILE;
  hsize_t event_offset[2] = {0, 0};
  hsize_t event_count[2] = {1, nevent_Vars};
  hsize_t cluster_offset[3] = {0, 0, 0};
  hsize_t cluster_count[3] = {1, ncluster_max, Ncluster_Vars};
  hsize_t mix_offset[3] = {0, 0, 0};
  hsize_t mix_count[3] = {1, 1, nmix};
  hsize_t jet_offset[3] = {0, 0, 0};
  hsize_t jet_count[3] = {1, njet_max, Njet_Vars};

  /* The Offset is how we iterate over the entire hdf5 file. */
  /* For example, To obtain data for event 68, set the */
  /* offset's to {68, njet_max, Njet_Vars}. */
  event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
  cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
  mix_dataspace.selectHyperslab( H5S_SELECT_SET, mix_count, mix_offset );
  jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");

  //Define the memory dataspace in which to place hyperslab
  DataSpace event_memspace(event_rank, eventdims );
  DataSpace cluster_memspace(cluster_rank, clusterdims );
  DataSpace mix_memspace(mix_rank, mixdims );
  DataSpace jet_memspace(jet_rank, jetdims );

  //Define DataSpace offset for hypreslab starting at begining:
  //We want to read the entire "chunk" of the DataSpace to an array
  //So we set the DataSpace offset to [0]. Can be important if low memory
  hsize_t event_offset_out[2] = {0};
  hsize_t cluster_offset_out[3] = {0};
  hsize_t mix_offset_out[3] = {0};
  hsize_t jet_offset_out[3] = {0};

  //define size of memory to be read into memory hyperslab (aka memspace)
  hsize_t event_count_out[2] = {1, nevent_Vars};
  hsize_t cluster_count_out[3] = {1, ncluster_max, Ncluster_Vars};
  hsize_t mix_count_out[3] = {1, 1, nmix};
  hsize_t jet_count_out[3] = {1, njet_max, Njet_Vars};

  //define space in memory for hyperslab, then write from file to memory
  event_memspace.selectHyperslab( H5S_SELECT_SET, event_count_out, event_offset_out );
  cluster_memspace.selectHyperslab( H5S_SELECT_SET, cluster_count_out, cluster_offset_out );
  mix_memspace.selectHyperslab( H5S_SELECT_SET, mix_count_out, mix_offset_out );
  jet_memspace.selectHyperslab( H5S_SELECT_SET, jet_count_out, jet_offset_out );

  //FINALLY read memory dataspace into arrays we can use easily. 
  event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );
  cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
  mix_dataset.read( mix_data_out, PredType::NATIVE_FLOAT, mix_memspace, mix_dataspace );
  jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );

  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "datasets succesfully read into array");

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

    //MAIN CORRELATION LOOP
    for (Long64_t ievent = 0; ievent < nentries ; ievent++) {
      // for some reason, loading these 2 events from this ntuple causes a segfault
      if (root_file == "/global/project/projectdirs/alice/NTuples/PbPb/15o_pass2_cluster15.root") {
        if (ievent == 7894 || ievent == 7895) {
          std::cout << std::endl << "skipping event " << ievent;
          if (ievent == 7895) {
            std::cout << std::endl;
          }
          continue;
        }
      }
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

        Isolated = GetIsIsolated(isolation, centrality_v0m, config["isolation"]);

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
          BR_purity_weight = (1.0 / getPurity(cluster_pt[n], centrality_v0m, config["purity"]) - 1); //(1-p)/p = 1/p - 1

          trigBR[0] = centrality_v0m;
          trigBR[1] = cluster_pt[n];
          hTrigBR->Fill(trigBR);
        }

        if (Signal and Isolated) {
          purity_weight = 1.0 / getPurity(cluster_pt[n], centrality_v0m, config["purity"]);

          trigSR[0] = centrality_v0m;
          trigSR[1] = cluster_pt[n];
          hTrigSR->Fill(trigSR);
        }

        //Jet Loop
