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

#include <omp.h>
#include "H5Cpp.h"

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_ISO_ITS_04_SUB, CLUSTER_ISO_TPC_02_SUB, CLUSTER_ISO_TPC_04_SUB, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

using namespace H5;

int main(int argc, char *argv[])
{
  if (argc < 6) {
    fprintf(stderr, "Format: [command] [config file] [root file] [MB hdf5 file] [PbPb or pp] [Mix Start] [Mix End]\n");
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

  YAML::Node purityconfig;
  YAML::Node isoconfig;

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
  }

  /*--------------------------------------------------------------
    Setting up THnSparses
hCorrSR: cluster-jet correlations for the signal region
hTrigSR: counting the number of clusters in each bin in the signal region
hCorrBR: cluster-jet correlations for the bkg region
hTrigBR: counting the number of clusters in each bin in the bkg region
--------------------------------------------------------------*/

  // dimensions: centrality, cluster pT
  const Int_t ndimTrig = 2;
  Int_t nbinsTrig[ndimTrig] = {10, 50};
  Double_t minbinsTrig[ndimTrig] = {0, 15};
  Double_t maxbinsTrig[ndimTrig] = {100, 40};
  THnSparseF* hTrigSR = new THnSparseF("hTrigSR", "Mixed Event Number of clusters (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
  THnSparseF* hTrigBR = new THnSparseF("hTrigBR", "Mixed Event Number of clusters (BR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

  THnSparseF* hnMixSR= new THnSparseF("hnMixSR", "Number of Mixed Events (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
  THnSparseF* hnMixBR= new THnSparseF("hnMixBR", "Number of Mixed Events (BR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

  Double_t trigSR[ndimTrig];
  Double_t trigBR[ndimTrig];
  Double_t nMixSR[ndimTrig];
  Double_t nMixBR[ndimTrig];


  // dimensions: centrality, cluster pT, delta phi, jet pT, pT ratio
  const Int_t ndimCorr = 5;
  Int_t nbinsCorr[ndimCorr] = {10, 50, 120, 120, 120};
  Double_t minbinsCorr[ndimCorr] = {0, 15, 0, 0, 0};
  Double_t maxbinsCorr[ndimCorr] = {100, 40, M_PI, 50, 2};
  THnSparseF* hCorrSR = new THnSparseF("hCorrSR", "Mixed Event Correlations (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
  THnSparseF* hCorrBR = new THnSparseF("hCorrBR", "Mixed Event Correlations (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
  hCorrSR->Sumw2();
  hCorrBR->Sumw2();

  Double_t corrSR[ndimCorr];
  Double_t corrBR[ndimCorr];

  //Debugging Histograms
  TH1F *hClusterpT = new TH1F("hClusterpT","Cluster pT Distribution",200,0,100);
  TH1F *hClusterLambda = new TH1F("hClusterLambda","Cluster Lambda Distribution",40,0,1.5);
  TH1F *hClustereta = new TH1F("hClustereta","Cluster eta Distribution",60,-3,3);
  TH1F *hClusterIso = new TH1F("hClusterIso","Cluster Iso Distribution",120,-10,50);

  /*---------------------------------------------------------------
    Mixed Event Book Keeping
    ---------------------------------------------------------------*/
  TH1I *SR_mixed_event_counter = new TH1I("SR_mixed_event_counter","Distribution of Mixed Event No.",1E6,0,1E6);
  TH1I *BR_mixed_event_counter = new TH1I("BR_mixed_event_counter","Distribution of Mixed Event No.",1E6,0,1E6);
  //GetEntries() of above is all that is needed. Distributions serve as additional check to mixing

  TH1D* z_vertices_triggered = new TH1D("Primary_Vertex_triggered", "Z-vertex (triggered)", 240, -12, 12);
  TH1D* z_vertices_MinBias = new TH1D("Primary_Vertex_MinBias", "Z-vertex (MinBias)", 240, -12, 12);
  TH1D* delta_z_vertices = new TH1D("Delta_Primary_Vertex", "#Delta V_z Distribution", 240, -12, 12);

  TH1D* multiplicity_triggered = new TH1D("multiplicity_triggered", "multiplicity (triggered)", 1000, 0, 1000);
  TH1D* multiplicity_MinBias = new TH1D("Multplicity_MinBias", "multiplicity (MinBias)", 500, 0, 1000);
  TH1D* delta_multiplicity = new TH1D("Delta_multiplicity", "#Delta Multiplicit Distribution", 500, 0, 1000);

  TH1D* flow_triggered = new TH1D("Flow_triggered", "Flow (triggered)", 500, -2, 2);
  TH1D* flow_MinBias = new TH1D("Flow_MinBias", "Flow (MinBias)", 500, -2, 2);
  TH1D* delta_flow= new TH1D("Delta_Flow", "#Delta Flow Distribution", 500, 0, 4);

  TH1D* centrality_triggered = new TH1D("centrality_triggered", "centrality (triggered)", 200, 0, 100);
  TH1D* centrality_MinBias = new TH1D("centrality_MinBias", "centrality (MinBias)", 200, 0, 100);
  TH1D* delta_centrality = new TH1D("Delta_centrality", "#Delta centrality Distribution", 200, 0, 100);


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
    Using low level hdf5 API for data 
    ---------------------------------------------------------------*/
  //Procedure goes something like this:
  //1. Get Datasets from the appropriate triggered or MB hdf5 file
  //2. Get Dataspaces from datasets, and record dimensions locally 
  //3. Define size and offset of hyperslab to be read from DISK
  //4. Define size and offset of hyperslab to read into MEMORY
  //5. Define local arrays to store information
  //6. DISK hyperslab -> MEMORY hyperslab -> local arrays
  //7. Iterate through hdf5 file by incrementing offset

  //Terms:
  //hyperslab: an n-dimensional "chunk" of the hdf5 data (disk or memory)
  //offset: the n-dimensional starting point of the hyperslab.
  //dataspace: defines size and shape of dataset or attribute

  //Documentation:
  //https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html
  //Instead of ranks as strict definitions (RANK_OUT for example), 
  //I use dimensions obtained by reading dataspace info the hdf5 file 

  //Triggered and MB hdf5 files
  const H5std_string triggered_hdf5_file_name(argv[2]);
  fprintf(stderr,(TString)argv[2]);
  H5File triggered_h5_file( triggered_hdf5_file_name, H5F_ACC_RDONLY ); 

  const H5std_string MB_hdf5_file_name(argv[3]);
  fprintf(stderr,(TString)argv[3]);
  H5File MB_h5_file( MB_hdf5_file_name, H5F_ACC_RDONLY );

  //get cluster dataset from TRIGGERED file
  const H5std_string cluster_ds_name( "cluster" );
  DataSet cluster_dataset = triggered_h5_file.openDataSet( cluster_ds_name );

  //get event dataset from TRIGGERED file
  const H5std_string event_ds_name( "event" );
  DataSet event_dataset = triggered_h5_file.openDataSet( event_ds_name );

  //Mixed Event pairings are found in the mix dateset, in the TRIGGERED file
  const H5std_string mix_ds_name( "mixing" );
  DataSet mix_dataset = triggered_h5_file.openDataSet( mix_ds_name );

  //get Jet Dateset from MIN-BIAS file
  const H5std_string jet_ds_name( "jet" );
  DataSet jet_dataset = MB_h5_file.openDataSet( jet_ds_name );
  //MINBIAS EVENT DATA
  DataSet mb_event_dataset = MB_h5_file.openDataSet( event_ds_name );

  //Get DataSpaces from datasets
  DataSpace cluster_dataspace = cluster_dataset.getSpace();
  DataSpace event_dataspace = event_dataset.getSpace();
  DataSpace mix_dataspace = mix_dataset.getSpace();
  DataSpace jet_dataspace = jet_dataset.getSpace();
  DataSpace mb_event_dataspace = mb_event_dataset.getSpace(); 

  //Load the dimensions of datasets from file, to be used in dataspace/hyperslab
  //first get the number of dimensions with ExtentNdims
  //Then get the length of each dimension with ExtentDims.
  const int cluster_rank = cluster_dataspace.getSimpleExtentNdims();
  hsize_t clusterdims[cluster_rank];
  cluster_dataspace.getSimpleExtentDims(clusterdims, NULL);
  UInt_t ncluster_max = clusterdims[1];
  UInt_t Ncluster_Vars = clusterdims[2];

  const int mix_rank = mix_dataspace.getSimpleExtentNdims();
  hsize_t mixdims[mix_rank];
  mix_dataspace.getSimpleExtentDims(mixdims, NULL);
  UInt_t nmix = mixdims[1];

  const int event_rank = event_dataspace.getSimpleExtentNdims();
  hsize_t eventdims[event_rank];
  event_dataspace.getSimpleExtentDims(eventdims, NULL);
  hsize_t nentries= eventdims[0]; //name break convention here to be similar same event corr. nentries = tree->GetEntries()
  UInt_t Nevent_Vars = eventdims[1];

  const int jet_rank = jet_dataspace.getSimpleExtentNdims();
  hsize_t jetdims[jet_rank];
  jet_dataspace.getSimpleExtentDims(jetdims, NULL);
  UInt_t njet_ak04tpc = jetdims[1];
  UInt_t Njet_Vars = jetdims[2];

  const int mb_event_rank = mb_event_dataspace.getSimpleExtentNdims();
  hsize_t mb_eventdims[mb_event_rank];
  mb_event_dataspace.getSimpleExtentDims(mb_eventdims, NULL);

  fprintf(stderr, "\n%s:%d: number of cluster variables = %i\n", __FILE__, __LINE__, Ncluster_Vars);
  fprintf(stderr, "\n%s:%d: number of mixed events from file = %i\n", __FILE__, __LINE__, nmix);
  fprintf(stderr, "\n%s:%d: number of event variables = %i\n", __FILE__, __LINE__, Nevent_Vars);
  fprintf(stderr, "\n%s:%d: number of jet variables = %i\n", __FILE__, __LINE__, Njet_Vars);

  //Block size should be determined by chunk size in to_hdf5. Usually 2000

  const int block_size = 4000;

  //Local arrays the hyperslabs will be fed into, and ultimatley used in LOOPS
  float cluster_data_out[block_size][ncluster_max][Ncluster_Vars];
  float mix_data_out[block_size][nmix];
  float event_data_out[block_size][Nevent_Vars];

  float jet_data_out[block_size][njet_ak04tpc][Njet_Vars];
  float mb_event_data_out[block_size][Nevent_Vars];

  //Define hyperslab size and offset to be read from FILE;
  hsize_t event_offset[2] = {0, 0};
  hsize_t event_count[2] = {block_size, Nevent_Vars};
  hsize_t cluster_offset[3] = {0, 0, 0};
  hsize_t cluster_count[3] = {block_size, ncluster_max, Ncluster_Vars};
  hsize_t mix_offset[2] = {0, 0};
  hsize_t mix_count[2] = {block_size, nmix};
  //MB
  //FIXME: Get MinBias to use block logic as well
  hsize_t jet_offset[3] = {0, 0, 0};
  hsize_t jet_count[3] = {block_size, njet_ak04tpc, Njet_Vars};
  hsize_t mb_event_offset[2] = {0, 0};
  hsize_t mb_event_count[2] = {block_size, Nevent_Vars};

  /* The Offset is how we iterate over the entire hdf5 file. */
  /* For example, To obtain data for event 68, set the */
  /* offset's to {68, njet_ak04tpc, Njet_Vars}. */
  event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
  cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
  mix_dataspace.selectHyperslab( H5S_SELECT_SET, mix_count, mix_offset );

  jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
  mb_event_dataspace.selectHyperslab( H5S_SELECT_SET, mb_event_count, mb_event_offset );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");

  //Define the memory dataspace in which to place hyperslab
  DataSpace event_memspace(event_rank, eventdims );
  DataSpace cluster_memspace(cluster_rank, clusterdims );
  DataSpace mix_memspace(mix_rank, mixdims );

  DataSpace jet_memspace(jet_rank, jetdims );
  DataSpace mb_event_memspace(mb_event_rank, mb_eventdims );

  //Define DataSpace offset for hypreslab starting at begining:
  //We want to read the entire "chunk" of the DataSpace to an array
  //So we set the DataSpace offset to [0]. Can be important if low memory
  hsize_t event_offset_out[2] = {0};
  hsize_t cluster_offset_out[3] = {0};
  hsize_t mix_offset_out[2] = {0};
  hsize_t jet_offset_out[3] = {0};
  hsize_t mb_event_offset_out[2] = {0};

  //define dimensions of hyperslab in memory (aka memspace)
  hsize_t event_count_out[2] = {block_size, Nevent_Vars};
  hsize_t cluster_count_out[3] = {block_size, ncluster_max, Ncluster_Vars};
  hsize_t mix_count_out[2] = {block_size, nmix};
  //FIXME: Get MinBias to Use block logic
  hsize_t jet_count_out[3] = {block_size, njet_ak04tpc, Njet_Vars};
  hsize_t mb_event_count_out[2] = {block_size, Nevent_Vars};

  //Apply the offset and dimensions from the previous two code blocks to the memspace
  event_memspace.selectHyperslab( H5S_SELECT_SET, event_count_out, event_offset_out );
  cluster_memspace.selectHyperslab( H5S_SELECT_SET, cluster_count_out, cluster_offset_out );
  mix_memspace.selectHyperslab( H5S_SELECT_SET, mix_count_out, mix_offset_out );
  jet_memspace.selectHyperslab( H5S_SELECT_SET, jet_count_out, jet_offset_out );
  mb_event_memspace.selectHyperslab( H5S_SELECT_SET, mb_event_count_out, mb_event_offset_out );

  //FINALLY use the well-defined memspace to read data from the dataspace, INTO the local array
  event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );
  cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
  mix_dataset.read( mix_data_out, PredType::NATIVE_FLOAT, mix_memspace, mix_dataspace );
  jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );
  mb_event_dataset.read( mb_event_data_out, PredType::NATIVE_FLOAT, mb_event_memspace, mb_event_dataspace );

  //First [block_size] number of events have just been read into local arrays

  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "datasets succesfully read into array");

  //FIXME: Add HDF5 Files to config files
  /* 
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
     */

  //IMPORTANT BOOLEAN VARIABLES
  Bool_t Signal = false;
  Bool_t Background = false;
  Bool_t Isolated = false;


  //MAIN CORRELATION LOOP
  /* nentries=100000; */
  /* #pragma omp parallel for */
  nmix = 40;
  for (Long64_t imix = 0; imix < nmix; imix++){
    std::cout<<std::endl<<"Mixed Event Number "<<imix<<" / "<<nmix<<std::endl;
    /* fprintf(stderr,"\n %s:%d: Mixed event = %lu, thread #%d", */
    /*     __FILE__,__LINE__,imix,omp_get_thread_num()); */

    //Grab 2000 Mixed Events at a time
    //Takes advantage of block structure used in pairing
    mb_event_offset[0]= imix * block_size;
    mb_event_dataspace.selectHyperslab( H5S_SELECT_SET, mb_event_count, mb_event_offset );
    mb_event_dataset.read( mb_event_data_out, PredType::NATIVE_FLOAT, mb_event_memspace, mb_event_dataspace );

    //FIXME: experiment with read outside of cluster loop. Would mean htrig histos need re-working as well
    jet_offset[0]= imix * block_size;
    jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
    jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );

    int offset = 0; //Offset for Triggered Events
    for (Long64_t ievent = 0; ievent < nentries; ievent++) {
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);
      int i = ievent % block_size;
      //reading from file is done every [block_size] number of events
      //this variable keeps track of the current increment within a block
      // as opposed to [ievent] which is looping through all events

      /* fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__,i,block_size-1); */
      if ((i == (block_size-1)) && (ievent!=nentries-1) && (ievent < nentries-block_size-1)) {//writes 1 block (2000 events) at a time. Faster/less memory

        offset += block_size;
        event_offset[0] = offset;
        event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
        event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );

        cluster_offset[0] = offset;
        cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
        cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );

        mix_offset[0] = offset;//need to get the mixed event pairing for this triggered event
        mix_dataspace.selectHyperslab( H5S_SELECT_SET, mix_count, mix_offset );
        mix_dataset.read( mix_data_out, PredType::NATIVE_FLOAT, mix_memspace, mix_dataspace );

      }

      //Set Event Variables
      float primary_vertex = event_data_out[i][0];//z-vertex
      float multiplicity = event_data_out[i][1];
      float v2 = event_data_out[i][2];
      float centrality_v0m = event_data_out[i][3];
      bool  is_pileup_from_spd_5_08 = event_data_out[i][4]; 
      float ue_estimate_its_const = event_data_out[i][5];
      float ue_estimate_tpc_const = event_data_out[i][6];

      Float_t purity_weight = 0;
      Float_t BR_purity_weight = 0;
      bool first_cluster = true;
      //if (not(first_cluster)) continue;

      //Event Selection
      if (TMath::Abs(primary_vertex) > 10) continue;
      if (primary_vertex == 0.00) continue;
      if (do_pile && is_pileup_from_spd_5_08) continue;

      z_vertices_triggered->Fill(primary_vertex);
      flow_triggered->Fill(v2);
      centrality_triggered->Fill(centrality_v0m);
      multiplicity_triggered->Fill(multiplicity);

      //Cluster Loop
      for (ULong64_t n = 0; n < ncluster_max; n++) {

        //See to_hdf5.cc for ROOT->HDF5 details
        if (std::isnan(cluster_data_out[i][n][0])) break; 
        float cluster_e = cluster_data_out[i][n][0];
        float cluster_pt = cluster_data_out[i][n][1];
        float cluster_eta = cluster_data_out[i][n][2];
        float cluster_phi = cluster_data_out[i][n][3];
        float cluster_lambda_square = cluster_data_out[i][n][4];
        float cluster_e_max = cluster_data_out[i][n][5];
        float cluster_e_cross = cluster_data_out[i][n][6];
        float cluster_iso_tpc_02 = cluster_data_out[i][n][7];
        float cluster_iso_tpc_03 = cluster_data_out[i][n][8];
        float cluster_iso_tpc_04 = cluster_data_out[i][n][9];
        float cluster_iso_its_02 = cluster_data_out[i][n][10];
        float cluster_iso_its_03 = cluster_data_out[i][n][11];
        float cluster_iso_its_04 = cluster_data_out[i][n][12];
        float cluster_frixione_tpc_04_02 = cluster_data_out[i][n][13];
        float cluster_frixione_its_04_02 = cluster_data_out[i][n][14];
        float cluster_s_nphoton0 = cluster_data_out[i][n][15];
        float cluster_mc_truth_index = cluster_data_out[i][n][16];//last dimension is length 32. Placeholder
        float cluster_ncell = cluster_data_out[i][n][17];
        float cluster_cell_id_max = cluster_data_out[i][n][18];
        float cluster_cell_e = cluster_data_out[i][n][19];//last dimension is length 17664. Placeholder
        float cluster_distance_to_bad_channel = cluster_data_out[i][n][20];
        float cluster_nlocal_maxima = cluster_data_out[i][n][21];
        float cluster_tof = cluster_data_out[i][n][22];
        float cluster_iso_its_02_ue = cluster_data_out[i][n][23];
        float cluster_iso_its_03_ue = cluster_data_out[i][n][24];
        float cluster_iso_its_04_ue = cluster_data_out[i][n][25];
        float cluster_iso_tpc_02_ue = cluster_data_out[i][n][26];
        float cluster_iso_tpc_03_ue = cluster_data_out[i][n][27];
        float cluster_iso_tpc_04_ue = cluster_data_out[i][n][28];
        float cluster_lambda1 = cluster_data_out[i][n][29];
        float cluster_s_nphoton = cluster_data_out[i][n][30];

        hClusterpT->Fill(cluster_pt);
        hClusterLambda->Fill(cluster_lambda_square);
        hClustereta->Fill(cluster_eta);
        hClusterIso->Fill(cluster_iso_tpc_04_ue);

        if ( not(cluster_pt > pT_min and cluster_pt < pT_max)) continue; //select pt of photons
        if ( not(TMath::Abs(cluster_eta) < Eta_max)) continue;           //cut edges of detector
        if ( not(cluster_ncell >= Cluster_min)) continue;                 //removes clusters with 1 or 2 cells
        if ( not(cluster_e_cross / cluster_e_max > EcrossoverE_min)) continue; //removes "spiky" clusters
        if ( not(cluster_distance_to_bad_channel >= Cluster_DtoBad)) continue; //removes clusters near bad channels
        if ( not(cluster_nlocal_maxima < 3)) continue; //require to have at most 2 local maxima.
        if (not(abs(cluster_tof) < cluster_time)) continue;

        float isolation;
        if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04;
        else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04;
        else if (determiner == CLUSTER_ISO_ITS_04_SUB)
          isolation = cluster_iso_its_04 + cluster_iso_its_04_ue - ue_estimate_its_const * 3.1416 * 0.4 * 0.4;
        else if (determiner == CLUSTER_ISO_TPC_04_SUB) {
          isolation = cluster_iso_tpc_04 + cluster_iso_tpc_04_ue - ue_estimate_tpc_const * 3.1416 * 0.4 * 0.4;
        }
        else if (determiner == CLUSTER_ISO_TPC_02_SUB) {
          isolation = cluster_iso_tpc_02 + cluster_iso_tpc_02_ue - ue_estimate_tpc_const * 3.1416 * 0.2 * 0.2;
        }
        else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02;
        else isolation = cluster_frixione_its_04_02;

        Isolated = GetIsIsolated(isolation, centrality_v0m, isoconfig);
        /* std::cout<<"ISOLATION "<<determiner<<" = "<<isolation<<std::endl; */

        float shower = -1;
        if (shower_shape == "cluster_Lambda") {
          shower = cluster_lambda_square;
        }
        else if (shower_shape == "cluster_NN1") {
          shower = cluster_s_nphoton;
        }
        else if (shower_shape == "cluster_emax_over_e") {
          shower = cluster_e_max / cluster_e;
        }

        Signal = (shower > srmin) and (shower < srmax);
        Background = (shower > brmin) and (shower < brmax);
        /* std::cout<<"Background Region = "<<brmin<<" to "<<brmax<<std::endl; */
        /* std::cout<<"Shower = "<<shower<<"; Shower Shape = "<<shower_shape<<"; Background BOOL = "<<Background<<std::endl; */
        /* std::cout<<"Centrality = "<<centrality_v0m<<"; Isolation Deteminer = "<<determiner<<"; Iso Val = "<<isolation<<std::endl; */

        float bkg_weight = 1.0;
        float track_weight = 1.0; //Fake Rate, smearing, efficiency

        if (Background and Isolated) {
          BR_purity_weight = (1.0 / getPurity(cluster_pt, centrality_v0m, purityconfig) - 1); //(1-p)/p = 1/p - 1

          trigBR[0] = centrality_v0m;
          trigBR[1] = cluster_pt;
          hTrigBR->Fill(trigBR);
        }

        if (Signal and Isolated) {
          purity_weight = 1.0 / getPurity(cluster_pt, centrality_v0m, purityconfig);

          trigSR[0] = centrality_v0m;
          trigSR[1] = cluster_pt;
          hTrigSR->Fill(trigSR);
        }

        //MIXING
        Long64_t mix_event =  mix_data_out[i][imix];
        if(mix_event  < 0) continue; //unpaired events set to negative numbers 
        int mix_index = mix_event % block_size;

        float mb_primary_vertex = mb_event_data_out[mix_index][0];//z-vertex
        float mb_multiplicity = mb_event_data_out[mix_index][1];
        float mb_v2 = mb_event_data_out[mix_index][2];
        float mb_centrality_v0m = mb_event_data_out[mix_index][3];

        z_vertices_MinBias->Fill(mb_primary_vertex);
        flow_MinBias->Fill(mb_v2);
        multiplicity_MinBias->Fill(mb_multiplicity);
        centrality_MinBias->Fill(mb_centrality_v0m);

        delta_z_vertices->Fill(TMath::Abs(primary_vertex-mb_primary_vertex));
        delta_flow->Fill(TMath::Abs(v2-mb_v2));
        delta_multiplicity->Fill(TMath::Abs(multiplicity-mb_multiplicity));
        delta_centrality->Fill(TMath::Abs(centrality_v0m-mb_centrality_v0m));

        //Add Delta Cuts on Event Pairing
        if (TMath::Abs(centrality_v0m-mb_centrality_v0m) > 10.) std::cout<<std::endl<<"SKIPPED MIXED EVENT"<<std::endl;
        if (TMath::Abs(centrality_v0m-mb_centrality_v0m) > 10.) continue;
        if (TMath::Abs(primary_vertex-mb_primary_vertex) > 2.) continue;
        if (TMath::Abs(v2-mb_v2) > 0.5) continue;

        if (Signal and Isolated) {
          nMixSR[0] = centrality_v0m;
          nMixSR[1] = cluster_pt;
          hnMixSR->Fill(nMixSR);
          SR_mixed_event_counter->Fill(mix_event);
        }

        if (Background and Isolated) {
          nMixBR[0] = centrality_v0m;
          nMixBR[1] = cluster_pt;
          hnMixBR->Fill(nMixBR);
          BR_mixed_event_counter->Fill(mix_event);
        }

        //Jet Loop
        for (ULong64_t ijet = 0; ijet < njet_ak04tpc; ijet++) {
          if (std::isnan(jet_data_out[mix_index][ijet][0])) break;
          float jet_ak04tpc_pt_raw = jet_data_out[mix_index][ijet][0];
          float jet_ak04tpc_eta= jet_data_out[mix_index][ijet][1];   
          float jet_ak04tpc_phi= jet_data_out[mix_index][ijet][2];   

          if (jet_ak04tpc_pt_raw < jet_pt_min) continue;
          if (jet_ak04tpc_pt_raw > jet_pt_max) continue;
          if (abs(jet_ak04tpc_eta) > jet_eta_max) continue;

          // Observables: delta phi, jet pT, pT ratio
          Float_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi - jet_ak04tpc_phi));
          Float_t jetpt = jet_ak04tpc_pt_raw;
          Float_t ptratio = jetpt / cluster_pt;

          if (Signal and Isolated) {
            corrSR[0] = centrality_v0m;
            corrSR[1] = cluster_pt;
            corrSR[2] = deltaphi;
            corrSR[3] = jetpt;
            corrSR[4] = ptratio;
            /* fprintf(stderr,"%i: cent = %f, pt = %f, deltaphi = %f, jetpt = %f, ptratio = %f\n",__LINE__,centrality_v0m,cluster_pt, deltaphi,jetpt,ptratio); */
            hCorrSR->Fill(corrSR, purity_weight);
          }

          if (Background and Isolated) {
            corrBR[0] = centrality_v0m;
            corrBR[1] = cluster_pt;
            corrBR[2] = deltaphi;
            corrBR[3] = jetpt;
            corrBR[4] = ptratio;
            hCorrBR->Fill(corrBR, BR_purity_weight);
          }

        }//jet loop
      }//cluster loop
    }//event loop
  }//mixed event loop

  TFile* fout;
  fout = new TFile("mixedEvent_h5.root", "RECREATE");
  // fout = new TFile((TString) configrunperiod["filelists"]["mixedevent"].as<std::string>(), "RECREATE");
  std::cout << "Writing to file" << std::endl;

  hTrigSR->Write();
  hCorrSR->Write();
  hnMixSR->Write();
  SR_mixed_event_counter->Write();

  hTrigBR->Write();
  hCorrBR->Write();
  hnMixBR->Write();
  BR_mixed_event_counter->Write();

  z_vertices_triggered->Write();
  z_vertices_MinBias->Write();
  delta_z_vertices->Write();

  multiplicity_triggered->Write();
  multiplicity_MinBias->Write();
  delta_multiplicity->Write();

  flow_triggered->Write();
  flow_MinBias->Write();
  delta_flow->Write();

  centrality_triggered->Write();
  centrality_MinBias->Write();
  delta_centrality->Write();

  hClusterpT->Write();
  hClusterLambda->Write();
  hClustereta->Write();
  hClusterIso->Write();

  fout->Close();
  std::cout << " ending " << std::endl;

  return EXIT_SUCCESS;
}
