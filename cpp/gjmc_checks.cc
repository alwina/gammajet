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
  double cluster_time = 20;

  double jet_pt_min = 0;
  double jet_pt_max = 500;
  double jet_eta_max = 10;

  bool do_pile = false;

  //double deta_max = 0;
  std::string jettype = "ak04tpc";
  std::string isovar = "cluster_iso_its_04";
  std::string shower_shape = "DNN";
  std::string purity_deviation = "None";

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

    if (config["jettype"]) {
      jettype = config["jettype"].as<std::string>();
    }

    if (config["jetcuts"]) {
      jet_pt_min = config["jetcuts"]["jet_pt_raw"]["min"].as<double>();
      jet_pt_max = config["jetcuts"]["jet_pt_raw"]["max"].as<double>();
      jet_eta_max = config["jetcuts"]["jet_eta"]["max"].as<double>();
    }

    if (config["do_pileup_cut"]) {
      do_pile = config["do_pileup_cut"].as<bool>();
    }

    if (config["isolation"]) {
      isoconfig = config["isolation"];
      isovar = config["isolation"]["isovar"].as<std::string>();
      std::cout << "Isolation variable: " << isovar << std::endl;
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
  hTrigPhoton: number of real photons
  hTrigSR: number of real photons in SR
  hCorrSRTruth: real photon-truth jet correlations for the signal region
  hCorrSRAll: real photon-all jet correlations for the signal region
  --------------------------------------------------------------*/

  // dimensions: centrality, photon pT, isolation
  Int_t ndimTrig = 3;
  Int_t nbinsTrig[ndimTrig] = {10, 50, 200};
  Double_t minbinsTrig[ndimTrig] = {0, 15, -20};
  Double_t maxbinsTrig[ndimTrig] = {100, 40, 80};
  THnSparseF* hTrigPhoton = new THnSparseF("hTrigPhoton", "Number of real photons", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
  THnSparseF* hTrigSR = new THnSparseF("hTrigSR", "Number of real photons in SR", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

  Double_t trig[ndimTrig];

  // dimensions: centrality, cluster pT, isolation, delta phi, jet pT, pT ratio
  Int_t ndimCorr = 6;
  Int_t nbinsCorr[ndimCorr] = {10, 50, 200, 120, 100, 120};
  Double_t minbinsCorr[ndimCorr] = {0, 15, -20, 0, 0, 0};
  Double_t maxbinsCorr[ndimCorr] = {100, 40, 80, M_PI, 50, 2};
  THnSparseF* hCorrSRTruth = new THnSparseF("hCorrSRTruth", "Correlations (Truth jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
  THnSparseF* hCorrSRAll = new THnSparseF("hCorrSRAll", "Correlations (All jets)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
  hCorrSRTruth->Sumw2();
  hCorrSRAll->Sumw2();

  Double_t corr[ndimCorr];

  /*--------------------------------------------------------------
  Setting up local variables to be linked with ROOT branches
  --------------------------------------------------------------*/
  //Events
  Bool_t is_pileup_from_spd_5_08;
  Double_t primary_vertex[3];
  Float_t ue_estimate_its_const;
  Float_t ue_estimate_tpc_const;
  Float_t centrality_v0m;
  Bool_t isINT7;

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
  unsigned int cluster_nmc_truth[NTRACK_MAX];
  unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
  Int_t cluster_ncell[NTRACK_MAX];
  UShort_t  cluster_cell_id_max[NTRACK_MAX];
  Float_t cluster_lambda_square[NTRACK_MAX][2];
  Float_t cell_e[17664];
  Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
  UChar_t cluster_nlocal_maxima[NTRACK_MAX];
  Float_t cluster_5x5all[NTRACK_MAX];
  bool cluster_is_prompt[NTRACK_MAX];

  Float_t cluster_tof[NTRACK_MAX];
  Float_t cluster_iso_its_04_ue[NTRACK_MAX];
  Float_t cluster_iso_tpc_02_ue[NTRACK_MAX];
  Float_t cluster_iso_tpc_04_ue[NTRACK_MAX];

  // Jets
  UInt_t njet;
  Float_t jet_pt_raw[NTRACK_MAX];
  Float_t jet_eta[NTRACK_MAX];
  Float_t jet_phi[NTRACK_MAX];

  // Truth jets
  UInt_t njet_charged_truth;
  Float_t jet_charged_truth_pt[NTRACK_MAX];
  Float_t jet_charged_truth_eta[NTRACK_MAX];
  Float_t jet_charged_truth_phi[NTRACK_MAX];

  //MC
  unsigned int nmc_truth;
  Float_t mc_truth_pt[NTRACK_MAX];
  Float_t mc_truth_eta[NTRACK_MAX];
  Float_t mc_truth_phi[NTRACK_MAX];
  short mc_truth_pdg_code[NTRACK_MAX];
  short mc_truth_first_parent_pdg_code[NTRACK_MAX];
  char mc_truth_charge[NTRACK_MAX];
  bool mc_truth_is_prompt_photon[NTRACK_MAX];

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
  fprintf(stderr, "%d: SR max = %f \n ", __LINE__, srmax);
  fprintf(stderr, "%d: ncell min = %f \n ", __LINE__, Cluster_min);
  fprintf(stderr, "%d: Ecross/Emax = %f \n ", __LINE__, EcrossoverE_min);
  fprintf(stderr, "%d: Dist. bad channel = %f \n ", __LINE__, Cluster_DtoBad);
  fprintf(stderr, "%d: cluster tof = %f \n ", __LINE__, cluster_time);
  std::cout << "Jet type: " << jettype << std::endl;

  YAML::Node filenames = configrunperiod["filelists"]["ntuples"]["gjmc"];
  for (YAML::const_iterator it = filenames.begin(); it != filenames.end(); it++) {
    std::string root_file = it->as<std::string>();
    std::cout << "Opening " << root_file << std::endl;
    TFile *file = TFile::Open((TString)root_file);

    if (file == NULL) {
      std::cout << "Failed to open file" << std::endl;
      exit(EXIT_FAILURE);
    }

    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

    if (_tree_event == NULL) {
      _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
      if (_tree_event == NULL) {
        std::cout << " fail " << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    std::string aux_filename = root_file.replace(root_file.find(".root"), 5, "_AUX.root");
    std::cout << "Opening " << aux_filename << std::endl;
    TFile *auxfile = TFile::Open((TString)aux_filename);

    // for these tests, the aux file is non-optional
    if (auxfile == NULL) {
      std::cout << "Failed to open aux file" << std::endl;
      exit(EXIT_FAILURE);
    }

    TTree *auxtree;
    // the aux file is only needed in certain situations, so check for those situations
    // things should still work even without the aux file
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


    // Set the branch addresses of the branches in the TTrees
    //event Addresses
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
    _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
    _tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
    _tree_event->SetBranchAddress("centrality_v0m", &centrality_v0m);
    if (auxfile != NULL) {
      auxtree->SetBranchAddress("isINT7", &isINT7);
    }

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
      auxtree->SetBranchAddress("cluster_is_prompt", cluster_is_prompt);
    }

    //_tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    //_tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);

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

    //IMPORTANT BOOLEAN VARIABLES
    Bool_t Signal = false;
    Bool_t Isolated = false;

    Long64_t nentries = _tree_event->GetEntries();
    // nentries = 1000;

    //MAIN CORRELATION LOOP
    for (Long64_t ievent = 0; ievent < nentries ; ievent++) {
      _tree_event->GetEntry(ievent);
      if (auxfile != NULL) {
        auxtree->GetEntry(ievent);
      }
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

      //Event Selection
      if (!isINT7) continue;
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
        if ( not(cluster_is_prompt[n])) continue;
        
        float isolation;
        if (isovar == "cluster_iso_tpc_04") isolation = cluster_iso_tpc_04[n];
        else if (isovar == "cluster_iso_its_04") isolation = cluster_iso_its_04[n];
        else if (isovar == "cluster_iso_its_04_sub") {
          isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n] - ue_estimate_its_const * 3.1416 * 0.4 * 0.4;
        }
        else if (isovar == "cluster_iso_tpc_04_sub") {
          isolation = cluster_iso_tpc_04[n] + cluster_iso_tpc_04_ue[n] - ue_estimate_tpc_const * 3.1416 * 0.4 * 0.4;
        }
        else if (isovar == "cluster_iso_tpc_02_sub") {
          isolation = cluster_iso_tpc_02[n] + cluster_iso_tpc_02_ue[n] - ue_estimate_tpc_const * 3.1416 * 0.2 * 0.2;
        }
        else if (isovar == "cluster_frixione_tpc_04_02") isolation = cluster_frixione_tpc_04_02[n];
        else if (isovar == "cluster_frixione_its_04_02") isolation = cluster_frixione_its_04_02[n];
        else {
          std::cout << "ERROR: Isolation variable " << isovar << " not recognized. Aborting" << std::endl;
          exit(EXIT_FAILURE);
        }

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
        else if (shower_shape == "cluster_5x5all") {
          shower = cluster_5x5all[n];
        }

        Signal = (shower > srmin) and (shower < srmax);

        // count all real photons
        trig[0] = centrality_v0m;
        trig[1] = cluster_pt[n];
        trig[2] = isolation;
        hTrigPhoton->Fill(trig);

        // count SR photons
        if (Signal) {
          hTrigSR->Fill(trig);
        } else {
          // we don't have to do anything else with the non-SR photons
          continue;
        }
        
        // correlations with all jets
        for (ULong64_t ijet = 0; ijet < njet; ijet++) {
          if (jet_pt_raw[ijet] < jet_pt_min) continue;
          if (jet_pt_raw[ijet] > jet_pt_max) continue;
          if (abs(jet_eta[ijet]) > jet_eta_max) continue;

          // Observables: delta phi, jet pT, pT ratio
          Float_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - jet_phi[ijet]));
          Float_t jetpt = jet_pt_raw[ijet];
          Float_t ptratio = jetpt / cluster_pt[n];

          corr[0] = centrality_v0m;
          corr[1] = cluster_pt[n];
          corr[2] = isolation;
          corr[3] = deltaphi;
          corr[4] = jetpt;
          corr[5] = ptratio;
          hCorrSRAll->Fill(corr);
        }//for ijets

        // correlations with truth jets
        for (ULong64_t ijet = 0; ijet < njet_charged_truth; ijet++) {
          if (jet_charged_truth_pt[ijet] < jet_pt_min) continue;
          if (jet_charged_truth_pt[ijet] > jet_pt_max) continue;
          if (abs(jet_charged_truth_eta[ijet]) > jet_eta_max) continue;

          // Observables: delta phi, jet pT, pT ratio
          Float_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - jet_charged_truth_phi[ijet]));
          Float_t jetpt = jet_charged_truth_pt[ijet];
          Float_t ptratio = jetpt / cluster_pt[n];

          corr[0] = centrality_v0m;
          corr[1] = cluster_pt[n];
          corr[2] = isolation;
          corr[3] = deltaphi;
          corr[4] = jetpt;
          corr[5] = ptratio;
          hCorrSRTruth->Fill(corr);
        }//for ijets      
      }//for nclusters
    } //for nevents
    file->Close();
    if (auxfile != NULL) {
      auxfile->Close();
    }
    std::cout << std::endl;
  }

  // Write to fout
  TFile* fout;
  fout = new TFile((TString) configrunperiod["filelists"]["correlations"]["gjmcchecks"].as<std::string>(), "RECREATE");
  std::cout << "Writing to file" << std::endl;

  hTrigPhoton->Write();
  hTrigSR->Write();
  hCorrSRAll->Write();
  hCorrSRTruth->Write();

  fout->Close();
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}