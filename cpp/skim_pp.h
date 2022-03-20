#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <THnSparse.h>

#define NTRACK_MAX (1U << 14)

/*--------------------------------------------------------------
Local variables
--------------------------------------------------------------*/
TFile *file;
TFile *auxfile;
TTree *_tree_event;
TTree *auxtree;
long nevents;

TTree *outtree;

/*--------------------------------------------------------------
Helper functions
--------------------------------------------------------------*/
void printCutSummary();
void openFilesAndGetTTrees(std::string root_filename);
void setBranchAddresses();

/*--------------------------------------------------------------
Variables from main tree
--------------------------------------------------------------*/
// Clusters
UInt_t ncluster;
Float_t cluster_pt[NTRACK_MAX];

/*--------------------------------------------------------------
Variables from aux tree
--------------------------------------------------------------*/
Float_t IN_cluster_5x5all[NTRACK_MAX];

// R=0.2 TPC jets
UInt_t IN_njet_ak02tpc;
Float_t IN_jet_ak02tpc_pt_raw[NTRACK_MAX];
Float_t IN_jet_ak02tpc_eta[NTRACK_MAX];
Float_t IN_jet_ak02tpc_phi[NTRACK_MAX];
Float_t IN_jet_ak02tpc_area[NTRACK_MAX];
UShort_t IN_jet_ak02tpc_multiplicity_raw[NTRACK_MAX];

// R=0.2 ITS jets
UInt_t IN_njet_ak02its;
Float_t IN_jet_ak02its_pt_raw[NTRACK_MAX];
Float_t IN_jet_ak02its_eta[NTRACK_MAX];
Float_t IN_jet_ak02its_phi[NTRACK_MAX];
Float_t IN_jet_ak02its_area[NTRACK_MAX];
UShort_t IN_jet_ak02its_multiplicity_raw[NTRACK_MAX];

// R=0.2 truth jets
UInt_t IN_njet_charged_truth_ak02;
Float_t IN_jet_charged_truth_ak02_pt_raw[NTRACK_MAX];
Float_t IN_jet_charged_truth_ak02_eta[NTRACK_MAX];
Float_t IN_jet_charged_truth_ak02_phi[NTRACK_MAX];
Float_t IN_jet_charged_truth_ak02_area[NTRACK_MAX];
UShort_t IN_jet_charged_truth_ak02_multiplicity_raw[NTRACK_MAX];

/*--------------------------------------------------------------
Variables to be written to output tree
--------------------------------------------------------------*/
UInt_t outncluster;
Float_t cluster_5x5all[NTRACK_MAX];

// R=0.2 TPC jets
UInt_t njet_ak02tpc;
Float_t jet_ak02tpc_pt_raw[NTRACK_MAX];
Float_t jet_ak02tpc_eta[NTRACK_MAX];
Float_t jet_ak02tpc_phi[NTRACK_MAX];
Float_t jet_ak02tpc_area[NTRACK_MAX];
UShort_t jet_ak02tpc_multiplicity_raw[NTRACK_MAX];

// R=0.2 ITS jets
UInt_t njet_ak02its;
Float_t jet_ak02its_pt_raw[NTRACK_MAX];
Float_t jet_ak02its_eta[NTRACK_MAX];
Float_t jet_ak02its_phi[NTRACK_MAX];
Float_t jet_ak02its_area[NTRACK_MAX];
UShort_t jet_ak02its_multiplicity_raw[NTRACK_MAX];

// R=0.2 truth jets
UInt_t njet_charged_truth_ak02;
Float_t jet_charged_truth_ak02_pt_raw[NTRACK_MAX];
Float_t jet_charged_truth_ak02_eta[NTRACK_MAX];
Float_t jet_charged_truth_ak02_phi[NTRACK_MAX];
Float_t jet_charged_truth_ak02_area[NTRACK_MAX];
UShort_t jet_charged_truth_ak02_multiplicity_raw[NTRACK_MAX];