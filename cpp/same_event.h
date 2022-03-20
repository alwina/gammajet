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
Long64_t nevents_max;
long nevents;

// THnSparses
THnSparseF* hTrigSR;
THnSparseF* hTrigBR;
THnSparseF* hTrigSRJet;
THnSparseF* hTrigBRJet;
THnSparseF* hTrigSRJetB2b;
THnSparseF* hTrigBRJetB2b;
THnSparseF* hCorrSR;
THnSparseF* hCorrBR;
THnSparseF* hCorr1ptSR;
THnSparseF* hCorr1ptBR;
int ndimTrig;
int ndimCorr;

// correlation variables
bool isSignal;
bool isBackground;
bool isIsolated;
float purity_weight;

/*--------------------------------------------------------------
Helper functions
--------------------------------------------------------------*/
void printCutSummary();
void initializeTHnSparses();
void openFilesAndGetTTrees(std::string root_filename);
void setBranchAddresses();
bool rejectCluster(int icluster);
float getIsolation(int icluster);
float getShower(int icluster);

/*--------------------------------------------------------------
Variables from TTrees
--------------------------------------------------------------*/
// Events
Bool_t is_pileup_from_spd_5_08;
Double_t primary_vertex[3];
Float_t ue_estimate_its_const;
Float_t ue_estimate_tpc_const;
Float_t centrality_v0m;

// Tracks
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

// Clusters
UInt_t ncluster;
Float_t cluster_e[NTRACK_MAX];
Float_t cluster_e_max[NTRACK_MAX];
Float_t cluster_e_cross[NTRACK_MAX];
Float_t cluster_pt[NTRACK_MAX];
Float_t cluster_eta[NTRACK_MAX];
Float_t cluster_phi[NTRACK_MAX];
Float_t cluster_iso_tpc_02[NTRACK_MAX];
Float_t cluster_iso_tpc_04[NTRACK_MAX];
Float_t cluster_iso_its_02[NTRACK_MAX];
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
Float_t cluster_5x5all[NTRACK_MAX];

Float_t cluster_tof[NTRACK_MAX];
Float_t cluster_iso_its_02_ue[NTRACK_MAX];
Float_t cluster_iso_its_04_ue[NTRACK_MAX];
Float_t cluster_iso_tpc_02_ue[NTRACK_MAX];
Float_t cluster_iso_tpc_04_ue[NTRACK_MAX];

// Jets
UInt_t njet;
Float_t jet_pt_raw[NTRACK_MAX];
Float_t jet_eta[NTRACK_MAX];
Float_t jet_phi[NTRACK_MAX];
Float_t jet_area[NTRACK_MAX];
UShort_t jet_multiplicity_raw[NTRACK_MAX];

// MC
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