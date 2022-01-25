#include <THStack.h>
#include <THnSparse.h>

using namespace H5;

/*--------------------------------------------------------------
Local variables
--------------------------------------------------------------*/
long nevents_max;
bool doprint;
std::string mixlabel;
int mix_start;
int mix_end;
int mix_index;

// THnSparses
THnSparseF* hTrigSR;
THnSparseF* hTrigBR;
THnSparseF* hnMixSR;
THnSparseF* hnMixBR;
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

// debugging histograms
TH1F *hClusterpT;
TH1F *hClusterLambda;
TH1F *hClustereta;
TH1F *hClusterIso;
TH1I ** hnJets;

// book-keeping histograms
TH1I *SR_mixed_event_counter;
TH1I *BR_mixed_event_counter;
TH1D* z_vertices_triggered;
TH1D* z_vertices_MinBias;
TH1D* delta_z_vertices;
TH1D* multiplicity_triggered;
TH1D* multiplicity_MinBias;
TH1D* delta_multiplicity;
TH1D* flow_triggered;
TH1D* flow_MinBias;
TH1D* delta_flow;
TH1D* centrality_triggered;
TH1D* centrality_MinBias;
TH1D* delta_centrality;

/*--------------------------------------------------------------
Helper functions
--------------------------------------------------------------*/
void parseInputs(int argc, char *argv[]);
void printCutSummary();
void initializeTHnSparses();
void initializeDebuggingHistograms();

void setTriggeredEventVariables(float event_data_values[]);
void setMBEventVariables(float mb_event_data_values[]);
void setClusterVariables(float cluster_data_values[]);
void setTrackVariables(float track_data_values[]);
void setJetVariables(float jet_data_values[]);

bool rejectCluster();
float getTriggeredIsolation();
float getShower();

/*--------------------------------------------------------------
Variables from HDF5 files
--------------------------------------------------------------*/
// triggered event variables
float primary_vertex; // z-vertex
float multiplicity;
float v2;
float centrality_v0m;
bool  is_pileup_from_spd_5_08;
float ue_estimate_its_const;
float ue_estimate_tpc_const;

// MB event variables
float mb_primary_vertex;
float mb_multiplicity;
float mb_v2;
float mb_centrality_v0m;
float mb_ue_estimate_its_const;
float mb_ue_estimate_tpc_const;

// cluster variables
float cluster_e;
float cluster_pt;
float cluster_eta;
float cluster_phi;
float cluster_lambda_square;
float cluster_e_max;
float cluster_e_cross;
float cluster_iso_tpc_02;
float cluster_iso_tpc_03;
float cluster_iso_tpc_04;
float cluster_iso_its_02;
float cluster_iso_its_03;
float cluster_iso_its_04;
float cluster_frixione_tpc_04_02;
float cluster_frixione_its_04_02;
float cluster_s_nphoton0;
float cluster_mc_truth_index;
float cluster_ncell;
float cluster_cell_id_max;
float cluster_cell_e;
float cluster_distance_to_bad_channel;
float cluster_nlocal_maxima;
float cluster_tof;
float cluster_iso_its_02_ue;
float cluster_iso_its_03_ue;
float cluster_iso_its_04_ue;
float cluster_iso_tpc_02_ue;
float cluster_iso_tpc_03_ue;
float cluster_iso_tpc_04_ue;
float cluster_lambda1;
float cluster_s_nphoton;
float cluster_5x5all;

// track variables
float track_pt;
float track_eta;
float track_phi;

// jet variables
float jet_pt_raw;
float jet_eta;
float jet_phi;
