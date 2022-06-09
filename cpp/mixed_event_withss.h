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
THnSparseF* hCorrSR_jetmultcorrected;
THnSparseF* hCorrBR_jetmultcorrected;
THnSparseF* hCorr1ptSR_jetmultcorrected;
THnSparseF* hCorr1ptBR_jetmultcorrected;
int ndimTrig;
int ndimCorr;

// correlation variables
bool isSignal;
bool isBackground;
bool isIsolated;
float purity_weight;
float jetmultcorrection;

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
Average njet_ak02tpc per 1% centrality bin
Access by taking the floor of centrality_v0m as the index
First bin corresponds to 0-1%, second bin corresponds to 1-2%,
etc up to the 90th bin corresponding to 89-90%
These numbers come from 18q and 18r EMCEGA triggered data
and a few runs of 18q MB data
--------------------------------------------------------------*/
float avg_njet_ak02tpc_trig[90] = {150.2069209039548, 149.71067415730337, 149.04982078853047, 148.81341209173036, 148.69140919366993, 148.07617328519856, 147.40874811463047, 147.17260940032415, 146.7401215805471, 146.35220883534137, 146.00083333333333, 145.32545141874462, 145.09604286892002, 144.50717488789238, 144.1675799086758, 143.74697110904006, 143.38401559454192, 142.51731996353692, 142.049115913556, 141.5728476821192, 140.81186094069528, 140.2983193277311, 139.57709011943538, 139.29662921348316, 138.54446978335233, 138.053216374269, 137.09101941747574, 136.75531914893617, 136.17989417989418, 135.11316872427983, 134.26266280752532, 133.4362962962963, 132.60476190476192, 131.6465256797583, 130.65472312703582, 129.2687074829932, 128.73519736842104, 127.50535714285714, 126.68619246861925, 125.25901328273245, 123.68585858585858, 122.42963752665246, 121.32471910112359, 119.72622107969151, 117.99738219895288, 116.3388746803069, 114.77567567567567, 113.36153846153846, 111.12058823529412, 109.28521126760563, 108.0, 105.26344086021506, 103.56574394463668, 101.35140562248996, 99.08796296296296, 96.37012987012987, 94.92857142857143, 92.67171717171718, 90.55882352941177, 86.71348314606742, 84.1829268292683, 82.44771241830065, 80.50813008130082, 77.30152671755725, 74.75471698113208, 70.70183486238533, 68.76262626262626, 65.94554455445545, 63.65, 59.485074626865675, 57.524096385542165, 53.398305084745765, 50.72727272727273, 49.1551724137931, 46.77272727272727, 43.59230769230769, 40.7037037037037, 39.87837837837838, 33.07575757575758, 32.794117647058826, 31.77027027027027, 29.772727272727273, 26.5, 25.9, 23.863636363636363, 21.81578947368421, 21.863636363636363, 18.055555555555557, 20.75, 15.5};

float avg_njet_ak02tpc_mb[90] = {150.04477767726672, 149.4501287237955, 148.9766263237519, 148.56956115779644, 148.21263819295316, 147.80663361258334, 147.30022661901856, 146.98398479913138, 146.4794952681388, 146.00159559036842, 145.7273264984227, 145.24744793793386, 144.7735368956743, 144.31488485927244, 143.9772595091311, 143.41913617276543, 142.98884036875305, 142.58207987551867, 141.91298527443107, 141.42793931731984, 140.91765972720748, 140.4634405339806, 139.80893416927898, 139.1440097799511, 138.54287389047062, 137.78113941967445, 137.17164457718872, 136.3475633528265, 135.70076137046684, 134.90983939825168, 134.1922422954304, 133.37264472190694, 132.38012689780194, 131.49027635619242, 130.54155707522358, 129.44607710973307, 128.52622767857142, 127.29150107000918, 126.30380154639175, 124.90228384991843, 123.87035775127768, 122.50073126142595, 121.12922230950511, 119.65687086092716, 118.09295154185023, 116.7731524789523, 114.91275978733688, 113.20318021201413, 111.59320175438596, 109.72808022922636, 107.78496042216359, 105.5937764010138, 103.66492954223297, 101.5948753462604, 99.32982111079892, 96.98756616398681, 94.70995171110644, 92.1623941405356, 89.64820143884891, 87.23704826889058, 84.67640633533588, 81.94522864284664, 79.30636237897649, 76.55570161025304, 73.839982190561, 71.03628724216959, 68.25115984816533, 65.52235636038453, 62.6437227338827, 59.683373429564284, 57.245575863069334, 54.22776025236593, 51.66211604095563, 48.82128365916636, 46.17395706763872, 43.89307897071872, 41.52376333656644, 38.633616118769886, 36.597631426920856, 34.34990378447723, 32.24633431085044, 30.30849358974359, 28.588101725703904, 26.359556494192187, 25.142276422764226, 23.094101123595507, 21.771103896103895, 20.6, 18.059808612440193, 16.94034090909091};

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
