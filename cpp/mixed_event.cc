#include <fstream>
#include <iostream>
#include <math.h>

#include <TFile.h>
#include <TLorentzVector.h>
#include "H5Cpp.h"

#include "mixed_event.h"
#include "config_parser.h"
#include "shared_defs.h"

using namespace H5;

int main(int argc, char *argv[])
{

	/*--------------------------------------------------------------
	Initial setup
	--------------------------------------------------------------*/
	parseInputs(argc, argv);
	YAML::Node configrunperiod = YAML::LoadFile(argv[1]);
	allconfigs.push_back(configrunperiod);
	YAML::Node configsystem = YAML::LoadFile(configrunperiod["systemconfig"].as<std::string>());
	allconfigs.push_back(configsystem);
	YAML::Node configglobal = YAML::LoadFile(configsystem["globalconfig"].as<std::string>());
	allconfigs.push_back(configglobal);
	parseConfig();
	if (doprint) printCutSummary();
	if (doprint) std::cout << "Mix label: " << mixlabel << std::endl;

	// set up THnSparses
	initializeDebuggingHistograms();
	initializeTHnSparses();
	Double_t trig[ndimTrig];
	Double_t corr[ndimCorr];

	// get filenames and read HDF5 files
	std::string triggered_hdf5_file_name = configrunperiod["filelists"]["mixing"][mixlabel]["triggered"].as<std::string>();
	if (doprint) std::cout << triggered_hdf5_file_name << std::endl;

	std::string MBhdf5_file_name = configrunperiod["filelists"]["mixing"][mixlabel]["mb"].as<std::string>();
	if (doprint) std::cout << MBhdf5_file_name << std::endl;

	std::string pairing_filename = configrunperiod["filelists"]["mixing"][mixlabel]["pairing"].as<std::string>();
	if (doprint) std::cout << pairing_filename << std::endl;

	/*---------------------------------------------------------------
	Read HDF5 files
	I tried to pull this into its own function, but the limitations
	of forward-declaring arrays are beyond my ability to handle.

	Using low level hdf5 API for data

	Procedure goes something like this:
	1. Get Datasets from the appropriate triggered or MB hdf5 file
	2. Get Dataspaces from datasets, and record dimensions locally
	3. Define size and offset of hyperslab to be read from DISK
	4. Define size and offset of hyperslab to read into MEMORY
	5. Define local arrays to store information
	6. DISK hyperslab -> MEMORY hyperslab -> local arrays
	7. Iterate through hdf5 file by incrementing offset

	Terms:
	hyperslab: an n-dimensional "chunk" of the hdf5 data (disk or memory)
	offset: the n-dimensional starting point of the hyperslab.
	dataspace: defines size and shape of dataset or attribute

	Documentation:
	https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html
	Instead of ranks as strict definitions (RANK_OUT for example),
	I use dimensions obtained by reading dataspace info the hdf5 file

	This is nearly unchanged from what was written by ftorales
	---------------------------------------------------------------*/
	// open files
	H5File triggered_h5_file((H5std_string) triggered_hdf5_file_name, H5F_ACC_RDONLY );
	H5File MB_h5_file((H5std_string) MBhdf5_file_name, H5F_ACC_RDONLY );

	// get DataSet names
	const H5std_string cluster_ds_name("cluster");
	const H5std_string event_ds_name("event");
	const H5std_string track_ds_name("track");

	H5std_string jet_ds_name;
	if (jettype == "ak04tpc" or jettype == "ak02tpc") {
		jet_ds_name = "jet_" + jettype;
	} else if (jettype == "none") {
		jet_ds_name = "jet";
	} else if (jettype == "ak04its") {
		std::cout << "ERROR: Jet type ak04its not yet supported. Aborting" << std::endl;
		exit(EXIT_FAILURE);
	} else {
		std::cout << "ERROR: Jet type " << jettype << " not recognized. Aborting" << std::endl;
		exit(EXIT_FAILURE);
	}

	// get DataSets
	DataSet cluster_dataset = triggered_h5_file.openDataSet(cluster_ds_name);
	DataSet event_dataset = triggered_h5_file.openDataSet(event_ds_name);
	DataSet jet_dataset = MB_h5_file.openDataSet(jet_ds_name);
	DataSet track_dataset = MB_h5_file.openDataSet(track_ds_name);
	DataSet mb_event_dataset = MB_h5_file.openDataSet(event_ds_name);

	// get DataSpaces
	DataSpace cluster_dataspace = cluster_dataset.getSpace();
	DataSpace event_dataspace = event_dataset.getSpace();
	DataSpace jet_dataspace = jet_dataset.getSpace();
	DataSpace track_dataspace = track_dataset.getSpace();
	DataSpace mb_event_dataspace = mb_event_dataset.getSpace();

	// Load the dimensions of datasets from file, to be used in dataspace/hyperslab
	// first get the number of dimensions with ExtentNdims
	// Then get the length of each dimension with ExtentDims.
	const int cluster_rank = cluster_dataspace.getSimpleExtentNdims();
	hsize_t clusterdims[cluster_rank];
	cluster_dataspace.getSimpleExtentDims(clusterdims, NULL);
	UInt_t ncluster_max = clusterdims[1];
	UInt_t Ncluster_Vars = clusterdims[2];

	const int event_rank = event_dataspace.getSimpleExtentNdims();
	hsize_t eventdims[event_rank];
	event_dataspace.getSimpleExtentDims(eventdims, NULL);
	long nevent = eventdims[0];
	UInt_t Nevent_Vars = eventdims[1];

	const int jet_rank = jet_dataspace.getSimpleExtentNdims();
	hsize_t jetdims[jet_rank];
	jet_dataspace.getSimpleExtentDims(jetdims, NULL);
	UInt_t njet = jetdims[1];
	UInt_t Njet_Vars = jetdims[2];

	const int track_rank = track_dataspace.getSimpleExtentNdims();
	hsize_t trackdims[track_rank];
	track_dataspace.getSimpleExtentDims(trackdims, NULL);
	UInt_t ntrack = trackdims[1];
	UInt_t Ntrack_Vars = trackdims[2];

	const int mb_event_rank = mb_event_dataspace.getSimpleExtentNdims();
	hsize_t mb_eventdims[mb_event_rank];
	mb_event_dataspace.getSimpleExtentDims(mb_eventdims, NULL);

	// Block size should be determined by chunk size in to_hdf5. Usually 2000 or its multiples.
	// A larger block size will speed things up, at the cost of more memory
	const int block_size = 2000;
	const hsize_t nmix_max = 300;
	hsize_t nmix = std::min(nmix_max, mb_eventdims[0] / block_size);

	// Local arrays the hyperslabs will be fed into, and ultimately used in loops
	float cluster_data_out[block_size][ncluster_max][Ncluster_Vars];
	float event_data_out[block_size][Nevent_Vars];
	float jet_data_out[block_size][njet][Njet_Vars];
	float track_data_out[block_size][ntrack][Ntrack_Vars];
	float mb_event_data_out[block_size][Nevent_Vars];

	//Define hyperslab size and offset to be read from file
	hsize_t event_offset[2] = {0, 0};
	hsize_t cluster_offset[3] = {0, 0, 0};
	hsize_t jet_offset[3] = {0, 0, 0};
	hsize_t track_offset[3] = {0, 0, 0};
	hsize_t mb_event_offset[2] = {0, 0};

	hsize_t event_count[2] = {block_size, Nevent_Vars};
	hsize_t cluster_count[3] = {block_size, ncluster_max, Ncluster_Vars};
	hsize_t jet_count[3] = {block_size, njet, Njet_Vars};
	hsize_t track_count[3] = {block_size, ntrack, Ntrack_Vars};
	hsize_t mb_event_count[2] = {block_size, Nevent_Vars};

	// The offset is how we iterate over the entire hdf5 file.
	// For example, To obtain data for event 68, set the
	// offset to {68, njet, Njet_Vars}.
	event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
	cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
	jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
	track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
	mb_event_dataspace.selectHyperslab( H5S_SELECT_SET, mb_event_count, mb_event_offset );

	// Define the memory dataspace in which to place hyperslab
	DataSpace event_memspace(event_rank, eventdims);
	DataSpace cluster_memspace(cluster_rank, clusterdims);
	DataSpace jet_memspace(jet_rank, jetdims);
	DataSpace track_memspace(track_rank, trackdims);
	DataSpace mb_event_memspace(mb_event_rank, mb_eventdims);

	// Define DataSpace offset for hyperslab starting at beginning:
	// We want to read the entire "chunk" of the DataSpace to an array
	// So we set the DataSpace offset to [0]. Can be important if low memory
	hsize_t event_offset_out[2] = {0};
	hsize_t cluster_offset_out[3] = {0};
	hsize_t jet_offset_out[3] = {0};
	hsize_t track_offset_out[3] = {0};
	hsize_t mb_event_offset_out[2] = {0};

	// Define dimensions of hyperslab in memory (aka memspace)
	hsize_t event_count_out[2] = {block_size, Nevent_Vars};
	hsize_t cluster_count_out[3] = {block_size, ncluster_max, Ncluster_Vars};
	hsize_t jet_count_out[3] = {block_size, njet, Njet_Vars};
	hsize_t track_count_out[3] = {block_size, ntrack, Ntrack_Vars};
	hsize_t mb_event_count_out[2] = {block_size, Nevent_Vars};

	// Apply the offset and dimensions from the previous two code blocks to the memspace
	event_memspace.selectHyperslab( H5S_SELECT_SET, event_count_out, event_offset_out );
	cluster_memspace.selectHyperslab( H5S_SELECT_SET, cluster_count_out, cluster_offset_out );
	jet_memspace.selectHyperslab( H5S_SELECT_SET, jet_count_out, jet_offset_out );
	track_memspace.selectHyperslab( H5S_SELECT_SET, track_count_out, track_offset_out );
	mb_event_memspace.selectHyperslab( H5S_SELECT_SET, mb_event_count_out, mb_event_offset_out );

	// FINALLY use the well-defined memspace to read data from the dataspace, INTO the local array
	event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );
	cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
	jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );
	track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
	mb_event_dataset.read( mb_event_data_out, PredType::NATIVE_FLOAT, mb_event_memspace, mb_event_dataspace );

	// First [block_size] number of events have just been read into local arrays

	/*--------------------------------------------------------------
	Main correlation loop
	Loop through number of events to mix with
	--------------------------------------------------------------*/
	for (Long64_t imix = mix_start; imix < mix_end; imix++) {
		if (imix > nmix) break;

		// Grab 2000 Mixed Events at a time
		// Takes advantage of block structure used in pairing
		if (imix * block_size < mb_eventdims[0] - block_size - 1) {
			mb_event_offset[0] = imix * block_size;
			mb_event_dataspace.selectHyperslab( H5S_SELECT_SET, mb_event_count, mb_event_offset );
			mb_event_dataset.read( mb_event_data_out, PredType::NATIVE_FLOAT, mb_event_memspace, mb_event_dataspace );

			//FIXME: experiment with read outside of cluster loop. Would mean htrig histos need re-working as well
			jet_offset[0] = imix * block_size;
			jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
			jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );

			track_offset[0] = imix * block_size;
			track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
			track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
		}

		// open pairing file and prepare to loop through it
		std::ifstream pairing_textfile;
		pairing_textfile.open(pairing_filename);

		/*--------------------------------------------------------------
		Loop through triggered events
		--------------------------------------------------------------*/
		int offset = 0; //Offset for Triggered Events
		nevent = std::min(nevent, nevents_max);
		for (Long64_t ievent = 0; ievent < nevent; ievent++) {
			if (doprint) fprintf(stderr, "\r%s:%d: mix %llu / %llu, event %llu / %llu", __FILE__, __LINE__, imix - mix_start, mix_end - mix_start, ievent, nevent);
			// reading from file is done every [block_size] number of events
			// this variable keeps track of the current increment within a block
			// as opposed to [ievent] which is looping through all events
			int ieventinblock = ievent % block_size;

			if ((ieventinblock == (block_size - 1)) && (ievent != nevent - 1) && (ievent < nevent - block_size - 1)) {
				// load 1 block (2000 events) at a time. Faster/less memory
				offset += block_size;

				event_offset[0] = offset;
				event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
				event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );

				cluster_offset[0] = offset;
				cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
				cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
			}

			setTriggeredEventVariables(event_data_out[ieventinblock]);

			// event selection
			if (abs(primary_vertex) > 10) continue;
			if (primary_vertex == 0.00) continue;
			if (do_pile && is_pileup_from_spd_5_08) continue;

			// fill triggered event debug histograms
			z_vertices_triggered->Fill(primary_vertex);
			flow_triggered->Fill(v2);
			centrality_triggered->Fill(centrality_v0m);
			multiplicity_triggered->Fill(multiplicity);

			// parse pairing file to get event number to mix with
			Long64_t mix_event = -1;
			std::string eventline;
			if (ievent > 0) {
				// skips \n that separates each triggered event's pairings list
				getline(pairing_textfile, eventline);
			}
			getline(pairing_textfile, eventline);

			if (eventline.size() == 0) {
				mix_event = -999;
			} else {
				std::string mixednum_string;
				std::istringstream parser[1];
				parser[0].str(eventline);
				// skip until we get to the mix number we're on
				// <= because imix starts at 0
				for (int jmix = 0; jmix <= imix; jmix++) {
					getline(parser[0], mixednum_string, '\t');
				}
				// account for any blanks (once imix gets too large)
				if (mixednum_string.size() == 0) {
					mix_event = -999;
				} else {
					mix_event = stoul(mixednum_string);
				}
			}

			if (mix_event < 0) continue; // unpaired events set to negative numbers
			mix_index = mix_event % block_size; // get the relevant index for this specific block
			setMBEventVariables(mb_event_data_out[mix_index]);

			// fill mixing debug histograms
			z_vertices_MinBias->Fill(mb_primary_vertex);
			flow_MinBias->Fill(mb_v2);
			multiplicity_MinBias->Fill(mb_multiplicity);
			centrality_MinBias->Fill(mb_centrality_v0m);

			delta_z_vertices->Fill(abs(primary_vertex - mb_primary_vertex));
			delta_flow->Fill(abs(v2 - mb_v2));
			delta_multiplicity->Fill(abs(multiplicity - mb_multiplicity));
			delta_centrality->Fill(abs(centrality_v0m - mb_centrality_v0m));

			// event pairing cuts - move this earlier so we can skip the cluster loop
			if (abs(centrality_v0m - mb_centrality_v0m) > 10.) continue;
			if (abs(primary_vertex - mb_primary_vertex) > 2.) continue;
			if (abs(v2 - mb_v2) > 0.5) continue;

			/*--------------------------------------------------------------
			Loop through clusters in triggered event
			--------------------------------------------------------------*/
			for (ULong64_t icluster = 0; icluster < ncluster_max; icluster++) {
				if (std::isnan(cluster_data_out[ieventinblock][icluster][0])) break;
				setClusterVariables(cluster_data_out[ieventinblock][icluster]);

				// fill cluster debug histograms
				hClusterpT->Fill(cluster_pt);
				hClusterLambda->Fill(cluster_lambda_square);
				hClustereta->Fill(cluster_eta);
				hClusterIso->Fill(cluster_iso_tpc_04_ue);

				// apply cluster cuts
				if (rejectCluster()) continue;
				
				// determine whether the cluster is isolated in the triggered event
				// otherwise, we're looking at a different sample of clusters
				float triggeredIsolation = getTriggeredIsolation();
				bool isTrigIsolated = GetIsIsolated(triggeredIsolation, centrality_v0m, isoconfig);
				if (not(isTrigIsolated)) continue;

				// determine whether the cluster is isolated in the MB event
				// calculate isolation based on MB tracks
				// for this, it only matters what the radius is, not the track type
				// since the making of the HDF5 file keeps track of which track type it wants
				// since we need the track array, we don't separate out getIsolation()
				float cone;
				float ue_estimate;
				if (isovar == "cluster_iso_tpc_04_sub") {
					cone = 0.4;
					ue_estimate = mb_ue_estimate_tpc_const;
				} else if (isovar == "cluster_iso_its_04_sub") {
					cone = 0.4;
					ue_estimate = mb_ue_estimate_its_const;
				} else if (isovar == "cluster_iso_tpc_02_sub") {
					cone = 0.2;
					ue_estimate = mb_ue_estimate_tpc_const;
				} else {
					std::cout << "ERROR: Isolation variable " << isovar << " not currently calculable from MB event. Aborting" << std::endl;
					exit(EXIT_FAILURE);
				}

				// start with the UE subtraction, because why not?
				float isolation = -M_PI * cone * cone * ue_estimate;
				for (int itrack = 0; itrack < ntrack; itrack++) {
					if (std::isnan(track_data_out[mix_index][itrack][0])) break;
					setTrackVariables(track_data_out[mix_index][itrack]);

					// if we're within the cone, add the track pT to the isolation energy
					float dist2 = (cluster_phi - track_phi) * (cluster_phi - track_phi) + (cluster_eta - track_eta) * (cluster_eta - track_eta);
					if (dist2 < (cone * cone)) {
						isolation += track_pt;
					}
				}

				isIsolated = GetIsIsolated(isolation, centrality_v0m, isoconfig);
				if (not(isIsolated)) continue;

				// determine whether it is SR or BR (or neither), calculate purity, and fill trigger THnSparse
				float shower = getShower();
				isSignal = (shower > srmin) and (shower < srmax);
				isBackground = (shower > brmin) and (shower < brmax);
				if (not(isSignal or isBackground)) continue;

				float purity = getPurity(cluster_pt, centrality_v0m, purityconfig);

				trig[0] = centrality_v0m;
				trig[1] = cluster_pt;

				// count number of cluster-ME pairs
				if (isSignal) {
					purity_weight = 1.0 / purity;
					hnMixSR->Fill(trig);
					SR_mixed_event_counter->Fill(mix_event);
					if (imix == 0) {
						// this should cause hTrigSR to be filled at most once
						// and therefore be a pretty good estimate of how many clusters were paired
						hTrigSR->Fill(trig);
					}
				}

				if (isBackground) {
					purity_weight = 1.0 / purity - 1;
					hnMixBR->Fill(trig);
					BR_mixed_event_counter->Fill(mix_event);
					if (imix == 0) {
						// this should cause hTrigBR to be filled at most once
						// and therefore be a pretty good estimate of how many clusters were paired
						hTrigBR->Fill(trig);
					}
				}

				/*--------------------------------------------------------------
				Loop through jets
				--------------------------------------------------------------*/
				for (ULong64_t ijet = 0; ijet < njet; ijet++) {
					if (std::isnan(jet_data_out[mix_index][ijet][0])) break;
					setJetVariables(jet_data_out[mix_index][ijet]);

					// apply jet cuts
					if (not(jet_pt_raw > jet_pt_min and jet_pt_raw < jet_pt_max)) continue;
					if (not(abs(jet_eta) < jet_eta_max)) continue;

					// Observables: delta phi, jet pT, pT ratio
					Float_t deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi - jet_phi));
					Float_t jetpt = jet_pt_raw;
					Float_t ptratio = jetpt / cluster_pt;

					corr[0] = centrality_v0m;
					corr[1] = cluster_pt;
					corr[2] = deltaphi;
					corr[3] = jetpt;
					corr[4] = ptratio;

					if (isSignal) {
						hCorrSR->Fill(corr, purity_weight);
						hCorr1ptSR->Fill(corr, purity_weight / jetpt);
					}

					if (isBackground) {
						hCorrBR->Fill(corr, purity_weight);
						hCorr1ptBR->Fill(corr, purity_weight / jetpt);
					}
				} // end jet loop
			} // end cluster loop
		} // end triggered event loop
		pairing_textfile.close();
	} // end number of mixed events loop

	/*--------------------------------------------------------------
	Write outputs to file
	--------------------------------------------------------------*/
	TFile* fout;
	/* fout = new TFile((TString) configrunperiod["filelists"]["correlations"]["mixedevent"].as<std::string>(), "RECREATE"); */
	std::string basename = configrunperiod["filelists"]["correlations"]["mixedevent"].as<std::string>();
	std::string extendedname = basename.replace(basename.find(".root"), 5, "_" + mixlabel + ".root");
	if (argc > 3) {
		fout = new TFile((TString) extendedname + Form("mix_%i.root", mix_start), "RECREATE");
	} else {
		fout = new TFile((TString) extendedname, "RECREATE");
	}
	/* std::cout << "Writing to file" << std::endl; */

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

	for (int h = 0; h < 10; h++) {
		hnJets[h]->Write();
	}

	fout->Close();
	if (doprint) std::cout << std::endl << "Ending" << std::endl;

	return EXIT_SUCCESS;
}

void parseInputs(int argc, char *argv[])
{
	if (argc < 3 || argc > 6) {
		fprintf(stderr, "Format: [command] [config file] [mixlabel] [mix start index (default=0)] [mix end index (default=10)] [nevents (optional)]\n");
		exit(EXIT_FAILURE);
	}

	mixlabel = argv[2];

	if (argc > 3) {
		if (argc < 5) {
			fprintf(stderr, "Must define both mix start index and mix end index\n");
			exit(EXIT_FAILURE);
		} else {
			mix_start = atoi(argv[3]);
			mix_end = atoi(argv[4]);
		}
	} else {
		mix_start = 0;
		mix_end = 10;
	}

	if (argc > 5) {
		nevents_max = std::stol(argv[5]);
	} else {
		nevents_max = 999999999999999;
	}

	// try to only print once
	doprint = (mix_start == 0);
	// std::cout << mix_start << "-" << mix_end << std::endl;
}

// Print cut summary
void printCutSummary()
{
	std::cout << "Cluster pT range: " << cluster_pt_min << "-" << cluster_pt_max << std::endl;
	std::cout << "Cluster eta max: " << cluster_eta_max << std::endl;
	std::cout << "Cluster ncell min: " << cluster_ncell_min << std::endl;
	std::cout << "Cluster Ecross/Emax min: " << cluster_ecross_emax_min << std::endl;
	std::cout << "Cluster dist to bad channel min: " << cluster_dbc_min << std::endl;
	std::cout << "Cluster nlocal maxima max: " << cluster_nlm_max << std::endl;
	std::cout << "Cluster TOF max: " << cluster_tof_max << std::endl;
	std::cout << "Shower shape SR range: " << srmin << "-" << srmax << std::endl;
	std::cout << "Shower shape BR range: " << brmin << "-" << brmax << std::endl;
	std::cout << "Jet type: " << jettype << std::endl;
	std::cout << "Jet pT range: " << jet_pt_min << "-" << jet_pt_max << std::endl;
	std::cout << "Jet eta max: " << jet_eta_max << std::endl;
}

void initializeTHnSparses()
{
	/*--------------------------------------------------------------
	Setting up THnSparses
	hCorrSR: cluster-jet correlations for the signal region
	hTrigSR: counting the number of clusters in each bin in the signal region
	hnMixSR: counting the number of cluster-ME pairs in the signal region
	hCorrBR: cluster-jet correlations for the bkg region
	hTrigBR: counting the number of clusters in each bin in the bkg region
	hnMixBR: counting the number of cluster-ME pairs in the bkg region
	--------------------------------------------------------------*/
	int nbinsClusterPt = 2 * (cluster_pt_max - cluster_pt_min);

	// dimensions: centrality, cluster pT
	ndimTrig = 2;
	Int_t nbinsTrig[ndimTrig];
	Double_t minbinsTrig[ndimTrig];
	Double_t maxbinsTrig[ndimTrig];

	// centrality
	nbinsTrig[0] = 10;
	minbinsTrig[0] = 0;
	maxbinsTrig[0] = 100;

	// cluster pT
	nbinsTrig[1] = nbinsClusterPt;
	minbinsTrig[1] = cluster_pt_min;
	maxbinsTrig[1] = cluster_pt_max;

	hTrigSR = new THnSparseF("hTrigSR", "Mixed Event Number of clusters (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
	hTrigBR = new THnSparseF("hTrigBR", "Mixed Event Number of clusters (BR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

	hnMixSR = new THnSparseF("hnMixSR", "Number of cluster-ME pairs (SR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);
	hnMixBR = new THnSparseF("hnMixBR", "Number of cluster-ME pairs (BR)", ndimTrig, nbinsTrig, minbinsTrig, maxbinsTrig);

	// dimensions: centrality, cluster pT, delta phi, jet pT, pT ratio
	ndimCorr = 5;
	Int_t nbinsCorr[ndimCorr];
	Double_t minbinsCorr[ndimCorr];
	Double_t maxbinsCorr[ndimCorr];

	// centrality
	nbinsCorr[0] = 10;
	minbinsCorr[0] = 0;
	maxbinsCorr[0] = 100;

	// cluster pT
	nbinsCorr[1] = nbinsClusterPt;
	minbinsCorr[1] = cluster_pt_min;
	maxbinsCorr[1] = cluster_pt_max;

	// deltaphi
	nbinsCorr[2] = 120;
	minbinsCorr[2] = deltaphi_min;
	maxbinsCorr[2] = deltaphi_max;

	// jetpt
	nbinsCorr[3] = 120;
	minbinsCorr[3] = jetpt_min;
	maxbinsCorr[3] = jetpt_max;

	// ptratio
	nbinsCorr[4] = 120;
	minbinsCorr[4] = ptratio_min;
	maxbinsCorr[4] = ptratio_max;

	hCorrSR = new THnSparseF("hCorrSR", "ME Correlations (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrBR = new THnSparseF("hCorrBR", "ME Correlations (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorr1ptSR = new THnSparseF("hCorr1ptSR", "ME Correlations with 1/jetpt weight (SR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorr1ptBR = new THnSparseF("hCorr1ptBR", "ME Correlations with 1/jetpt weight (BR)", ndimCorr, nbinsCorr, minbinsCorr, maxbinsCorr);
	hCorrSR->Sumw2();
	hCorrBR->Sumw2();
	hCorr1ptSR->Sumw2();
	hCorr1ptBR->Sumw2();
}

void initializeDebuggingHistograms()
{
	hClusterpT = new TH1F("hClusterpT", "Cluster pT Distribution", 200, 0, 100);
	hClusterLambda = new TH1F("hClusterLambda", "Cluster Lambda Distribution", 40, 0, 1.5);
	hClustereta = new TH1F("hClustereta", "Cluster eta Distribution", 60, -3, 3);
	hClusterIso = new TH1F("hClusterIso", "Cluster Iso Distribution", 120, -10, 50);
	/* TH1I *hnJets = new TH1I("hnJets","Number of jets that pass cuts",20,0,20); */

	hnJets = new TH1I*[10];
	for (int h = 0; h < 10; h++)
		hnJets[h] = new TH1I(Form("hnJetMix_%i", h), Form("Number of Jets in mixed event %i", h), 20, 0, 20);

	/*---------------------------------------------------------------
		Mixed Event Book Keeping
		---------------------------------------------------------------*/
	SR_mixed_event_counter = new TH1I("SR_mixed_event_counter", "Distribution of Mixed Event No.", 1E6, 0, 1E6);
	BR_mixed_event_counter = new TH1I("BR_mixed_event_counter", "Distribution of Mixed Event No.", 1E6, 0, 1E6);
	//GetEntries() of above is all that is needed. Distributions serve as additional check to mixing

	z_vertices_triggered = new TH1D("Primary_Vertex_triggered", "Z-vertex (triggered)", 240, -12, 12);
	z_vertices_MinBias = new TH1D("Primary_Vertex_MinBias", "Z-vertex (MinBias)", 240, -12, 12);
	delta_z_vertices = new TH1D("Delta_Primary_Vertex", "#Delta V_z Distribution", 240, -12, 12);

	multiplicity_triggered = new TH1D("multiplicity_triggered", "multiplicity (triggered)", 1000, 0, 1000);
	multiplicity_MinBias = new TH1D("Multplicity_MinBias", "multiplicity (MinBias)", 500, 0, 1000);
	delta_multiplicity = new TH1D("Delta_multiplicity", "#Delta Multiplicit Distribution", 500, 0, 1000);

	flow_triggered = new TH1D("Flow_triggered", "Flow (triggered)", 500, -2, 2);
	flow_MinBias = new TH1D("Flow_MinBias", "Flow (MinBias)", 500, -2, 2);
	delta_flow = new TH1D("Delta_Flow", "#Delta Flow Distribution", 500, 0, 4);

	centrality_triggered = new TH1D("centrality_triggered", "centrality (triggered)", 200, 0, 100);
	centrality_MinBias = new TH1D("centrality_MinBias", "centrality (MinBias)", 200, 0, 100);
	delta_centrality = new TH1D("Delta_centrality", "#Delta centrality Distribution", 200, 0, 100);
}

void setTriggeredEventVariables(float event_data_values[])
{
	primary_vertex = event_data_values[0];
	multiplicity = event_data_values[1];
	v2 = event_data_values[2];
	centrality_v0m = event_data_values[3];
	is_pileup_from_spd_5_08 = event_data_values[4];
	ue_estimate_its_const = event_data_values[5];
	ue_estimate_tpc_const = event_data_values[6];
}

void setMBEventVariables(float mb_event_data_values[])
{
	mb_primary_vertex = mb_event_data_values[0];
	mb_multiplicity = mb_event_data_values[1];
	mb_v2 = mb_event_data_values[2];
	mb_centrality_v0m = mb_event_data_values[3];
	mb_ue_estimate_its_const = mb_event_data_values[5];
	mb_ue_estimate_tpc_const = mb_event_data_values[6];
}

void setClusterVariables(float cluster_data_values[])
{
	//See to_hdf5.cc for ROOT->HDF5 details
	cluster_e = cluster_data_values[0];
	cluster_pt = cluster_data_values[1];
	cluster_eta = cluster_data_values[2];
	cluster_phi = cluster_data_values[3];
	cluster_lambda_square = cluster_data_values[4];
	cluster_e_max = cluster_data_values[5];
	cluster_e_cross = cluster_data_values[6];
	cluster_iso_tpc_02 = cluster_data_values[7];
	cluster_iso_tpc_03 = cluster_data_values[8];
	cluster_iso_tpc_04 = cluster_data_values[9];
	cluster_iso_its_02 = cluster_data_values[10];
	cluster_iso_its_03 = cluster_data_values[11];
	cluster_iso_its_04 = cluster_data_values[12];
	cluster_frixione_tpc_04_02 = cluster_data_values[13];
	cluster_frixione_its_04_02 = cluster_data_values[14];
	cluster_s_nphoton0 = cluster_data_values[15];
	cluster_mc_truth_index = cluster_data_values[16];//last dimension is length 32. Placeholder
	cluster_ncell = cluster_data_values[17];
	cluster_cell_id_max = cluster_data_values[18];
	cluster_cell_e = cluster_data_values[19];//last dimension is length 17664. Placeholder
	cluster_distance_to_bad_channel = cluster_data_values[20];
	cluster_nlocal_maxima = cluster_data_values[21];
	cluster_tof = cluster_data_values[22];
	cluster_iso_its_02_ue = cluster_data_values[23];
	cluster_iso_its_03_ue = cluster_data_values[24];
	cluster_iso_its_04_ue = cluster_data_values[25];
	cluster_iso_tpc_02_ue = cluster_data_values[26];
	cluster_iso_tpc_03_ue = cluster_data_values[27];
	cluster_iso_tpc_04_ue = cluster_data_values[28];
	cluster_lambda1 = cluster_data_values[29];
	cluster_s_nphoton = cluster_data_values[30];
	cluster_5x5all = cluster_data_values[31];
}

void setTrackVariables(float track_data_values[])
{
	track_pt = track_data_values[1];
	track_eta = track_data_values[2];
	track_phi = track_data_values[3];
}

void setJetVariables(float jet_data_values[])
{
	jet_pt_raw = jet_data_values[0];
	jet_eta = jet_data_values[1];
	jet_phi = jet_data_values[2];
}

bool rejectCluster()
{
	if (not(cluster_pt > cluster_pt_min and cluster_pt < cluster_pt_max)) return true;
	if (not(abs(cluster_eta) < cluster_eta_max)) return true;
	if (not(cluster_ncell >= cluster_ncell_min)) return true;
	if (not(cluster_e_cross / cluster_e_max > cluster_ecross_emax_min)) return true;
	if (not(cluster_distance_to_bad_channel >= cluster_dbc_min)) return true;
	if (not(cluster_nlocal_maxima < cluster_nlm_max)) return true;
	if (not(abs(cluster_tof) < cluster_tof_max)) return true;
	return false;
}

float getTriggeredIsolation()
{
	float isolation;
	if (isovar == "cluster_iso_tpc_04") {
		isolation = cluster_iso_tpc_04;
	} else if (isovar == "cluster_iso_its_04") {
		isolation = cluster_iso_its_04;
	} else if (isovar == "cluster_iso_its_04_sub") {
		isolation = cluster_iso_its_04 + cluster_iso_its_04_ue - ue_estimate_its_const * M_PI * 0.4 * 0.4;
	} else if (isovar == "cluster_iso_tpc_02_sub") {
		isolation = cluster_iso_tpc_02 + cluster_iso_tpc_02_ue - ue_estimate_tpc_const * M_PI * 0.2 * 0.2;
	} else if (isovar == "cluster_iso_tpc_04_sub") {
		isolation = cluster_iso_tpc_04 + cluster_iso_tpc_04_ue - ue_estimate_tpc_const * M_PI * 0.4 * 0.4;
	} else if (isovar == "cluster_frixione_tpc_04_02") {
		isolation = cluster_frixione_tpc_04_02;
	} else if (isovar == "cluster_frixione_its_04_02") {
		isolation = cluster_frixione_its_04_02;
	} else {
		std::cout << "ERROR: Isolation variable " << isovar << " not recognized. Aborting" << std::endl;
		exit(EXIT_FAILURE);
	}
	return isolation;
}

float getShower()
{
	float shower;
	if (shower_shape == "cluster_Lambda") {
		shower = cluster_lambda_square;
	} else if (shower_shape == "cluster_NN1") {
		shower = cluster_s_nphoton;
	} else if (shower_shape == "cluster_emax_over_e") {
		shower = cluster_e_max / cluster_e;
	} else if (shower_shape == "cluster_5x5all") {
		shower = cluster_5x5all;
	} else {
		std::cout << "ERROR: Shower shape variable " << shower_shape << " not recognized. Aborting" << std::endl;
		exit(EXIT_FAILURE);
	}
	return shower;
}
