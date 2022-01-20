#include "yaml-cpp/yaml.h"

std::vector<YAML::Node> allconfigs;

double srmin = 0;
double srmax = 0;
double brmin = 0;
double brmax = 0;
double cluster_pt_min = 0;
double cluster_pt_max = 0;
double cluster_eta_max = 0;
double cluster_ncell_min = 0;
double cluster_dbc_min = 0;
int cluster_nlm_max = 3;
double cluster_ecross_emax_min = 0;
double cluster_tof_max = 20;

double jet_pt_min = 0;
double jet_pt_max = 500;
double jet_eta_max = 10;

double deltaphi_min = 0;
double deltaphi_max = M_PI;
double jetpt_min = 5;
double jetpt_max = 50;
double ptratio_min = 0;
double ptratio_max = 2;

bool do_pile = false;

std::string jettype = "ak04tpc";
std::string isovar = "cluster_iso_tpc_02_sub";
std::string shower_shape = "DNN";
std::string purity_deviation = "None";

YAML::Node purityconfig;
YAML::Node isoconfig;
YAML::Node centralityranges;

int ncentralityranges;
bool keepFakes;
bool keepMisses;

void parseConfig()
{
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
		}

		if (config["clustercuts"]["all"]["cluster_pt"]) {
			cluster_pt_min = config["clustercuts"]["all"]["cluster_pt"]["min"].as<double>();
			cluster_pt_max = config["clustercuts"]["all"]["cluster_pt"]["max"].as<double>();
		}

		if (config["clustercuts"]["all"]["cluster_eta"]) {
			cluster_eta_max = config["clustercuts"]["all"]["cluster_eta"]["absmax"].as<double>();
		}

		if (config["clustercuts"]["all"]["cluster_ncell"]) {
			cluster_ncell_min = config["clustercuts"]["all"]["cluster_ncell"]["incmin"].as<double>();
		}

		if (config["clustercuts"]["all"]["cluster_distance_to_bad_channel"]) {
			cluster_dbc_min = config["clustercuts"]["all"]["cluster_distance_to_bad_channel"]["incmin"].as<double>();
		}

		if (config["clustercuts"]["all"]["cluster_nlocal_maxima"]) {
			cluster_nlm_max = config["clustercuts"]["all"]["cluster_nlocal_maxima"]["max"].as<int>();
		}

		if (config["clustercuts"]["all"]["cluster_ecross_emax"]) {
			cluster_ecross_emax_min = config["clustercuts"]["all"]["cluster_ecross_emax"]["min"].as<double>();
		}

		if (config["clustercuts"]["data"]["cluster_tof"]) {
			cluster_tof_max = config["clustercuts"]["data"]["cluster_tof"]["absmax"].as<double>();
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
		}

		if (config["Purity_Dev"]) {
			purity_deviation = config["Purity_Dev"].as<std::string>();
		}

		if (config["purity"]) {
			purityconfig = config["purity"];
		}

		if (config["responsematrix"]) {
		  keepMisses = config["responsematrix"]["keepmisses"].as<bool>();
		  keepFakes = config["responsematrix"]["keepfakes"].as<bool>();
		}

		if (config["centralityranges"]) {
		  centralityranges = config["centralityranges"];
		  ncentralityranges = centralityranges.size();
		}

		if (config["observables"]) {
			YAML::Node observables = config["observables"];
			for (YAML::const_iterator it = observables.begin(); it != observables.end(); it++) {
				YAML::Node observableInfo = *it;
				if (observableInfo["name"].as<std::string>() == "deltaphi") {
					deltaphi_min = observableInfo["min"].as<double>();
					deltaphi_max = observableInfo["max"].as<double>();
				}
				if (observableInfo["name"].as<std::string>() == "jetpt") {
					jetpt_min = observableInfo["min"].as<double>();
					jetpt_max = observableInfo["max"].as<double>();
				}
				if (observableInfo["name"].as<std::string>() == "ptratio") {
					ptratio_min = observableInfo["min"].as<double>();
					ptratio_max = observableInfo["max"].as<double>();
				}
			}
		}
	}
}
