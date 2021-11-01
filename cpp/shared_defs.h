#include <vector>
#include "emcal.h"
#include "yaml-cpp/yaml.h"

/*--------------------------------------------------------------
Isolation calculations with and without centrality-dependence
--------------------------------------------------------------*/
bool GetIsIsolated(float isolation, float centrality, YAML::Node isoConfig) {
	// centrality-dependent isolation
	if (isoConfig["centralityvalues"]) {
		YAML::Node centvalues = isoConfig["centralityvalues"];
		for (YAML::const_iterator it = centvalues.begin(); it != centvalues.end(); it++) {
			const YAML::Node& centvalue = *it;
			std::pair<float, float> centrange = centvalue["centralityrange"].as<std::pair<float, float>>();
			if (centrality >= centrange.first && centrality < centrange.second) {
				return isolation < centvalue["isocut"].as<float>();
			}
		}
		// out of centrality range
		std::cout << "centrality " << centrality << " out of range in isolation calculation" << std::endl;
		return false;
	// flat isolation
	} else {
		float isocut = isoConfig["isocut"].as<float>();
		return isolation < isocut;
	}
}

/*--------------------------------------------------------------
Purity function definitions with and without centrality-dependence
--------------------------------------------------------------*/
float GetPurityBinned(float pt, std::vector<float> binEdges, std::vector<float> values) {
	for (int binNumber = 0; binNumber < values.size(); ++binNumber) {
		if (pt >= binEdges[binNumber] && pt < binEdges[binNumber + 1]) {
			return values[binNumber];
		}
	}
	std::cout << "pT " << pt << " out of range in purity calculation" << std::endl;
	return 0;
}

float GetPurityBinnedCentrality(float centrality, float pt, std::vector<float> binEdges, std::map<std::pair<float, float>, std::vector<float>> allValues) {
	for (auto centValues : allValues) {
		std::pair<float, float> centrange = centValues.first;
		if (centrality >= centrange.first && centrality < centrange.second) {
			std::vector<float> values = centValues.second;
			return GetPurityBinned(pt, binEdges, values);
		}
	}
	std::cout << "centrality " << centrality << " out of range in purity calculation" << std::endl;
	return 0;
}

float GetPurityErf(float pt, std::vector<float> params) {
	float a = params[0];
	float b = params[1];
	float c = params[2];

	return a * TMath::Erf((pt - b) / c);
}

float GetPurityErfCentrality(float centrality, float pt, std::map<std::pair<float, float>, std::vector<float>> allParams) {
	for (auto centParams : allParams) {
		std::pair<float, float> centrange = centParams.first;
		if (centrality >= centrange.first && centrality < centrange.second) {
			std::vector<float> params = centParams.second;
			return GetPurityErf(pt, params);
		}
	}
	std::cout << "centrality " << centrality << " out of range in purity calculation" << std::endl;
	return 0;
}

// figure out which function to call based on the config
float getPurity(float pt, float centrality, YAML::Node purityConfig)
{
	YAML::Node params = purityConfig["functionparams"];
	if (purityConfig["function"].as<std::string>() == "binned") {
		std::vector<float> binEdges = params["binEdges"].as<std::vector<float>>();
		if (params["centralityvalues"]) {
			std::map<std::pair<float, float>, std::vector<float>> allValues;
			YAML::Node centvalues = params["centralityvalues"];

			for (YAML::const_iterator it = centvalues.begin(); it != centvalues.end(); it++) {
				const YAML::Node& centvalue = *it;
				std::pair<float, float> centrange = centvalue["centralityrange"].as<std::pair<float, float>>();
				std::vector<float> values = centvalue["values"].as<std::vector<float>>();
				allValues[centrange] = values;
			}

			return GetPurityBinnedCentrality(centrality, pt, binEdges, allValues);
		} else {
			std::vector<float> values = params["values"].as<std::vector<float>>();
			return GetPurityBinned(pt, binEdges, values);
		}
	}
	else if (purityConfig["function"].as<std::string>() == "erf") {
		if (params["centralityvalues"]) {
			std::map<std::pair<float, float>, std::vector<float>> allValues;
			YAML::Node centvalues = params["centralityvalues"];

			for (YAML::const_iterator it = centvalues.begin(); it != centvalues.end(); it++) {
				const YAML::Node& centvalue = *it;
				std::pair<float, float> centrange = centvalue["centralityrange"].as<std::pair<float, float>>();
				std::vector<float> erfparams;
				erfparams.push_back(params["a"].as<float>());
				erfparams.push_back(params["b"].as<float>());
				erfparams.push_back(params["c"].as<float>());
				allValues[centrange] = erfparams;
			}

			return GetPurityErfCentrality(centrality, pt, allValues);

		} else {
			std::vector<float> erfparams;
			erfparams.push_back(params["a"].as<float>());
			erfparams.push_back(params["b"].as<float>());
			erfparams.push_back(params["c"].as<float>());

			return GetPurityErf(pt, erfparams);
		}
	}
}

// calculate 5x5all shower shape
float get5x5all(const unsigned int cellMaxId, float cluster_e, float cell_e[17664])
{
	unsigned int cells5x5[25];
	cell_5_5(cells5x5, cellMaxId);

	float wtot = 0.0;
    float w = 0.0;
	float x = 0.0;
	float z = 0.0;
	float dxx = 0.0;
	float dzz = 0.0;
	float dxz = 0.0;

	unsigned int sm_max;
	unsigned int nphi_max;
	to_sm_nphi(sm_max, nphi_max, cellMaxId);

	for (int i = 0; i < 25; i++) {
		unsigned int cellId = cells5x5[i];

        if (cellId < 0 || cellId > 17763) continue;
		if (cell_e[cellId] < 0.1 || TMath::IsNaN(cell_e[cellId])) continue;

		unsigned int sm;
		unsigned int ieta;
		unsigned int iphi;
		to_sm_ieta_iphi(sm, ieta, iphi, cellId);

		if (sm % 2) {
			ieta = ieta + 48;
		}

		float w = TMath::Max(0.0, 4.5 + TMath::Log(cell_e[cellId] / cluster_e));
		dxx = dxx + w * ieta * ieta;
		x = x + w * ieta;
		dzz = dzz + w * iphi * iphi;
		z = z + w * iphi;
		dxz = dxz + w * ieta * iphi;
		wtot = wtot + w;
	}

	if (wtot > 0) {
		dxx /= wtot ;
		x   /= wtot ;
		dxx -= x * x ;
		dzz /= wtot ;
		z   /= wtot ;
		dzz -= z * z ;
		dxz /= wtot ;
		dxz -= x * z ;
	}

	return 0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz );
}
