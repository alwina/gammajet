#include <vector>
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
