// Copyright Nicolas Legroux 2015

#include <memory>
#include <regex>
#include "parameters.hpp"
#include "TCGAData.hpp"
#include "TCGADataLoader.hpp"
#include "TCGADataClusterer.hpp"
#include "TCGADataNormalizer.hpp"
#include "TCGADataDistanceMatrixAnalyzer.hpp"

int main() {

//	std::string to_match = "595 [label=\"XBP1\\n1\\n595\"]";
//	std::string reg = "([0-9]+) \\[label=\"([A-Z0-1]+)";
//	std::smatch match;
//	std::regex e(reg);
//	std::regex_search(to_match, match, e);
//	std::cout << match[1] << " " << match[2];

	/* Read Data */
	TCGAData data;
	TCGADataLoader loader(&data, CANCERS, MAX_CONTROL_SAMPLES,
			MAX_TUMOR_SAMPLES, VERBOSE);
	loader.loadData();

	data.keepOnlyGenesInGraph(GRAPH_NODE_FILE);

	/* Normalize */
	std::shared_ptr<Normalizer> normalizer = std::make_shared<
			BinaryQuantileNormalizer>(BINARY_QUANTILE_NORMALIZATION_PARAM);
	TCGADataNormalizer tcgaNormalizer(&data, normalizer, VERBOSE);
	tcgaNormalizer.normalize();

	for(double d = 0.5; d<=8.0; d+= 0.25){
		tcgaNormalizer.exportToFile(1, -1.0*d);
	}
//
//
//	/* Output distance matrix */
//	auto metric = ClusterXX::buildMetric(DEFAULT_METRIC);
//	TCGADataDistanceMatrixAnalyser distanceMetricAnalyzer(&data, metric,
//			VERBOSE);
//	distanceMetricAnalyzer.computeDistanceMatrix();
//	distanceMetricAnalyzer.exportClassStats();
//	distanceMetricAnalyzer.exportHeatMap();
////
////	TCGADataKMeansClusterer kMeansClusterer(&data, 0, K_MEANS_MAX_ITERATIONS, VERBOSE);
////	kMeansClusterer.computeClustering();
////	kMeansClusterer.printClusteringInfo();
////
////	TCGADataHierarchicalClusterer hierarchicalClusterer(&data,
////			distanceMetricAnalyzer.getDistanceMatrixHandler(), metric, 0, DEFAULT_LINKAGE_METHOD, false);
////	hierarchicalClusterer.computeClustering();
////	hierarchicalClusterer.printClusteringInfo();
////
////	TCGADataSpectralClusterer spectralClusterer(&data,
////			distanceMetricAnalyzer.getDistanceMatrixHandler(), metric, 0, DEFAULT_GRAPH_TRANSFORMATION, VERBOSE);
////	spectralClusterer.computeClustering();
////	spectralClusterer.printClusteringInfo();

}
