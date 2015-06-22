// Copyright Nicolas Legroux 2015

#include <memory>
#include "config.hpp"
#include "TCGAData.hpp"
#include "TCGADataLoader.hpp"
#include "TCGADataClusterer.hpp"
#include "TCGADataNormalizer.hpp"
#include "TCGADataDistanceMatrixAnalyzer.hpp"

int main() {
	/* Read Data */
	TCGAData data;
	TCGADataLoader loader(&data, CANCERS, MAX_CONTROL_SAMPLES,
			MAX_TUMOR_SAMPLES, VERBOSE);
	loader.loadData();

	/* Normalize */
	std::shared_ptr<Normalizer> normalizer = std::make_shared<
			BinaryQuantileNormalizer>(BINARY_QUANTILE_NORMALIZATION_PARAM);
	TCGADataNormalizer tcgaNormalizer(&data, normalizer, VERBOSE);
	tcgaNormalizer.normalize();

	/* Output distance matrix */
	auto metric = ClusterXX::buildMetric(DEFAULT_METRIC);
	TCGADataDistanceMatrixAnalyser distanceMetricAnalyzer(&data, metric,
			VERBOSE);
	distanceMetricAnalyzer.computeDistanceMatrix();
	distanceMetricAnalyzer.exportClassStats();
	distanceMetricAnalyzer.exportHeatMap();

	TCGADataKMeansClusterer kMeansClusterer(&data);
	kMeansClusterer.computeClustering();
	kMeansClusterer.printClusteringInfo();

	TCGADataHierarchicalClusterer hierarchicalClusterer(&data,
			distanceMetricAnalyzer.getDistanceMatrixHandler(), metric);
	hierarchicalClusterer.computeClustering();
	hierarchicalClusterer.printClusteringInfo();

	TCGADataSpectralClusterer spectralClusterer(&data,
			distanceMetricAnalyzer.getDistanceMatrixHandler(), metric);
	spectralClusterer.computeClustering();
	spectralClusterer.printClusteringInfo();
}
