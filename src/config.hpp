/*
 * config.hpp
 *
 *  Created on: Jun 17, 2015
 *      Author: nicolas
 */

#ifndef SRC_CONFIG_HPP_
#define SRC_CONFIG_HPP_

#include <vector>

#include <ClusterXX/clustering/clusterer_parameters.hpp>
#include <ClusterXX/metrics/metrics.hpp>

const std::string DATA_ROOT_DIRECTORY = "data/";
const std::string ROOT_EXPORT_DIRECTORY = "export/";
const std::string SAMPLE_DATA_FILE = DATA_ROOT_DIRECTORY
		+ "BRCA-normalized/TCGA-3C-AAAU-01.genes.normalized.results";

const unsigned int MAX_CONTROL_SAMPLES = 100;
const unsigned int MAX_TUMOR_SAMPLES = 200;
const std::vector<std::string> CANCERS = { "BRCA", "LUAD",  "KIRC", "THCA", "OV"};

const int K_MEANS_NORMALIZATION_PARAM = 2;
const double BINARY_QUANTILE_NORMALIZATION_PARAM = 0.45;

const bool VERBOSE = true;

const int K_MEANS_MAX_ITERATIONS = 1000;

const ClusterXX::MetricName::MetricName DEFAULT_METRIC =
		ClusterXX::MetricName::MANHATTAN_DISTANCE;

const ClusterXX::HierarchicalParameters::LinkageMethod DEFAULT_LINKAGE_METHOD =
		ClusterXX::HierarchicalParameters::COMPLETE;

const ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName SPECTRAL_GRAPH_K_NEAREST_NEIGHBORS =
		ClusterXX::SpectralParameters::GraphTransformationMethod::K_NEAREST_NEIGHBORS;

const ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName SPECTRAL_GRAPH_GAUSSIAN_MIXTURE =
		ClusterXX::SpectralParameters::GraphTransformationMethod::GAUSSIAN_MIXTURE;

const ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName SPECTRAL_GRAPH_NO_TRASNFORMATION =
		ClusterXX::SpectralParameters::GraphTransformationMethod::NO_TRANSFORMATION;

const int SPECTRAL_K_NEAREST_NEIGHBORS = 6;
const double SPECTRAL_GAUSSIAN_MIXTURE_STDDEV = 1.0;

const std::pair<
		ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName,
		double> DEFAULT_GRAPH_TRANSFORMATION = {
		SPECTRAL_GRAPH_K_NEAREST_NEIGHBORS, SPECTRAL_K_NEAREST_NEIGHBORS };

#endif /* SRC_CONFIG_HPP_ */
