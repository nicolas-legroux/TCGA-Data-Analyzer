/*
 * parameters.hpp
 *
 *  Created on: Jun 24, 2015
 *      Author: nicolas
 */

#ifndef PARAMETERS_HPP_
#define PARAMETERS_HPP_

#include <vector>
#include <ClusterXX/metrics/metrics.hpp>
#include <ClusterXX/clustering/clusterer_parameters.hpp>

const std::string GRAPH_NODE_FILE = "biogrid-nodes.txt";

const bool VERBOSE = true;

const unsigned int MAX_CONTROL_SAMPLES = 20;
const unsigned int MAX_TUMOR_SAMPLES = 20;
const std::vector<std::string> CANCERS = { "LUSC", "LUAD" };

const int K_MEANS_NORMALIZATION_PARAM = 2;
const double BINARY_QUANTILE_NORMALIZATION_PARAM = 0.995;

const int K_MEANS_MAX_ITERATIONS = 1000;

const ClusterXX::MetricName::MetricName DEFAULT_METRIC =
		ClusterXX::MetricName::JACCARD_SIMILARITY;

const ClusterXX::HierarchicalParameters::LinkageMethod DEFAULT_LINKAGE_METHOD =
		ClusterXX::HierarchicalParameters::COMPLETE;

const ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName SPECTRAL_GRAPH_K_NEAREST_NEIGHBORS =
		ClusterXX::SpectralParameters::GraphTransformationMethod::K_NEAREST_NEIGHBORS;

const ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName SPECTRAL_GRAPH_GAUSSIAN_MIXTURE =
		ClusterXX::SpectralParameters::GraphTransformationMethod::GAUSSIAN_MIXTURE;

const ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName SPECTRAL_GRAPH_NO_TRASNFORMATION =
		ClusterXX::SpectralParameters::GraphTransformationMethod::NO_TRANSFORMATION;

const int SPECTRAL_K_NEAREST_NEIGHBORS = 10;
const double SPECTRAL_GAUSSIAN_MIXTURE_STDDEV = 150.0;

const std::pair<
		ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName,
		double> DEFAULT_GRAPH_TRANSFORMATION = {
				SPECTRAL_GRAPH_K_NEAREST_NEIGHBORS, SPECTRAL_K_NEAREST_NEIGHBORS };


#endif /* PARAMETERS_HPP_ */
