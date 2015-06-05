#ifndef SRC_CLUSTERING_HPP_
#define SRC_CLUSTERING_HPP_

#include <vector>
#include <map>
#include "typedefs.hpp"
#include "hierarchical_clustering.hpp"
#include "distanceMatrix.hpp"

// Function to get real labels
std::vector<int> getRealClusters(
		std::vector<SampleIdentifier> &sampleIdentifiers);

std::map<int, std::string> getRealLabelsMap(
		std::vector<SampleIdentifier> &sampleIdentifiers);

//Functions to compute clusters
std::vector<int> cluster_KMeans(const std::vector<std::vector<double>> &data,
		int K, int Nmax, bool verbose = false);
std::vector<int> cluster_Hierarchical(const std::vector<double> &matrix,
		const DistanceMetric &distanceMetric,
		const LinkageMethod &linkageMethod, int K, bool verbose = false);
std::vector<int> cluster_Spectral(const std::vector<double> &matrix, unsigned int K);

void printClustering(const std::map<int, std::string> &labelsMap,
		const std::vector<int> &realClusters,
		const std::vector<int> &computedClusters);
/*
 *
 * CLUSTERING EVALUATION
 *
 */

double randIndex(const std::vector<int> &clustering1,
		const std::vector<int> &clustering2);

double adjustedRandIndex(const std::vector<int> &clustering1,
		const std::vector<int> &clustering2);

#endif /* SRC_CLUSTERING_HPP_ */
