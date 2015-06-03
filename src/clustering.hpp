#ifndef SRC_CLUSTERING_HPP_
#define SRC_CLUSTERING_HPP_

#include <vector>
#include "typedefs.hpp"
#include "hierarchical_clustering.hpp"
#include "distanceMatrix.hpp"

// Function to get real labels
std::vector<int> getRealClusters(
		std::vector<SampleIdentifier> &sampleIdentifiers);

//Functions to compute clusters
std::vector<int> cluster_KMeans(const std::vector<std::vector<double>> &data,
		int K, int Nmax);
std::vector<int> cluster_Hierarchical(const std::vector<double> &matrix,
		const DistanceMetric &distanceMetric,
		const LinkageMethod &linkageMethod, int K);

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
