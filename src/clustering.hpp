#ifndef SRC_CLUSTERING_HPP_
#define SRC_CLUSTERING_HPP_

#include <vector>
#include "typedefs.hpp"

std::vector<int> getRealClusters(std::vector<SampleIdentifier> &sampleIdentifiers);
std::vector<int> clusterWithKMeans(const std::vector<std::vector<double>> &data, int K, int Nmax);

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
