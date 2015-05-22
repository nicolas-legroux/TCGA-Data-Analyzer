/*
 * clustering.hpp
 *
 *  Created on: May 22, 2015
 *      Author: nicolas
 */

#ifndef SRC_CLUSTERING_HPP_
#define SRC_CLUSTERING_HPP_

#include <vector>
#include "typedefs.hpp"

std::vector<int> getRealClusters(std::vector<SampleIdentifier> &sampleIdentifiers);
std::vector<int> clusterWithKMeans(const std::vector<std::vector<double>> &data, int K, int Nmax);



#endif /* SRC_CLUSTERING_HPP_ */
