/*
 * correlationMatrix.hpp
 *
 *  Created on: May 5, 2015
 *      Author: nicolas
 */

#ifndef SRC_CORRELATIONMATRIX_HPP_
#define SRC_CORRELATIONMATRIX_HPP_

#include <vector>

#include "typedefs.hpp"

std::vector<double> pearson(std::vector<std::vector<double>> &data);
std::vector<double> spearman(std::vector<std::vector<double>> &data);

std::vector<double> euclidean(std::vector<std::vector<double>> &data);
std::vector<double> manhattan(std::vector<std::vector<double>> &data);

void exportCorrelationMatrix(const std::vector<double> &correlationMatrix,
		const std::vector<SampleIdentifier> &sampleIdentifiers,
		const std::string &filemaneMatrix,
		const std::string &filenamePatientsIds,
		const std::string &filenameHeatMapLabels);

void exportClassStats(const std::vector<double> &correlationMatrix,
		const CancerPatientIDList &cancerPatientIDList,
		const std::vector<SampleIdentifier> &sampleIdentifiers,
		const std::string &filemaneCorrelationMeans);

#endif /* SRC_CORRELATIONMATRIX_HPP_ */
