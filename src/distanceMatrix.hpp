/*
 * correlationMatrix.hpp
 *
 *  Created on: May 5, 2015
 *      Author: nicolas
 */

#ifndef SRC_DISTANCEMATRIX_HPP_
#define SRC_DISTANCEMATRIX_HPP_

#include <vector>

#include "typedefs.hpp"

enum DistanceMetric {
	PEARSON_CORRELATION,
	SPEARMAN_CORRELATION,
	EUCLIDEAN_DISTANCE,
	MANHATTAN_DISTANCE,
	COSINE_SIMILARITY
};

enum MatrixType {
	DISTANCE, SIMILARITY
};

std::string distanceMetricName(const DistanceMetric &distanceMetric);

std::vector<double> computeDistanceMatrix(
		const std::vector<std::vector<double>> &data, DistanceMetric method);

MatrixType getMatrixType(const DistanceMetric &method);

void exportDistanceMatrix(const std::vector<double> &distanceMatrix,
		const std::vector<SampleIdentifier> &sampleIdentifiers,
		const std::string &filemaneMatrix,
		const std::string &filenamePatientsIds,
		const std::string &filenameHeatMapLabels);

void exportClassStats(const std::vector<double> &distanceMatrix,
		const CancerPatientIDList &cancerPatientIDList,
		const std::vector<SampleIdentifier> &sampleIdentifiers,
		const std::string &filemaneCorrelationMeans);

#endif /* SRC_DISTANCEMATRIX_HPP_ */
