/*
 * correlationMatrix.hpp
 *
 *  Created on: May 5, 2015
 *      Author: nicolas
 */

#ifndef SRC_TCGADATADISTANCEMATRIXANALYZER_HPP_
#define SRC_TCGADATADISTANCEMATRIXANALYZER_HPP_

#include <memory>
#include <array>
#include <ClusterXX/metrics/metrics.hpp>
#include "TCGAData.hpp"

class TCGADataDistanceMatrixAnalyser {
public:
	TCGADataDistanceMatrixAnalyser(TCGAData *_ptrToData,
			const std::shared_ptr<ClusterXX::Metric> &_metric,
			bool _verbose = false);
	void computeDistanceMatrix();
	void exportDistanceMatrix();
	void exportClassStats();
	void exportHeatMap(bool withClassDivision = true,
			std::array<unsigned char, 3> separatorColor = std::array<
					unsigned char, 3> { static_cast<unsigned char>(255),
					static_cast<unsigned char>(155), 0 });
	Eigen::MatrixXd &getDistanceMatrixHandler() {
		return distanceMatrix;
	}
private:
	TCGAData *ptrToData;
	std::shared_ptr<ClusterXX::Metric> metric;
	Eigen::MatrixXd distanceMatrix;
	bool verbose;
	bool matrixIsComputed;

	std::vector<unsigned int> buildClassDivisionForHeatmap();
};

#endif /* SRC_TCGADATADISTANCEMATRIXANALYZER_HPP_ */
