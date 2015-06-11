/*
 * patientUnsupervisedNormalization.hpp
 *
 *  Created on: Apr 30, 2015
 *      Author: nicolas
 */

#ifndef SRC_UNSUPERVISEDNORMALIZATION_HPP_
#define SRC_UNSUPERVISEDNORMALIZATION_HPP_

#include "typedefs.hpp"

enum UnsupervisedNormalizationMethod {
	KMEANS, BINARY_QUANTILE, NO_NORMALIZATION
};

struct UnsupervisedNormalizationParameters {
	//For KMEANS
	unsigned int K;
	int Nmax;

	//For BINARY_QUANTILE
	double cutPercentage;

	void setKMeansParameters(unsigned int _K, int _Nmax) {
		K = _K;
		Nmax = _Nmax;
	}

	void setBinaryQuantileParameters(double _cutPercentage) {
		cutPercentage = _cutPercentage;
	}
};

void unsupervisedNormalization(Data &data,
		const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters);

void printMaxExpressedGenes(const Data &data, unsigned int maxNumberGenes,
		const std::string &filename);

#endif /* SRC_UNSUPERVISEDNORMALIZATION_HPP_ */
