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
	KMEANS, BINARY_ITERATED_KMEANS, BINARY_QUANTILE, NO_NORMALIZATION
};

struct UnsupervisedNormalizationParameters {
	//For KMEANS
	unsigned int K;
	int Nmax;

	//For BINARY_ITERATED_KMEANS
	int Niteration;

	//For BINARY_QUANTILE
	double cutPercentage;

	void setKMeansParameters(unsigned int _K, int _Nmax) {
		K = _K;
		Nmax = _Nmax;
	}

	void setBinaryIteratedKMeansParameters(int _Niteration) {
		Niteration = _Niteration;
	}

	void setBinaryQuantileParameters(double _cutPercentage) {
		cutPercentage = _cutPercentage;
	}
};

void unsupervisedNormalization(Data &data,
		const UnsupervisedNormalizationMethod &method,
		UnsupervisedNormalizationParameters &parameters);

void printMaxExpressedGenes(const Data &data, unsigned int maxNumberGenes,
		const std::string &filename);

#endif /* SRC_UNSUPERVISEDNORMALIZATION_HPP_ */
