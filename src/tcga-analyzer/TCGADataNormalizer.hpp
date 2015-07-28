/*
 * patientUnsupervisedNormalization.hpp
 *
 *  Created on: Apr 30, 2015
 *      Author: nicolas
 */

#ifndef SRC_TCGADATANORMALIZER_HPP_
#define SRC_TCGADATANORMALIZER_HPP_

#include <memory>
#include <fstream>
#include "../tcga-analyzer/TCGAData.hpp"
#include "../config.hpp"

class Normalizer {
public:
	virtual void normalize(std::vector<double> *v) = 0;
	virtual ~Normalizer() = default;
};

class NoOperationNormalizer: public Normalizer {
public:
	NoOperationNormalizer() = default;
	void normalize(std::vector<double> *v) {
	}
};

class KMeansNormalizer: public Normalizer {
public:
	KMeansNormalizer(unsigned int _K, unsigned int _maxIterations) :
			K(_K), maxIterations(_maxIterations) {

	}
	void normalize(std::vector<double> *v);
private:
	unsigned int K;
	unsigned int maxIterations;
};

class BinaryQuantileNormalizer: public Normalizer {
public:
	BinaryQuantileNormalizer(double _cutPercentage) :
			cutPercentage(_cutPercentage) {
	}
	void normalize(std::vector<double> *v);
private:
	double cutPercentage;
};

class TCGADataNormalizer {
public:
	TCGADataNormalizer(TCGAData *_ptrToData,
			const std::shared_ptr<Normalizer> &_ptrToNormalizer, bool _verbose =
					true);
	void normalize();
	void exportToFile(double positiveValue, double negativeValuess);
private:
	TCGAData *ptrToData;
	std::shared_ptr<Normalizer> ptrToNormalizer;
	bool verbose;
	void normalizeIndividualSample(unsigned int sampleId);
};

#endif /* SRC_TCGADATANORMALIZER_HPP_ */
