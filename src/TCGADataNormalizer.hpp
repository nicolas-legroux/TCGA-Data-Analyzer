/*
 * patientUnsupervisedNormalization.hpp
 *
 *  Created on: Apr 30, 2015
 *      Author: nicolas
 */

#ifndef SRC_TCGADATANORMALIZER_HPP_
#define SRC_TCGADATANORMALIZER_HPP_

#include "TCGAData.hpp"
#include <memory>
#include <fstream>

enum UnsupervisedNormalizationMethod {
	KMEANS_NORMALIZATION, BINARY_QUANTILE_NORMALIZATION, NO_NORMALIZATION
};

class Normalizer {
public:
	virtual void normalize(std::vector<double> *v) = 0;
	virtual ~Normalizer() = default;
};

class NoOperationNormalizer: Normalizer {
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
	void printMostExpressedGenesByClass(unsigned int maxNumberGenes,
			const std::string &filename);
private:
	TCGAData *ptrToData;
	std::shared_ptr<Normalizer> ptrToNormalizer;
	bool verbose;
	void normalizeIndividualSample(const std::string &cancer, bool isTumor,
			unsigned int patientID);
	void printMostExpressedGenesByClassUtility(std::ofstream &outputStream, unsigned int maxNumberGenes,
			std::string &cancer, bool isTumor);
};

#endif /* SRC_TCGADATANORMALIZER_HPP_ */
