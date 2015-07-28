#include "TCGADataNormalizer.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <Eigen/Dense>
#include <ClusterXX/clustering/kmeans_clusterer.hpp>
#include <ClusterXX/utils/utils.hpp>
#include "../config.hpp"
#include "../utilities.hpp"

void KMeansNormalizer::normalize(std::vector<double> *v) {
	Eigen::Map<Eigen::MatrixXd> mapToData((*v).data(), 1, (*v).size());
	std::shared_ptr<ClusterXX::ClustererParameters> kMeansParams =
			std::make_shared<ClusterXX::KMeansParameters>(K, maxIterations);
	ClusterXX::KMeans_Clusterer clusterer(mapToData, kMeansParams);
	clusterer.compute();
	std::vector<int> clusters = clusterer.getClusters();
	for (unsigned int i = 0; i < v->size(); ++i) {
		(*v)[i] = (double) clusters[i];
	}
}

void BinaryQuantileNormalizer::normalize(std::vector<double> *v) {
	unsigned int numberOfGenes = v->size();
	std::vector<size_t> rank = ClusterXX::Utilities::get_rank_increasing(*v);
	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		double p = (double) rank[i] / (double) numberOfGenes;
		if (p >= cutPercentage) {
			(*v)[i] = 1.0;
		} else {
			(*v)[i] = 0.0;
		}
	}
}

TCGADataNormalizer::TCGADataNormalizer(TCGAData *_ptrToData,
		const std::shared_ptr<Normalizer> &_ptrToNormalizer, bool _verbose) :
		ptrToData(_ptrToData), ptrToNormalizer(_ptrToNormalizer), verbose(
				_verbose) {
}

void TCGADataNormalizer::normalizeIndividualSample(unsigned int sampleId) {
	unsigned int numberOfGenes = ptrToData->getNumberOfGenes();
	std::vector<double> dataToNormalize(numberOfGenes);
	//Copy the data in a vector
	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		dataToNormalize[i] = ptrToData->getDataHandler()[i][sampleId];
	}
	//Normalize
	ptrToNormalizer->normalize(&dataToNormalize);
	//Copy the normalized data back into the structure
	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		ptrToData->getDataHandler()[i][sampleId] = dataToNormalize[i];
	}
}

void TCGADataNormalizer::normalize() {
	if (verbose) {
		std::cout << "Normalizing data... " << std::flush;
	}
	for(unsigned int i=0; i<ptrToData->getNumberOfSamples(); ++i){
		normalizeIndividualSample(i);
	}
	if (verbose) {
		std::cout << "Done." << std::endl;
	}
}

void TCGADataNormalizer::exportToFile(double positiveValue,
		double negativeValue) {
	ptrToData->buildDataMatrix();
	std::ofstream outputStreamSamples(HEINZ_SAMPLES_LIST);
	const auto &dataMatrix = ptrToData->getDataMatrixHandler();
	unsigned int N = dataMatrix.cols();
	for (unsigned int i = 0; i < N; ++i) {

		std::string outputFilename =
				ptrToData->getPatientsHandler()[i].toString();
		outputStreamSamples << outputFilename << std::endl;

		std::ofstream outputStream(
				HEINZ_INPUT_DIRECTORY
						+ removeTrailingZeros(
								std::to_string(std::fabs(negativeValue))) + '_'
						+ outputFilename + ".txt");

		for (unsigned int j = 0; j < dataMatrix.rows(); ++j) {
			outputStream << ptrToData->getGeneListHandler()[j].first << " ";
			if (dataMatrix(j, i) > 0.5) {
				outputStream << positiveValue << std::endl;
			} else {
				outputStream << negativeValue << std::endl;
			}
		}

		ClusterXX::Utilities::printAdvancement(i, N);
	}
}
