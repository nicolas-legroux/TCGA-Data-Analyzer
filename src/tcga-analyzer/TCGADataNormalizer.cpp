#include "../tcga-analyzer/TCGADataNormalizer.hpp"

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

void TCGADataNormalizer::normalizeIndividualSample(const std::string &cancer,
		bool isTumor, unsigned int patientId) {
	unsigned int numberOfGenes = ptrToData->getNumberOfGenes();
	std::vector<double> dataToNormalize(numberOfGenes);
	//Copy the data in a vector
	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		dataToNormalize[i] =
				ptrToData->getRNASeqDataHandler(isTumor)[cancer][i][patientId];
	}
	//Normalize
	ptrToNormalizer->normalize(&dataToNormalize);
	//Copy the normalized data back into the structure
	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		ptrToData->getRNASeqDataHandler(isTumor)[cancer][i][patientId] =
				dataToNormalize[i];
	}
}

void TCGADataNormalizer::normalize() {
	if (verbose) {
		std::cout << "Normalizing control data..." << std::endl;
	}
	for (auto &kv : ptrToData->getRNASeqDataHandler(false)) {
		std::string cancer = kv.first;
		if (verbose) {
			std::cout << "\t" << cancer << "... " << std::flush;
		}
		for (unsigned int patientId = 0; patientId < kv.second[0].size();
				++patientId) {
			normalizeIndividualSample(cancer, false, patientId);
		}
		if (verbose) {
			std::cout << "Done." << std::endl;
		}
	}
	if (verbose) {
		std::cout << "Normalizing tumor data..." << std::endl;
	}
	for (auto &kv : ptrToData->getRNASeqDataHandler(true)) {
		std::string cancer = kv.first;
		if (verbose) {
			std::cout << "\t" << cancer << "... " << std::flush;
		}
		for (unsigned int patientId = 0; patientId < kv.second[0].size();
				++patientId) {
			normalizeIndividualSample(cancer, true, patientId);
		}
		if (verbose) {
			std::cout << "Done." << std::endl;
		}
	}
}

/*
 *
 * PRINTS MOST EXPRESSED GENES PER CLASS
 *
 */

void TCGADataNormalizer::printMostExpressedGenesByClassUtility(
		std::ofstream &outputStream, unsigned int maxNumberGenes,
		std::string &cancer, bool isTumor) {

	unsigned int numberOfPatients = ptrToData->getRNASeqDataHandler(isTumor).at(
			cancer)[0].size();
	if (numberOfPatients > 0) {
		outputStream << cancer;
		if (isTumor) {
			outputStream << "-Tumor : ";
		} else {
			outputStream << "-Control : ";
		}
		outputStream << numberOfPatients << " patients. " << std::endl;

		unsigned int numberOfGenes = ptrToData->getNumberOfGenes();
		std::vector<double> aggregation(numberOfGenes, 0.0);
		for (unsigned int i = 0; i < numberOfGenes; i++) {
			for (unsigned int j = 0; j < numberOfPatients; ++j) {
				aggregation[i] += ptrToData->getRNASeqDataHandler(isTumor).at(
						cancer)[i][j];
			}
		}

		outputStream << "Most expressed genes in the class : {";
		std::vector<size_t> sortedIndexes =
				ClusterXX::Utilities::sort_indexes_decreasing(aggregation);
		for (unsigned int i = 0; i < maxNumberGenes; ++i) {
			std::string geneSymbol =
					ptrToData->getGeneListHandler()[sortedIndexes[i]].first;
			outputStream << " " << geneSymbol << "("
					<< 100.0 * aggregation[sortedIndexes[i]]
							/ (double) numberOfPatients << "%) ";
		}
		outputStream << "}" << std::endl << std::endl;
	}
}

void TCGADataNormalizer::printMostExpressedGenesByClass(
		unsigned int maxNumberGenes, const std::string &filename) {
	std::ofstream outputStream(EXPORT_DIRECTORY + filename);
	if (verbose) {
		std::cout << std::endl << "****** FINDING MOST EXPRESSED GENES ******"
				<< std::endl;
	}
	for (const auto &kv : ptrToData->getRNASeqDataHandler(true)) {
		std::string cancer = kv.first;
		printMostExpressedGenesByClassUtility(outputStream, maxNumberGenes,
				cancer, false);
		printMostExpressedGenesByClassUtility(outputStream, maxNumberGenes,
				cancer, true);
	}
}

void TCGADataNormalizer::exportToFile(double positiveValue,
		double negativeValue) {
	ptrToData->transposeData(false);
	std::ofstream outputStreamSamples(HEINZ_SAMPLES_LIST);
	const auto &dataMatrix = ptrToData->getDataMatrixHandler();
	unsigned int N = dataMatrix.cols();
	for (unsigned int i = 0; i < N; ++i) {

		std::string outputFilename =
				ptrToData->getSamplesHandler()[i].toFullString();
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
