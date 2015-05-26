#include "unsupervisedNormalization.hpp"

#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>

#include "typedefs.hpp"
#include "k_means.hpp"
#include "utilities.hpp"

using namespace std;

void individualNormalization(Data &data, const string &cancer, int patientId,
		bool isTumor, const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters) {

	unsigned int numberOfGenes = data.getNumberOfProteins();
	vector<double> dataToNormalize(numberOfGenes);
	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		if (isTumor) {
			dataToNormalize[i] = data.tumorRNASeqData[cancer][i][patientId];
		} else {
			dataToNormalize[i] = data.controlRNASeqData[cancer][i][patientId];
		}
	}

	if (method == UnsupervisedNormalizationMethod::KMEANS) {
		vector<int> clusters(numberOfGenes, 0);
		K_Means<double> kMeans(dataToNormalize, clusters, parameters.K,
				parameters.Nmax, distanceDouble, addToDouble,
				divideDoubleByConstant, 0.0);
		kMeans.compute();
		for (unsigned int i = 0; i < numberOfGenes; ++i) {
			if (isTumor) {
				data.tumorRNASeqData[cancer][i][patientId] =
						(double) clusters[i];
			} else {
				data.controlRNASeqData[cancer][i][patientId] =
						(double) clusters[i];
			}
		}
	}

	else if (method
			== UnsupervisedNormalizationMethod::BINARY_ITERATED_KMEANS) {
		vector<int> clusters(numberOfGenes, 0);
		K_Means<double> kMeans(dataToNormalize, clusters, 2, 100,
				distanceDouble, addToDouble, divideDoubleByConstant, 0.0);
		kMeans.computeIteratedBinaryKMeans(parameters.Niteration);
		for (unsigned int i = 0; i < numberOfGenes; ++i) {
			if (isTumor) {
				data.tumorRNASeqData[cancer][i][patientId] =
						(double) clusters[i];
			} else {
				data.controlRNASeqData[cancer][i][patientId] =
						(double) clusters[i];
			}
		}
	}

	else if (method == UnsupervisedNormalizationMethod::BINARY_QUANTILE) {
		vector<size_t> rank = get_rank_increasing(dataToNormalize);
		for (unsigned int i = 0; i < numberOfGenes; ++i) {
			double p = (double) rank[i] / (double) numberOfGenes;
			if (p >= parameters.cutPercentage) {
				if (isTumor) {
					data.tumorRNASeqData[cancer][i][patientId] = 1.0;
				} else {
					data.controlRNASeqData[cancer][i][patientId] = 1.0;
				}
			} else {
				if (isTumor) {
					data.tumorRNASeqData[cancer][i][patientId] = 0.0;
				} else {
					data.controlRNASeqData[cancer][i][patientId] = 0.0;
				}
			}
		}
	}
}

void unsupervisedNormalization(Data &data,
		const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters) {
	cout << "Normalizing control Data..." << endl;
	for (auto &kv : data.controlRNASeqData) {
		string cancer = kv.first;
		cout << "\t" << cancer << "... " << flush;
		for (unsigned int patientId = 0; patientId < kv.second[0].size();
				++patientId) {
			individualNormalization(data, cancer, patientId, false, method,
					parameters);
		}
		cout << "Done." << endl;
	}
	cout << "Normalizing tumor Data..." << endl;
	for (auto &kv : data.tumorRNASeqData) {
		string cancer = kv.first;
		cout << "\t" << cancer << "... " << flush;
		for (unsigned int patientId = 0; patientId < kv.second[0].size();
				++patientId) {
			individualNormalization(data, cancer, patientId, true, method,
					parameters);
		}
		cout << "Done." << endl;
	}
}

/*
 *
 * PRINTS MOST EXPRESSED GENES PER CLASS
 *
 */

void RNASeqPrintMostExpressedGenes(const RNASeqData &rnaData,
		const string &cancer, bool isTumor, const GeneList &geneList,
		unsigned int maxNumberGenes, ofstream &outputStream) {

	unsigned int numberOfPatients = rnaData.at(cancer)[0].size();

	if (numberOfPatients > 0) {
		outputStream << cancer;
		if (isTumor) {
			outputStream << "-Tumor : ";
		} else {
			outputStream << "-Control : ";
		}

		outputStream << numberOfPatients << " patients. " << endl;

		unsigned int numberOfGenes = rnaData.at(cancer).size();
		vector<double> aggregation(numberOfGenes);
		fill(aggregation.begin(), aggregation.end(), 0.0);
		for (unsigned int i = 0; i < numberOfGenes; i++) {
			for (unsigned int j = 0; j < numberOfPatients; ++j) {
				aggregation[i] += rnaData.at(cancer)[i][j];
			}
		}

		outputStream << "Most expressed genes in the class : {";
		vector<size_t> sortedIndexes = sort_indexes_decreasing(aggregation);
		for (unsigned int i = 0; i < maxNumberGenes; ++i) {
			string geneSymbol = geneList[sortedIndexes[i]].first;
			outputStream << " " << geneSymbol << "("
					<< 100.0 * aggregation[sortedIndexes[i]]
							/ (double) numberOfPatients << "%) ";
		}
		outputStream << "}" << endl << endl;
	}
}

//Assumes binary normalization
void printMaxExpressedGenes(const Data &data, unsigned int maxNumberGenes,
		const string &filename) {

	ofstream outputStream("export/" + filename);

	cout << endl << "****** FINDING MOST EXPRESSED GENES ******" << endl;
	for (const auto &kv : data.tumorRNASeqData) {
		string cancer = kv.first;
		RNASeqPrintMostExpressedGenes(data.controlRNASeqData, cancer, false,
				data.geneList, maxNumberGenes, outputStream);
		RNASeqPrintMostExpressedGenes(data.tumorRNASeqData, cancer, true,
				data.geneList, maxNumberGenes, outputStream);
	}
}

