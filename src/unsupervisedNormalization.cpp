#include "unsupervisedNormalization.hpp"

#include <unordered_map>
#include <string>
#include <iostream>

#include "typedefs.hpp"
#include "k_means.hpp"
#include "utilities.hpp"

using namespace std;

/*
 *
 *  METHOD 1 : K MEANS
 *
 */

void individualNormalizationKMeans(const string &cancer, RNASeqData &rnaData,
		int patientId, int K, int Nmax) {
	unsigned int numberOfGenes = rnaData[cancer].size();
	vector<double> data(numberOfGenes);
	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		data[i] = rnaData[cancer][i][patientId];
	}

	vector<int> clusters(numberOfGenes);

	computeKMeans(data, clusters, K, Nmax);

	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		rnaData[cancer][i][patientId] = (double) clusters[i];
	}
}

void RNASeqDataNormalizationKMeans(RNASeqData &rnaData, int K, int Nmax) {
	for (auto &mappedData : rnaData) {
		string cancer = mappedData.first;
		cout << "\t" << cancer << "... " << flush;
		for (unsigned int i = 0; i < mappedData.second[0].size(); ++i) {
			individualNormalizationKMeans(cancer, rnaData, i, K, Nmax);
		}
		cout << "Done." << endl;
	}
}

void normalizeKMeans(RNASeqData &controlData, RNASeqData &tumorData, int K,
		int Nmax) {
	cout << "Normalizing control Data..." << endl;
	RNASeqDataNormalizationKMeans(controlData, K, Nmax);
	cout << "Normalizing tumor Data..." << endl;
	RNASeqDataNormalizationKMeans(tumorData, K, Nmax);
}

/*
 *
 *  METHOD 2 : QUANTILE
 *
 */

void individualNormalizationQuantile(const string &cancer, RNASeqData &rnaData,
		int patientId, double cutPercentage) {

	int numberOfGenes = rnaData[cancer].size();
	vector<double> data(numberOfGenes);
	for (int i = 0; i < numberOfGenes; ++i) {
		data[i] = rnaData[cancer][i][patientId];
	}

	vector<size_t> rank = get_rank_increasing(data);
	for (int i = 0; i < numberOfGenes; ++i) {
		double p = (double) rank[i] / (double) numberOfGenes;
		if (p >= cutPercentage) {
			rnaData[cancer][i][patientId] = 1.0;
		} else {
			rnaData[cancer][i][patientId] = 0.0;
		}
	}
}

void RNASeqDataNormalizationQuantile(RNASeqData &rnaData,
		double cutPercentage) {
	for (auto &mappedData : rnaData) {
		string cancer = mappedData.first;
		cout << "\t" << cancer << "... " << flush;
		for (unsigned int i = 0; i < mappedData.second[0].size(); ++i) {
			individualNormalizationQuantile(cancer, rnaData, i, cutPercentage);
		}
		cout << "Done." << endl;
	}
}

void normalizeQuantile(RNASeqData &controlData, RNASeqData &tumorData,
		double cutPercentage) {
	cout << endl << "Normalizing control Data..." << endl;
	RNASeqDataNormalizationQuantile(controlData, cutPercentage);

	cout << endl << "Normalizing tumor Data..." << endl;
	RNASeqDataNormalizationQuantile(tumorData, cutPercentage);
}

/*
 *
 * PRINTS MOST EXPRESSED GENES PER CLASS
 *
 */

void RNASeqPrintMostExpressedGenes(const RNASeqData &rnaData,
		const string &cancer, bool isTumor, const GeneList &geneList,
		unsigned int maxNumberGenes) {

	unsigned int numberOfPatients = rnaData.at(cancer)[0].size();

	if (numberOfPatients > 0) {
		cout << cancer;
		if (isTumor) {
			cout << "-Tumor : ";
		} else {
			cout << "-Control : ";
		}

		cout << numberOfPatients << " patients. " << endl;

		unsigned int numberOfGenes = rnaData.at(cancer).size();
		vector<double> aggregation(numberOfGenes);
		fill(aggregation.begin(), aggregation.end(), 0.0);
		for (unsigned int i = 0; i < numberOfGenes; i++) {
			for (unsigned int j = 0; j < numberOfPatients; ++j) {
				aggregation[i] += rnaData.at(cancer)[i][j];
			}
		}

		cout << "Most expressed genes in the class : {";
		vector<size_t> sortedIndexes = sort_indexes_decreasing(aggregation);
		for (unsigned int i = 0; i < maxNumberGenes; ++i) {
			string geneSymbol = geneList[sortedIndexes[i]].first;
			cout << " " << geneSymbol << "("
					<< 100.0 * aggregation[sortedIndexes[i]]
							/ (double) numberOfPatients << "%) ";
		}
		cout << "}" << endl << endl;
	}
}

//Assumes binary normalization
void printMaxExpressedGenes(const RNASeqData &controlNormalized,
		const RNASeqData &tumorNormalized, const GeneList &geneList,
		unsigned int maxNumberGenes) {

	cout << endl << "****** FINDING MOST EXPRESSED GENES ******" << endl;
	for(const auto &kv : tumorNormalized){
		string cancer = kv.first;
		RNASeqPrintMostExpressedGenes(controlNormalized, cancer, false, geneList,
				maxNumberGenes);
		RNASeqPrintMostExpressedGenes(tumorNormalized, cancer, true, geneList,
				maxNumberGenes);
	}
}
