#include <unordered_map>
#include <string>
#include <iostream>

#include "typedefs.hpp"
#include "k_means.hpp"
#include "patientUnsupervisedNormalization.hpp"
#include "utilities.hpp"

using namespace std;

void normalizePatientKMeans(const string &cancer, RNASeqData &rnaData,
		int patientId, int K, int Nmax) {
	int numberOfGenes = rnaData[cancer].size();
	vector<double> data(numberOfGenes);
	for (int i = 0; i < numberOfGenes; ++i) {
		data[i] = rnaData[cancer][i][patientId];
	}

	vector<int> clusters(numberOfGenes);

	computeKMeans(data, clusters, K, Nmax);

	for (int i = 0; i < numberOfGenes; ++i) {
		rnaData[cancer][i][patientId] = (double) clusters[i];
	}
}

void normalizeKMeans(RNASeqData &controlData, RNASeqData &tumorData, int K,
		int Nmax) {
	cout << "Normalizing control Data..." << endl;
	for (auto &mappedData : controlData) {
		string cancer = mappedData.first;
		cout << "\t" << cancer << "... " << flush;
		for (unsigned int i = 0; i < mappedData.second[0].size(); ++i) {
			normalizePatientKMeans(cancer, controlData, i, K, Nmax);
		}
		cout << "Done." << endl;
	}

	cout << "Normalizing tumor Data..." << endl;
	for (auto &mappedData : tumorData) {
		string cancer = mappedData.first;
		cout << "\t" << cancer << "... " << flush;
		for (unsigned int i = 0; i < mappedData.second[0].size(); ++i) {
			normalizePatientKMeans(cancer, tumorData, i, K, Nmax);
		}
		cout << "Done." << endl;
	}
}

void normalizePatientQuantile(const string &cancer, RNASeqData &rnaData,
		int patientId, double cutPercentage) {

	int numberOfGenes = rnaData[cancer].size();
	vector<double> data(numberOfGenes);
	for (int i = 0; i < numberOfGenes; ++i) {
		data[i] = rnaData[cancer][i][patientId];
	}

	vector<size_t> rank = get_rank_increasing(data);
	for(int i=0; i< numberOfGenes; ++i){
		double p = (double)rank[i]/(double)numberOfGenes;
		if(p>=cutPercentage){
			rnaData[cancer][i][patientId] = 1.0;
		}
		else{
			rnaData[cancer][i][patientId] = 0.0;
		}
	}
}

void normalizeQuantile(RNASeqData &controlData, RNASeqData &tumorData,
		double cutPercentage) {
	cout << "Normalizing control Data..." << endl;
	for (auto &mappedData : controlData) {
		string cancer = mappedData.first;
		cout << "\t" << cancer << "... " << flush;
		for (unsigned int i = 0; i < mappedData.second[0].size(); ++i) {
			normalizePatientQuantile(cancer, controlData, i, cutPercentage);
		}
		cout << "Done." << endl;
	}

	cout << "Normalizing tumor Data..." << endl;
	for (auto &mappedData : tumorData) {
		string cancer = mappedData.first;
		cout << "\t" << cancer << "... " << flush;
		for (unsigned int i = 0; i < mappedData.second[0].size(); ++i) {
			normalizePatientQuantile(cancer, tumorData, i, cutPercentage);
		}
		cout << "Done." << endl;
	}
}

void printMaxExpressedGenes(const RNASeqData &controlNormalized,
		const RNASeqData &tumorNormalized, const GeneList &geneList) {

	cout << endl << "****** FINDING MOST EXPRESSED GENES ******" << endl;

	for (const auto &kv : controlNormalized) {
		string cancer = kv.first;
		int numberOfPatients = kv.second[0].size();

		if (numberOfPatients > 0) {
			cout << cancer << "-Control : " << numberOfPatients
					<< " patients. Aggregating the scores... " << flush;

			int numberOfGenes = kv.second.size();
			vector<double> aggregation(numberOfGenes);
			fill(aggregation.begin(), aggregation.end(), 0.0);
			for (int i = 0; i < numberOfGenes; i++) {
				for (int j = 0; j < numberOfPatients; ++j) {
					aggregation[i] += kv.second[i][j];
				}
			}

			cout << "Done." << endl;
			cout << "Most expressed genes in the class : {";
			vector<size_t> sortedIndexes = sort_indexes_decreasing(aggregation);
			for (int i = 0; i < 10; ++i) {
				string geneSymbol = geneList[sortedIndexes[i]].first;
				cout << " " << geneSymbol << "("
						<< aggregation[sortedIndexes[i]] << ") ";
			}
			cout << "}" << endl << endl;
		}
	}

	for (const auto &kv : tumorNormalized) {
		string cancer = kv.first;
		int numberOfPatients = kv.second[0].size();

		if (numberOfPatients > 0) {
			cout << cancer << "-Tumor : " << numberOfPatients
					<< " patients. Aggregating the scores... " << flush;

			int numberOfGenes = kv.second.size();
			vector<double> aggregation(numberOfGenes);
			fill(aggregation.begin(), aggregation.end(), 0.0);
			for (int i = 0; i < numberOfGenes; i++) {
				for (int j = 0; j < numberOfPatients; ++j) {
					aggregation[i] += kv.second[i][j];
				}
			}

			cout << "Done." << endl;
			cout << "Most expressed genes in the class : {";
			vector<size_t> sortedIndexes = sort_indexes_decreasing(aggregation);
			for (int i = 0; i < 10; ++i) {
				string geneSymbol = geneList[sortedIndexes[i]].first;
				cout << " " << geneSymbol << "("
						<< aggregation[sortedIndexes[i]] << ") ";
			}
			cout << "}" << endl << endl;
		}
	}
}
