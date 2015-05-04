#include <unordered_map>
#include <string>
#include <iostream>

#include "k_means.hpp"
#include "patientUnsupervisedNormalization.hpp"

using namespace std;

typedef unordered_map<string, vector<vector<double>>> RNASeqData;

void normalizePatientKMeans(const string &cancer, RNASeqData &rnaData, int patientId, int K){
	int numberOfGenes = rnaData[cancer].size();
	vector<double> data(numberOfGenes);
	for(int i=0; i<numberOfGenes; ++i){
		data[i] = rnaData[cancer][i][patientId];
	}

	vector<int> clusters(numberOfGenes);
	int Nmax = 100;

	computeKMeans(data, clusters, K, Nmax);

	for(int i=0; i<numberOfGenes; ++i){
		rnaData[cancer][i][patientId] = (double) clusters[i];
	}
}

void normalizeKMeans(RNASeqData &controlData, RNASeqData &tumorData, int K){
	cout << "Normalizing control Data..." << endl;
	for(auto &mappedData : controlData){
		string cancer = mappedData.first;
		for(unsigned int i=0; i<mappedData.second[0].size(); ++i){
			normalizePatientKMeans(cancer, controlData, i, K);
		}
	}

	cout << "Normalizing tumor Data..." << endl;
	for(auto &mappedData : tumorData){
		string cancer = mappedData.first;
		for(unsigned int i=0; i<mappedData.second[0].size(); ++i){
			normalizePatientKMeans(cancer, tumorData, i, K);
		}
	}
}
