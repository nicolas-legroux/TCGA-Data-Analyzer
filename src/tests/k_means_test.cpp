#include <vector>
#include <string>
#include <iostream>
#include "../dataReader.hpp"
#include "../k_means.hpp"
#include "k_means_test.hpp"

using namespace std;

void kMeansTest1(int K, int Nmax, string cancerName, int patientId) {

	vector<string> cancers { cancerName };
	PatientList patientControlList;
	PatientList patientTumorList;
	RNASeqData controlData;
	RNASeqData tumorData;
	GeneList geneMapping(
			makeGeneMapping(
					"data/BRCA-normalized/TCGA-A1-A0SJ-01.genes.normalized.results"));
	readPatientData(cancers, patientControlList, patientTumorList);
	readRNASeqData(patientControlList, patientTumorList, geneMapping,
			controlData, tumorData, patientId + 1);

	int numberOfProteins = geneMapping.size();

	vector<double> data(numberOfProteins);

	for (int i = 0; i < numberOfProteins; ++i) {
		data[i] = tumorData.at(cancerName).at(i).at(patientId);
	}

	std::vector<int> clusters(data.size(), 0);

	std::vector<double> means = computeKMeans(data, clusters, K, Nmax);

	std::vector<int> clusterCount(K, 0);
	for (int i : clusters) {
		clusterCount[i]++;
	}

	for (int i = 0; i != K; ++i) {
		std::cout << "Cluster " << (i + 1) << ": " << means[i] << ", size="
				<< clusterCount[i] << std::endl;
	}
}

void iteratedBinaryKMeans_test(int N_iter, string cancerName, int patientId) {
	vector<string> cancers { cancerName };
	PatientList patientControlList;
	PatientList patientTumorList;
	RNASeqData controlData;
	RNASeqData tumorData;
	GeneList geneMapping(
			makeGeneMapping(
					"data/BRCA-normalized/TCGA-A1-A0SJ-01.genes.normalized.results"));
	readPatientData(cancers, patientControlList, patientTumorList);
	readRNASeqData(patientControlList, patientTumorList, geneMapping,
			controlData, tumorData, patientId + 1);

	int numberOfProteins = geneMapping.size();

	std::vector<double> data(numberOfProteins);

	for (int i = 0; i < numberOfProteins; ++i) {
		data[i] = tumorData.at(cancerName).at(i).at(patientId);
	}
	std::vector<int> clusters(data.size(), 0);

	iteratedBinaryKMeans(data, clusters, N_iter);

	std::vector<int> clusterCount(2, 0);
	for (int i : clusters) {
		clusterCount[i]++;
	}

	for (int i = 0; i != 2; ++i) {
		std::cout << "Cluster " << (i + 1) << ": size=" << clusterCount[i]
				<< std::endl;
	}
}
