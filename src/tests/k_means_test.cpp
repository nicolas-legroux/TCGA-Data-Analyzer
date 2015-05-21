#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
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

void twodimensionalKmeans_test() {
	unsigned int K = 2;
	vector<double> vec1 { 0, 0 };
	vector<double> vec2 { 1, 1 };
	vector<double> vec3 { 0, 1 };
	vector<double> vec4 { 8, 7 };
	vector<double> vec5 { 8, 8 };
	vector<vector<double>> data { vec1, vec2, vec3, vec4, vec5 };
	vector<int> clusters(data.size(), 0);
	vector<vector<double>> means = computeKMeans(data, clusters, K, 10,
			euclidianNorm);
	std::vector<int> clusterCount(K, 0);
	for (int i : clusters) {
		clusterCount[i]++;
	}

	for (unsigned int i = 0; i != K; ++i) {
		cout << "Cluster " << (i + 1) << ": { ";
		for_each(means[i].cbegin(), means[i].cend(),
				[](double d) {cout << d << " ";});
		cout << "}, size=" << clusterCount[i] << endl;
	}
}

/* R Code :
 r1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
 r2 <- c(1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1)
 library(mclust)
 adjustedRandIndex(r1,r2)
 R output is 0.2737606
 */

void adjustedRandIndex_test() {
	vector<int> clustering1 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1 };
	vector<int> clustering2 { 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1 };
	cout << "AdjustedRandIndex = " << adjustedRandIndex(clustering1, clustering2);
}
