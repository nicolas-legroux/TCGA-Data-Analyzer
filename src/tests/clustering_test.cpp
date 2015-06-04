#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "clustering_test.hpp"
#include "../dataReader.hpp"
#include "../clustering.hpp"
#include "../distanceMatrix.hpp"
#include "../unsupervisedNormalization.hpp"
#include "../utilities.hpp"

using std::vector;
using std::map;
using std::string;
using std::cout;
using std::endl;

void clustering_KMeans_test(const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters) {
	//STEP 1 : READ DATA
	string filenameCancers = "cancer.list";
	vector<string> cancers { "BRCA", "KIRC", "OV", "LUAD", "GBM", "THCA"};
	Data data;
	readData(cancers, data, 50, 300);

	//STEP 2 : NORMALIZE
	unsupervisedNormalization(data, method, parameters);

	//STEP3 : PREPARE DATA FOR CLUSTERING
	std::vector<std::vector<double>> transposedData;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;
	data.transposeData(transposedData, sampleIdentifiers, cancerPatientIDList);

	//STEP4 : CLUSTER
	map<int, string> labelsMap = getRealLabelsMap(sampleIdentifiers);
	cout << endl << "Labels Map : ";
	print_map(labelsMap);

	vector<int> realClusters = getRealClusters(sampleIdentifiers);
	print_vector(realClusters);

	vector<int> computedClusters = cluster_KMeans(transposedData, labelsMap.size(), 1000);
	print_vector(computedClusters);

	cout << "Adjusted Rand Index : "
			<< adjustedRandIndex(realClusters, computedClusters) << endl;
}

void clustering_Hierarchical_test(const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters,
		const DistanceMetric &distanceMetric,
		const LinkageMethod &linkageMethod) {
	//STEP 1 : READ DATA
	string filenameCancers = "cancer.list";
	vector<string> cancers { "BRCA", "KIRC", "OV", "LUAD", "GBM", "THCA" };
	Data data;
	readData(cancers, data, 50, 300);

	//STEP 2 : NORMALIZE
	unsupervisedNormalization(data, method, parameters);

	//STEP3 : PREPARE DATA FOR CLUSTERING
	std::vector<std::vector<double>> transposedData;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;
	data.transposeData(transposedData, sampleIdentifiers, cancerPatientIDList);

	//STEP4: COMPUTE DISTANCE MATRIX
	std::vector<double> distanceMatrix = computeDistanceMatrix(transposedData,
			distanceMetric);

	//STEP5 : CLUSTER
	vector<int> realClusters = getRealClusters(sampleIdentifiers);
	for_each(realClusters.cbegin(), realClusters.cend(),
			[](int i) {cout << i << " ";});
	cout << endl;

	vector<int> computedClusters = cluster_Hierarchical(distanceMatrix,
			distanceMetric, linkageMethod, 10);

	for_each(computedClusters.cbegin(), computedClusters.cend(),
			[](int i) {cout << i << " ";});

	cout << endl;

	cout << "Adjusted Rand Index : "
			<< adjustedRandIndex(realClusters, computedClusters) << endl;
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
	cout << "AdjustedRandIndex = "
			<< adjustedRandIndex(clustering1, clustering2);
}
