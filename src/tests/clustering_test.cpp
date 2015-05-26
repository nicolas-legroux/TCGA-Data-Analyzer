#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "clustering_test.hpp"
#include "../dataReader.hpp"
#include "../heatMap.hpp"
#include "../clustering.hpp"

#include "../distanceMatrix.hpp"
#include "../unsupervisedNormalization.hpp"
#include "../k_means.hpp"

using namespace std;

void clustering_test() {
	/*
	string filenameCancers = "cancer.list";
	vector<string> cancers {"BRCA"};
	PatientList patientControlList;
	PatientList patientTumorList;
	RNASeqData controlData;
	RNASeqData tumorData;
	GeneList geneMapping(
			makeGeneMapping(
					"data/BRCA-normalized/TCGA-A1-A0SJ-01.genes.normalized.results"));
	readPatientData(cancers, patientControlList, patientTumorList);
	readRNASeqData(patientControlList, patientTumorList, geneMapping,
			controlData, tumorData, 500);

	vector<double> COL1A1;

	cout << geneMapping[4041].first << endl;

	//normalizeKMeans(controlData, tumorData, 2, 1000);

	std::vector<std::vector<double>> data;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;

	prepareData(data, sampleIdentifiers, cancerPatientIDList,
			patientControlList, patientTumorList, controlData, tumorData);

	for(int i=0; i<data.size(); ++i){
		COL1A1.push_back(data[i][4041]);
	}

	vector<double> col1a1Control;
	vector<double> col1a1Tumor;



	vector<int> clustersCOL1A1(data.size(), 0);
	computeKMeans(COL1A1, clustersCOL1A1, 2, 1000);
	for_each(clustersCOL1A1.cbegin(), clustersCOL1A1.cend(),
			[](int i) {cout << i << " ";});
	cout << endl;

	vector<int> realClusters = getRealClusters(sampleIdentifiers);
	for_each(realClusters.cbegin(), realClusters.cend(),
			[](int i) {cout << i << " ";});
	cout << endl;

	vector<int> computedClusters = clusterWithKMeans(data, 2, 1000);

	for_each(computedClusters.cbegin(), computedClusters.cend(),
			[](int i) {cout << i << " ";});

	cout << endl;

	cout << "Adjusted Rand Index : " << adjustedRandIndex(realClusters, computedClusters) << endl;
	cout << "Adjusted Rand Index : " << adjustedRandIndex(realClusters, clustersCOL1A1) << endl;
	cout << "Adjusted Rand Index : " << adjustedRandIndex(computedClusters, clustersCOL1A1) << endl;
	*/
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
