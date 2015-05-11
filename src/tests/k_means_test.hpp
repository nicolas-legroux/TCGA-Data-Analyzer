/*
 * k_means_test.hpp
 *
 *  Created on: May 4, 2015
 *      Author: nicolas
 */

#ifndef SRC_TESTS_K_MEANS_TEST_HPP_
#define SRC_TESTS_K_MEANS_TEST_HPP_

#include <vector>
#include <string>
#include <iostream>
#include "../dataReader.hpp"
#include "../k_means.hpp"

void kMeansTest1(int K, int Nmax) {

	std::string cancerName = "BRCA";
	std::vector<std::string> cancers{cancerName};
	PatientList patientControlList;
	PatientList patientTumorList;
	RNASeqData controlData;
	RNASeqData tumorData;
	GeneList geneMapping(
				makeGeneMapping(
						"data/BRCA-normalized/TCGA-A1-A0SJ-01.genes.normalized.results"));
	readPatientData(cancers, patientControlList, patientTumorList);
	readRNASeqData(patientControlList, patientTumorList, geneMapping,
				controlData, tumorData, 100);

	int numberOfProteins = geneMapping.size();

	std::vector<double> data(numberOfProteins);

	for (int i = 0; i < numberOfProteins; ++i) {
		data[i] = tumorData.at(cancerName).at(i).at(54);
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

#endif /* SRC_TESTS_K_MEANS_TEST_HPP_ */
