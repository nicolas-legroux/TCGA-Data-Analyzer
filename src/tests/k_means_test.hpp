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

void kMeansTest1(const RNASeqData &rnaSeqData, std::string cancerName, int patientID) {
	int numberOfProteins = rnaSeqData.at(cancerName).size();
	std::vector<double> data(numberOfProteins);
	for (int i = 0; i < numberOfProteins; ++i) {
		data[i] = rnaSeqData.at(cancerName).at(i).at(patientID);
	}
	std::vector<int> clusters(data.size(), 0);

	int K = 8;
	double Nmax = 1000;

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
