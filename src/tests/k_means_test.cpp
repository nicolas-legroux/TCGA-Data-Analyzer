#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "k_means_test.hpp"
#include "../dataReader.hpp"
#include "../k_means.hpp"

#include "../distanceMatrix.hpp"
#include "../utilities.hpp"
#include "../unsupervisedNormalization.hpp"
#include "../typedefs.hpp"

using namespace std;

void kMeansTest1(int K, int Nmax, string cancerName, int patientId) {
	vector<string> cancers { cancerName };
	Data data;
	readData(cancers, data, 0, patientId + 1);
	vector<double> dataToClusterSTL(
			data.getPatientTumorData(cancerName, patientId));

	MatrixX dataToCluster(1, dataToClusterSTL.size());
	for (unsigned int i = 0; i < dataToClusterSTL.size(); ++i) {
		dataToCluster(0, i) = dataToClusterSTL[i];
	}

	ClusteringParameters kMeansParameters;
	kMeansParameters.setKMeansParameters(K, Nmax, true);

	K_Means kMeans(dataToCluster, kMeansParameters);

	std::vector<int> clusters = kMeans.compute();
	MatrixX medoids = kMeans.getMedoids();

	std::vector<int> clusterCount(K, 0);
	for (int i : clusters) {
		clusterCount[i]++;
	}

	for (int i = 0; i != K; ++i) {
		double d = medoids(0, i);
		std::cout << "Cluster " << (i + 1) << ": " << d << ", size="
				<< clusterCount[i] << std::endl;
	}
}

void twodimensionalKmeans_test() {
	unsigned int K = 2;
	int Nmax = 10;
	unsigned int dim = 2;
	unsigned int N = 5;
	MatrixX dataToCluster(dim, N);
	dataToCluster << 0, 1, 0, 8, 8, 0, 1, 1, 7, 8;

	std::cout << "Clustering the following matrix : " << std::endl
			<< dataToCluster << std::endl;

	ClusteringParameters kMeansParameters;
	kMeansParameters.setKMeansParameters(K, Nmax, false);

	K_Means kMeans(dataToCluster, kMeansParameters);
	vector<int> clusters = kMeans.compute();
	MatrixX medoids = kMeans.getMedoids();
	std::cout << "The medoids are : " << std::endl << medoids;
	std::cout << std::endl << "The cluster asignments are : " << std::endl;
	print_vector(clusters);
}
