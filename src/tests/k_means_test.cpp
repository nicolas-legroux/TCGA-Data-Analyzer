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
	std::vector<int> clusters(dataToClusterSTL.size(), 0);

	MatrixX dataToCluster(1, dataToClusterSTL.size());
	for (unsigned int i = 0; i < dataToClusterSTL.size(); ++i) {
		dataToCluster(0, i) = dataToClusterSTL[i];
	}

	ClusteringParameters kMeansParameters;
	kMeansParameters.setKMeansParameters(K, Nmax, true);

	K_Means kMeans(dataToCluster, kMeansParameters, clusters);

	MatrixX medoids = kMeans.compute();

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

void iteratedBinaryKMeans_test(int N_iter, string cancerName, int patientId) {
	vector<string> cancers { cancerName };
	Data data;
	readData(cancers, data, 0, patientId + 1);
	vector<double> dataToClusterSTL(
			data.getPatientTumorData(cancerName, patientId));
	std::vector<int> clusters(dataToClusterSTL.size(), 0);

	MatrixX dataToCluster(1, dataToClusterSTL.size());
	for (unsigned int i = 0; i < dataToClusterSTL.size(); ++i) {
		dataToCluster(0, i) = dataToClusterSTL[i];
	}

	ClusteringParameters kMeansParameters;
	kMeansParameters.setKMeansParameters(2, 100, true);

	K_Means kMeans(dataToCluster, kMeansParameters, clusters);

	kMeans.computeIteratedBinaryKMeans(N_iter);

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
	int Nmax = 10;
	unsigned int dim = 2;
	unsigned int N = 5;
	MatrixX dataToCluster(dim, N);
	dataToCluster << 0, 1, 0, 8, 8, 0, 1, 1, 7, 8;

	std::cout << "Clustering the following matrix : " << std::endl
			<< dataToCluster << std::endl;
	vector<int> clusters(N, 0);

	ClusteringParameters kMeansParameters;
	kMeansParameters.setKMeansParameters(K, Nmax, true);

	K_Means kMeans(dataToCluster, kMeansParameters, clusters);

	MatrixX medoids = kMeans.compute();
	std::cout << "The medoids are : " << medoids;
	std::cout << std::endl << "The cluster asignments are : " << std::endl;
	print_vector(clusters);
}
