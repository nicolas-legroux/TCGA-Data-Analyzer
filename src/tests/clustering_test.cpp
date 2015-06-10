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
#include "../spectral_clustering.hpp"

using std::vector;
using std::map;
using std::string;
using std::cout;
using std::endl;

void clustering_KMeans_test(const vector<string> &cancers, int maxControl,
		int maxTumor, const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters) {
	//STEP 1 : READ DATA
	Data data;
	readData(cancers, data, maxControl, maxTumor);

	//STEP 2 : NORMALIZE
	unsupervisedNormalization(data, method, parameters);

	//STEP3 : PREPARE DATA FOR CLUSTERING
	data.transposeData();

	//STEP4 : CLUSTER
	map<int, string> labelsMap = getRealLabelsMap(data.sampleIdentifiers);
	vector<int> realClusters = getRealClusters(data.sampleIdentifiers);
	ClusteringParameters kMeansParameters;
	kMeansParameters.setKMeansParameters(labelsMap.size(), 1000, true);
	vector<int> computedClusters = cluster_KMeans(data, kMeansParameters);

	printClustering(labelsMap, realClusters, computedClusters);

	cout << endl << "Adjusted Rand Index : "
			<< adjustedRandIndex(realClusters, computedClusters) << endl;
}

void clustering_Hierarchical_test(const vector<string> &cancers, int maxControl,
		int maxTumor, const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters,
		const DistanceMetric &distanceMetric,
		const LinkageMethod &linkageMethod) {

	// STEP 1 : READ DATA
	Data data;
	readData(cancers, data, maxControl, maxTumor);

	// STEP 2 : NORMALIZE
	unsupervisedNormalization(data, method, parameters);

	// STEP3 : PREPARE DATA FOR CLUSTERING
	data.transposeData();

	// STEP4: COMPUTE DISTANCE MATRIX
	MatrixX distanceMatrix = computeDistanceMatrix(
			data.transposedData, distanceMetric);

	// STEP5 : CLUSTER
	map<int, string> labelsMap = getRealLabelsMap(data.sampleIdentifiers);
	vector<int> realClusters = getRealClusters(data.sampleIdentifiers);
	vector<int> computedClusters = cluster_Hierarchical(distanceMatrix,
			distanceMetric, linkageMethod, labelsMap.size(), true);

	// STEP6 : PRINT RESULTS
	printClustering(labelsMap, realClusters, computedClusters);
	cout << endl << "Adjusted Rand Index : "
			<< adjustedRandIndex(realClusters, computedClusters) << endl;
}

void clustering_Spectral_test(const vector<string> &cancers, int maxControl,
		int maxTumor, const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters,
		const DistanceMetric &distanceMetric) {
	// STEP 1 : READ DATA
	Data data;
	readData(cancers, data, maxControl, maxTumor);

	// STEP 2 : NORMALIZE
	unsupervisedNormalization(data, method, parameters);

	// STEP3 : PREPARE DATA FOR CLUSTERING
	data.transposeData();

	// STEP4: COMPUTE DISTANCE MATRIX
	MatrixX distanceMatrix = computeDistanceMatrix(
			data.transposedData, distanceMetric);

	// STEP5 : CLUSTER
	SimilarityGraphTransformation similarityGraphTransformation =
			SimilarityGraphTransformation::K_NEAREST_NEIGHBORS;
	SpectralClusteringParameters similarityGraphTransformationParameters;
	similarityGraphTransformationParameters.setKNearestNeighborsParameters(3);
	map<int, string> labelsMap = getRealLabelsMap(data.sampleIdentifiers);
	vector<int> realClusters = getRealClusters(data.sampleIdentifiers);

	Spectral_Clustering spectralClustering(distanceMatrix, distanceMetric,
			similarityGraphTransformation,
			similarityGraphTransformationParameters);
	vector<int> computedClusters = spectralClustering.compute(labelsMap.size());
	;

	// STEP6 : PRINT RESULTS
	printClustering(labelsMap, realClusters, computedClusters);
	cout << endl << "Adjusted Rand Index : "
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
