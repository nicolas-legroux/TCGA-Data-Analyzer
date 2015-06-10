#ifndef SRC_CLUSTERING_HPP_
#define SRC_CLUSTERING_HPP_

#include <vector>
#include <map>
#include "typedefs.hpp"
#include "hierarchical_clustering.hpp"
#include "distanceMatrix.hpp"

enum ClusteringMethod {
	KMEANS_CLUSTERING, HIERARCHICAL_CLUSTERING, SPECTRAL_CLUSTERING
};

struct ClusteringParameters {
	//Number of clusters
	unsigned int K = 2;

	//For K-Means : maximum number of iterations
	unsigned int maxIterations = 1000;

	//For Hierarchical and spectral clustering : matrix type
	DistanceMetric distanceMetric;

	//For hierarchical clustering: Linkage method
	LinkageMethod linkageMethod = LinkageMethod::COMPLETE;

	bool verbose = false;

	void setKMeansParameters(unsigned int _K,
			unsigned int _maxIterations = 1000, bool _verbose = false) {
		K = _K;
		maxIterations = _maxIterations;
		verbose = _verbose;
	}

	void setHierarchicalParameters(unsigned int _K,
			DistanceMetric _distanceMetric, LinkageMethod _linkageMethod =
					LinkageMethod::COMPLETE, bool _verbose = false) {
		K = _K;
		distanceMetric = _distanceMetric;
		linkageMethod = _linkageMethod;
		verbose = _verbose;
	}

	void setSpectralParameters(unsigned int _K, DistanceMetric _distanceMetric,
			bool _verbose = false) {
		K = _K;
		distanceMetric = _distanceMetric;
		verbose = _verbose;
	}
};

// Function to get real labels
std::vector<int> getRealClusters(
		std::vector<SampleIdentifier> &sampleIdentifiers);

std::map<int, std::string> getRealLabelsMap(
		std::vector<SampleIdentifier> &sampleIdentifiers);

//Functions to compute clusters

std::vector<int> clusterRNASeqData(const Data &data,
		ClusteringMethod clusteringMethod,
		ClusteringParameters clusteringParameters);

std::vector<int> cluster_KMeans(const Data &data, ClusteringParameters clusteringParameters);
std::vector<int> cluster_Hierarchical(const MatrixX &matrix,
		const DistanceMetric &distanceMetric,
		const LinkageMethod &linkageMethod, int K, bool verbose = false);
std::vector<int> cluster_Spectral(const MatrixX &matrix,
		const DistanceMetric &distanceMetric, unsigned int K);

void printClustering(const std::map<int, std::string> &labelsMap,
		const std::vector<int> &realClusters,
		const std::vector<int> &computedClusters);

/*
 *
 * CLUSTERING EVALUATION
 *
 */

double randIndex(const std::vector<int> &clustering1,
		const std::vector<int> &clustering2);

double adjustedRandIndex(const std::vector<int> &clustering1,
		const std::vector<int> &clustering2);

#endif /* SRC_CLUSTERING_HPP_ */
