#include <vector>
#include "clustering.hpp"
#include "utilities.hpp"
#include "typedefs.hpp"
#include "k_means.hpp"
#include "hierarchical_clustering.hpp"

using std::vector;
using std::map;
using std::string;
using std::cout;
using std::endl;

vector<int> getRealClusters(vector<SampleIdentifier> &sampleIdentifiers) {
	vector<int> realClusters;
	SampleIdentifier prev = sampleIdentifiers[0];
	int currentCluster = 0;
	realClusters.push_back(currentCluster);

	for (auto it = sampleIdentifiers.begin() + 1; it != sampleIdentifiers.end();
			++it) {
		SampleIdentifier next = *it;
		if (prev.cancerName != next.cancerName
				|| prev.isTumor != next.isTumor) {
			++currentCluster;
		}
		realClusters.push_back(currentCluster);
		prev = next;
	}

	return realClusters;
}

map<int, string> getRealLabelsMap(
		std::vector<SampleIdentifier> &sampleIdentifiers) {

	map<int, string> labelsMap;
	SampleIdentifier prev = sampleIdentifiers[0];
	string label = prev.cancerName + "-"
			+ ((prev.isTumor) ? "Tumor" : "Control");
	int currentCluster = 0;
	labelsMap.insert(std::make_pair(currentCluster, label));

	for (auto it = sampleIdentifiers.begin() + 1; it != sampleIdentifiers.end();
			++it) {
		SampleIdentifier next = *it;
		if (prev.cancerName != next.cancerName
				|| prev.isTumor != next.isTumor) {
			++currentCluster;
			string label = next.cancerName + "-"
					+ ((next.isTumor) ? "Tumor" : "Control");
			labelsMap.insert(std::make_pair(currentCluster, label));
		}
		prev = next;
	}
	return labelsMap;
}

vector<int> cluster_KMeans(const vector<vector<double>> &data, int K,
		int Nmax, bool verbose) {
	vector<int> clusters(data.size(), 0);
	unsigned int n = data[0].size();

	K_Means<vector<double>> kMeans(data, clusters, K, Nmax,
			EuclideanSpace<double>(n), verbose);

	kMeans.compute();
	return clusters;
}

vector<int> cluster_Hierarchical(const vector<double> &matrix,
		const DistanceMetric &distanceMetric,
		const LinkageMethod &linkageMethod, int K, bool verbose) {

	MatrixType matrixType;
	switch (distanceMetric) {
	case DistanceMetric::PEARSON_CORRELATION:
		matrixType = MatrixType::SIMILARITY;
		break;
	case DistanceMetric::SPEARMAN_CORRELATION:
		matrixType = MatrixType::SIMILARITY;
		break;
	case DistanceMetric::EUCLIDEAN_DISTANCE:
		matrixType = MatrixType::DISTANCE;
		break;
	case DistanceMetric::MANHATTAN_DISTANCE:
		matrixType = MatrixType::DISTANCE;
		break;
	case DistanceMetric::COSINE_SIMILARITY:
		matrixType = MatrixType::SIMILARITY;
		break;
	default:
		throw std::invalid_argument("Unknown distance measure.");
	}

	Hierarchical_Clustering hierarchicalClustering(matrix, linkageMethod,
			matrixType, verbose);
	vector<int> clusters = hierarchicalClustering.compute(K);

	return clusters;
}

/*
 *
 * CLUSTERING EVALUATION
 *
 */

double randIndex(const vector<int> &clustering1,
		const vector<int> &clustering2) {
	assert(clustering1.size() == clustering2.size());
	unsigned int n = clustering1.size();
	int count = 0;
	for (unsigned int i = 0; i < n - 1; ++i) {
		for (unsigned int j = i + 1; j < n; ++j) {
			if ((clustering1[i] == clustering1[j])
					&& clustering2[i] == clustering2[j]) {
				count++;
			} else if ((clustering1[i] != clustering1[j])
					&& clustering2[i] != clustering2[j]) {
				count++;
			}
		}
	}

	double number_of_pairs = (double) n * ((double) n - 1.0) * 0.5;
	return (double) count / number_of_pairs;
}

//Assumes the cluster labels range from 0 to (K-1)
double adjustedRandIndex(const std::vector<int> &clustering1,
		const std::vector<int> &clustering2) {

	assert(clustering1.size() == clustering2.size());
	unsigned int n = clustering1.size();

	int K1 = *max_element(clustering1.cbegin(), clustering1.cend()) + 1;
	int K2 = *max_element(clustering2.cbegin(), clustering2.cend()) + 1;

	vector<int> contingencyTable(K1 * K2, 0);
	vector<int> count1(K1, 0);
	vector<int> count2(K2, 0);

	for (unsigned int i = 0; i != n; ++i) {
		int c1 = clustering1[i];
		int c2 = clustering2[i];
		++contingencyTable[c1 * K2 + c2];
		++count1[c1];
		++count2[c2];
	}

	transform(contingencyTable.begin(), contingencyTable.end(),
			contingencyTable.begin(), numberOfPairs);
	transform(count1.begin(), count1.end(), count1.begin(), numberOfPairs);
	transform(count2.begin(), count2.end(), count2.begin(), numberOfPairs);

	double index = (double) accumulate(contingencyTable.cbegin(),
			contingencyTable.cend(), 0);
	double temp1 = (double) accumulate(count1.cbegin(), count1.cend(), 0.0);
	double temp2 = (double) accumulate(count2.cbegin(), count2.cend(), 0.0);
	double maxIndex = 0.5 * (temp1 + temp2);
	double expectedIndex = temp1 * temp2 / (double) numberOfPairs(n);

	return (index - expectedIndex) / (maxIndex - expectedIndex);
}

void printClustering(const map<int, string> &labelsMap,
		const vector<int> &realClusters, const vector<int> &computedClusters) {

	cout << endl << "****** Showing clustering results : ******" << endl << endl;;
	unsigned int realClustersN = labelsMap.size();
	unsigned int computedClustersN = *std::max_element(
			computedClusters.cbegin(), computedClusters.cend()) + 1;
	vector<int> clusteringGraph(realClustersN * computedClustersN, 0);

	assert(realClusters.size() == computedClusters.size());

	for (unsigned int i = 0; i < realClusters.size(); ++i) {
		int realCluster = realClusters[i];
		int computedCluster = computedClusters[i];
		clusteringGraph[realCluster * computedClustersN + computedCluster]++;
	}

	cout << "-----------\t";
	for (unsigned int j = 0; j < computedClustersN; ++j) {
		cout << "#" << j << "\t";
	}
	cout << "SUM" << endl;

	for (unsigned int i = 0; i < realClustersN; ++i) {
		cout << labelsMap.at(i) << "\t";
		unsigned int sum = 0;
		for (unsigned int j = 0; j < computedClustersN; ++j) {
			sum += clusteringGraph[i * computedClustersN + j];
			cout << clusteringGraph[i * computedClustersN + j] << "\t";
		}
		cout << sum << endl;
	}

	unsigned int global_sum = 0;
	cout << "SUM" << "\t\t";
	for (unsigned int j = 0; j < computedClustersN; ++j) {
		unsigned int sum = 0;

		for (unsigned int i = 0; i < realClustersN; ++i) {
			sum += clusteringGraph[i * computedClustersN + j];
		}
		global_sum += sum;
		cout << sum << "\t";
	}

	cout << global_sum << endl;
}
