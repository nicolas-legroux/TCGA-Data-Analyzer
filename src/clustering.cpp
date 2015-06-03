#include <vector>
#include "clustering.hpp"
#include "utilities.hpp"
#include "typedefs.hpp"
#include "k_means.hpp"
#include "hierarchical_clustering.hpp"

using namespace std;

vector<int> getRealClusters(vector<SampleIdentifier> &sampleIdentifiers) {
	vector<int> realClusters;
	SampleIdentifier prev = sampleIdentifiers[0];
	int currentCluster = 0;
	cout << currentCluster << " -> " << prev.toString() << endl;
	realClusters.push_back(currentCluster);

	for (auto it = sampleIdentifiers.begin() + 1; it != sampleIdentifiers.end();
			++it) {
		SampleIdentifier next = *it;
		if (prev.cancerName != next.cancerName
				|| prev.isTumor != next.isTumor) {
			++currentCluster;
			cout << currentCluster << " -> " << next.toString() << endl;
		}
		realClusters.push_back(currentCluster);
		prev = next;
	}

	return realClusters;
}

vector<int> cluster_KMeans(const vector<vector<double>> &data, int K,
		int Nmax) {
	vector<int> clusters(data.size(), 0);
	K_Means<vector<double>> kMeans(data, clusters, K, Nmax, EuclideanSpace<double>(data[0].size()));

	kMeans.compute();
	return clusters;
}

vector<int> cluster_Hierarchical(const vector<double> &matrix,
		const DistanceMetric &distanceMetric,
		const LinkageMethod &linkageMethod, int K) {

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
		throw invalid_argument("Unknown distance measure.");
	}

	Hierarchical_Clustering hierarchicalClustering(matrix, linkageMethod, matrixType);
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
