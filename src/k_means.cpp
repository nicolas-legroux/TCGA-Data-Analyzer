#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <iterator>
#include <assert.h>

#include "k_means.hpp"
#include "stats.hpp"
#include "utilities.hpp"

using namespace std;

/*
 *
 * ONE DIMENSIONAL K-MEANS
 *
 */

void initializeClusters(const vector<double> &data, vector<double> &means,
		int K, const vector<bool> &dataToCluster) {
	vector<double> copyData;
	copy_if_two_ranges(data.cbegin(), data.cend(), dataToCluster.cbegin(),
			back_inserter(copyData));
	auto it = unique(copyData.begin(), copyData.end());
	int N_unique = it - copyData.begin();
	int step = N_unique / (K + 1);
	for (int i = 0; i < K; ++i) {
		means[i] = copyData[(i + 1) * step];
	}
}

void initializeClustersRandomly(const vector<double> &data,
		vector<double> &means, int K, const vector<bool> &dataToCluster) {
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0, data.size() - 1);
	for (int i = 0; i < K; i++) {
		while (true) {
			int randomInt = distribution(generator);
			if (dataToCluster[randomInt]) {
				double randomData = data[randomInt];
				if (find(means.begin(), means.begin() + i, randomData)
						== (means.begin() + i)) {
					means[i] = randomData;
					break;
				}
			}
		}
	}
}

int findClosestClusterFromDataPoint(const vector<double> &means,
		double dataPoint, int K) {
	int closestCluster = 0;
	double closestDistance = abs(dataPoint - means[0]);
	for (int i = 1; i < K; ++i) {
		double d = abs(dataPoint - means[i]);
		if (d < closestDistance) {
			closestCluster = i;
			closestDistance = d;
		}
	}

	return closestCluster;
}

void recalculateMeans(const vector<double> &data, vector<double> &means,
		const vector<int> &clusters, int K) {
	vector<int> clusterSize(K, 0);
	fill(means.begin(), means.end(), 0.0);
	for (unsigned int i = 0; i != clusters.size(); ++i) {
		int cluster = clusters[i];
		if (cluster != -1) {
			++clusterSize[cluster];
			means[cluster] += data[i];
		}
	}
	for (int i = 0; i < K; ++i) {
		means[i] /= (double) clusterSize[i];
	}
}

bool kMeansIteration(const vector<double> &data, vector<double> &means,
		vector<int> &clusters, int K) {
	bool clustersChanged = false;

	//Assign clusters
	for (unsigned int i = 0; i != data.size(); ++i) {
		int oldCluster = clusters[i];
		if (oldCluster != -1) {
			int newCluster = findClosestClusterFromDataPoint(means, data[i], K);
			if (newCluster != clusters[i]) {
				clustersChanged = true;
				clusters[i] = newCluster;
			}
		}
	}

	//Recalculate means
	recalculateMeans(data, means, clusters, K);

	return clustersChanged;
}

void assignSortedClusters(vector<int> &clusters,
		const vector<size_t> &clusterRanks) {
	for (unsigned int i = 0; i != clusters.size(); ++i) {
		int currentCluster = clusters[i];
		if (currentCluster != -1) {
			clusters[i] = clusterRanks.at(currentCluster);
		}
	}
}

vector<double> computeKMeans(const vector<double> &data, vector<int> &clusters,
		int K, int Nmax) {

	//Initialize data Structures
	vector<double> means(K);
	vector<bool> dataToCluster(data.size());
	transform(clusters.cbegin(), clusters.cend(), dataToCluster.begin(),
			[](int cluster) {
				return (cluster != -1);
			});

	initializeClustersRandomly(data, means, K, dataToCluster);

	//Iterate
	int i = 0;
	while (i < Nmax && kMeansIteration(data, means, clusters, K)) {
		++i;
	}

	//Sort the clusters
	vector<size_t> clusterRanks(get_rank_increasing(means));
	assignSortedClusters(clusters, clusterRanks);
	sort(means.begin(), means.end());

	if (i == Nmax) {
		cout << "K-Means did not converge." << endl;
	}

	//cout << "K-Means finished after " << i << " iterations." << endl;

	return means;
}

void iteratedBinaryKMeans(const vector<double> &data, vector<int> &clusters,
		int N_iter) {
	vector<int> temporaryClusters(clusters.size(), 0);
	for (int i = 0; i < N_iter; ++i) {
		computeKMeans(data, temporaryClusters, 2, 1000); //1000 should be enough
		transform(temporaryClusters.cbegin(), temporaryClusters.cend(),
				temporaryClusters.begin(), [](int cluster) {
					return (cluster == 0)? 0 : -1;
				});
	}
	transform(temporaryClusters.cbegin(), temporaryClusters.cend(),
			clusters.begin(), [](int cluster) {
				return (cluster == 0)? 0 : 1;
			});
}

/*
 *
 * MULTIDIMENSIONAL K-MEANS
 *
 */

void initializeClustersRandomly(const vector<vector<double>> &data,
		vector<vector<double>> &means, int K,
		const vector<bool> &dataToCluster) {
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0, data.size() - 1);
	for (int i = 0; i < K; i++) {
		while (true) {
			int randomInt = distribution(generator);
			if (dataToCluster[randomInt]) {
				vector<double> randomData(data[randomInt]);
				if (find(means.begin(), means.begin() + i, randomData)
						== (means.begin() + i)) {
					means[i] = randomData;
					break;
				}
			}
		}
	}
}

int findClosestClusterFromDataPoint(const vector<vector<double>> &means,
		const vector<double> &dataPoint, int K, DistanceFunction distance) {
	int closestCluster = 0;
	double closestDistance = distance(dataPoint, means[0]);
	for (int i = 1; i < K; ++i) {
		double d = distance(dataPoint, means[i]);
		if (d < closestDistance) {
			closestCluster = i;
			closestDistance = d;
		}
	}

	return closestCluster;
}

void recalculateMeans(const vector<vector<double>> &data,
		vector<vector<double>> &means, const vector<int> &clusters, int K) {
	vector<int> clusterSize(K, 0);
	vector<double> nullVector(data[0].size(), 0.0);
	fill(means.begin(), means.end(), nullVector);
	for (unsigned int i = 0; i != clusters.size(); ++i) {
		int cluster = clusters[i];
		if (cluster != -1) {
			++clusterSize[cluster];
			transform(means[cluster].begin(), means[cluster].end(),
					data[i].begin(), means[cluster].begin(), plus<double>());
		}
	}
	for (int i = 0; i < K; ++i) {
		double s = (double) clusterSize[i];
		transform(means[i].begin(), means[i].end(), means[i].begin(),
				[s](double d) {return d/s;});
	}
}

bool kMeansIteration(const vector<vector<double>> &data,
		vector<vector<double>> &means, vector<int> &clusters, int K,
		DistanceFunction distance) {
	bool clustersChanged = false;

	//Assign clusters
	for (unsigned int i = 0; i != data.size(); ++i) {
		int oldCluster = clusters[i];
		if (oldCluster != -1) {
			int newCluster = findClosestClusterFromDataPoint(means, data[i], K,
					distance);
			if (newCluster != clusters[i]) {
				clustersChanged = true;
				clusters[i] = newCluster;
			}
		}
	}

	//Recalculate means
	recalculateMeans(data, means, clusters, K);

	return clustersChanged;
}

vector<vector<double>> computeKMeans(const vector<vector<double>> &data,
		vector<int> &clusters, int K, int Nmax, DistanceFunction distance) {

	//Initialize data Structures
	vector<vector<double>> means(K);
	vector<bool> dataToCluster(data.size());
	transform(clusters.cbegin(), clusters.cend(), dataToCluster.begin(),
			[](int cluster) {
				return (cluster != -1);
			});

	initializeClustersRandomly(data, means, K, dataToCluster);

	//Iterate
	int i = 0;
	while (i < Nmax && kMeansIteration(data, means, clusters, K, distance)) {
		++i;
	}

	if (i == Nmax) {
		cout << "K-Means did not converge." << endl;
	}

	return means;
}

/*
 *
 * DISTANCE MEASURES
 *
 */

double euclidianNorm(const std::vector<double> &a,
		const std::vector<double> &b) {
	double dist = 0.0;
	for_each_two_ranges(a.cbegin(), a.cend(), b.cbegin(),
			[&dist](double x, double y) {
				dist += (x-y)*(x-y);
			});
	return sqrt(dist);
}

double norm1(const std::vector<double> &a, const std::vector<double> &b) {
	double dist = 0.0;
	for_each_two_ranges(a.cbegin(), a.cend(), b.cbegin(),
			[&dist](double x, double y) {
				dist += abs(x-y);
			});
	return dist;
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
