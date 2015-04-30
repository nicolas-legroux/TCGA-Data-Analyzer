#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

#include "k_means.hpp"
#include "stats.hpp"

using namespace std;

void initializeClusters(const vector<double> &data, vector<double> &means, int K) {
	vector<double> copyData(data);
	auto it = unique(copyData.begin(), copyData.end());
	int N_unique = it-copyData.begin();
	int step = N_unique/(K+1);
	for(int i=0; i<K; ++i){
		means[i] = copyData[(i+1)*step];
	}
}

void initializeClustersRandomly(const vector<double> &data, vector<double> &means, int K){
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0,data.size()-1);
	for(int i=0; i<K; i++){
		while(true){
			int randomInt = distribution(generator);
			double randomData = data[randomInt];
			if(find(means.begin(), means.begin()+i, randomData) == (means.begin()+i)){
				cout << "OK" << endl;
				means[i] = randomData;
				break;
			}
		}
	}
}

int findClosestCluster(const vector<double> &means, double dataPoint, int K){
	int closestCluster = 0;
	double closestDistance = abs(dataPoint - means[0]);
	for(int i=1; i<K; ++i){
		double d = abs(dataPoint-means[i]);
		if(d < closestDistance){
			closestCluster = i;
			closestDistance = d;
		}
	}

	return closestCluster;
}

void computeNewMeans(const vector<double> &data, vector<double> &means, const vector<int> &clusters, int K){
	vector<int> clusterSize(K, 0);
	fill(means.begin(), means.end(), 0.0);
	for(unsigned int i=0; i != clusters.size(); ++i){
		int cluster = clusters[i];
		++clusterSize[cluster];
		means[cluster] += data[i];
	}
	for(int i=0; i<K; ++i){
		means[i] /= (double)clusterSize[i];
	}
}

bool recomputeClusters(const vector<double> &data, vector<double> &means, vector<int> &clusters, int K){
	bool clustersChanged = false;

	//Assign clusters
	for(unsigned int i=0; i != data.size(); ++i){
		int newCluster = findClosestCluster(means, data[i], K);
		if(newCluster != clusters[i]){
			clustersChanged = true;
			clusters[i] = newCluster;
		}
	}

	//Compute new means
	computeNewMeans(data, means, clusters, K);

	return clustersChanged;
}

void sortClusters(vector<int> &clusters, const vector<int> &mapping){
	for(unsigned int i=0; i != clusters.size(); ++i){
		int currentCluster = clusters[i];
		clusters[i] = mapping.at(currentCluster);
	}
}

vector<double> computeKMeans(const vector<double> &data, vector<int> &clusters, int K, int N){

	vector<double> means(K);
	initializeClustersRandomly(data, means, K);

	int i = 0;
	while(i<N && recomputeClusters(data, means, clusters, K)){
		++i;
	}

	vector<size_t> sortedClusterIndexes(sort_indexes_increasing<double>(means));

	vector<int> sortingMapping(K);
	for(int i=0; i!=K; ++i){
		sortingMapping[sortedClusterIndexes[i]] = i;
	}

	sortClusters(clusters, sortingMapping);
	sort(means.begin(), means.end());

	cout << "K-Means finished after " << i << " iterations." << endl;

	return means;
}
