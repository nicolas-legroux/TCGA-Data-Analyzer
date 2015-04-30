#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "k_means.hpp"

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

vector<double> computeKMeans(const vector<double> &data, vector<int> &clusters, int K, int N){

	vector<double> means(K);
	initializeClusters(data, means, K);

	int i = 0;
	while(i<N && !recomputeClusters(data, means, clusters, K)){
		++i;
	}

	cout << "K-Means finished after " << i << " iterations." << endl;

	return means;
}
