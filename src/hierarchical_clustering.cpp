#include "hierarchical_clustering.hpp"
#include <limits>
#include <utility>
#include <cmath>
#include <cassert>
#include <map>
#include "utilities.hpp"

using namespace std;

Hierarchical_Clustering::Hierarchical_Clustering(const vector<double> &matrix,
		LinkageMethod _linkageMethod, MatrixType _matrixType) :
		linkageMethod(_linkageMethod), matrixType(_matrixType) {
	n = sqrt(matrix.size());

	unionFindDataStructure.resize(n);
	fill(unionFindDataStructure.begin(), unionFindDataStructure.end(), -1);
	for (unsigned int i = 0; i < n; ++i) {
		clusterRepresentatives.insert(i);
	}

	if (linkageMethod == LinkageMethod::AVERAGE) {
		clusterSizes.resize(n);
		fill(clusterSizes.begin(), clusterSizes.end(), 1);
	}

	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j <= i; ++j) {
			data.push_back(matrix[i * n + j]);
		}
	}

	assert(data.size() == n * (n + 1) / 2);
}

double& Hierarchical_Clustering::getDistance(int i, int j) {
	if (j > i) {
		swap(i, j);
	}
	return data[i * (i + 1) / 2 + j];
}

double Hierarchical_Clustering::worstPossibleDistance() {
	switch (matrixType) {
	case MatrixType::SIMILARITY:
		return -1.0*std::numeric_limits<double>::infinity();
	case MatrixType::DISTANCE:
		return std::numeric_limits<double>::infinity();
	default:
		throw invalid_argument("Unknown matrix type");
	}
}

bool Hierarchical_Clustering::isBetterDistance(double oldDistance,
		double newDistance) {
	switch (matrixType) {
	case MatrixType::SIMILARITY:
		return newDistance > oldDistance;
	case MatrixType::DISTANCE:
		return newDistance < oldDistance;
	default:
		throw invalid_argument("Unknown matrix type");
	}
}

pair<int, int> Hierarchical_Clustering::findClustersToMerge() {
	pair<int, int> clustersToMerge { -1, -1 };
	double bestDistance = worstPossibleDistance();
	for (int i : clusterRepresentatives) {
		for (int j : clusterRepresentatives) {
			if(i != j){
				double d = getDistance(i, j);
				if (isBetterDistance(bestDistance, d)) {
					bestDistance = d;
					clustersToMerge = make_pair(i, j);
				}
			}
		}
	}

	return clustersToMerge;
}

void Hierarchical_Clustering::updateDistances(int deletedCluster,
		int mergedCluster) {
	for (int i : clusterRepresentatives) {
		if (i != mergedCluster) {
			double dist1 = getDistance(i, deletedCluster);
			double dist2 = getDistance(i, mergedCluster);
			if (linkageMethod == LinkageMethod::COMPLETE) {
				switch (matrixType) {
				case MatrixType::DISTANCE:
					getDistance(i, mergedCluster) = max(dist1, dist2);
					break;
				case MatrixType::SIMILARITY:
					getDistance(i, mergedCluster) = min(dist1, dist2);
					break;
				default:
					throw invalid_argument("Unknown matrix type");
				}
			} else if (linkageMethod == LinkageMethod::SINGLE) {
				switch (matrixType) {
				case MatrixType::DISTANCE:
					getDistance(i, mergedCluster) = min(dist1, dist2);
					break;
				case MatrixType::SIMILARITY:
					getDistance(i, mergedCluster) = max(dist1, dist2);
					break;
				default:
					throw invalid_argument("Unknown matrix type");
				}
			} else if (linkageMethod == LinkageMethod::AVERAGE) {
				double deletedClusterSize = clusterSizes[deletedCluster];
				double mergedClusterSize = clusterSizes[mergedCluster];
				getDistance(i, mergedCluster) = (deletedClusterSize * dist1
						+ mergedClusterSize * dist2)
						/ (deletedClusterSize + mergedClusterSize);
			} else {
				throw invalid_argument("Unknown linkage method");
			}
		}
	}
}

void Hierarchical_Clustering::mergeClusters(int i, int j) {
	clusterRepresentatives.erase(i);
	unionFindDataStructure[i] = j;
	updateDistances(i, j);
	if (linkageMethod == LinkageMethod::AVERAGE) {
		clusterSizes[j] += clusterSizes[i];
	}
}

int Hierarchical_Clustering::findClusterRepresentative(int i) {
	int next = unionFindDataStructure[i];
	while (next != -1) {
		i = next;
		next = unionFindDataStructure[next];
	}
	return i;
}

vector<int> Hierarchical_Clustering::compute(unsigned int k) {
	assert(k <= clusterRepresentatives.size() && k >= 2);
	while (clusterRepresentatives.size() > k) {
		pair<int, int> clusterToMerge = findClustersToMerge();
		int i = clusterToMerge.first;
		int j = clusterToMerge.second;
		assert(i >= 0 && j >= 0);
		//cout << "Merging " << i << " and " << j << endl;
		mergeClusters(i, j);
	}

	vector<int> clusters(n);
	map<int, unsigned int> map = buildIndexMap<int>(
			clusterRepresentatives.cbegin(), clusterRepresentatives.cend());
	for (unsigned int i = 0; i < n; ++i) {
		int j = findClusterRepresentative(i);
		clusters[i] = map[j];
	}
	return clusters;
}
