#include <vector>
#include "clustering.hpp"
#include "utilities.hpp"
#include "typedefs.hpp"
#include "k_means.hpp"

using namespace std;

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

vector<int> clusterWithKMeans(const vector<vector<double>> &data, int K, int Nmax){
	vector<int> clusters(data.size(), 0);
	computeKMeans(data, clusters, K, Nmax, euclideanDistance);
	return clusters;
}

