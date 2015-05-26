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
*/
