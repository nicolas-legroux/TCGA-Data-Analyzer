// Copyright Nicolas Legroux 2015

#include <limits>
#include <iostream>
#include <vector>
#include <climits>
#include "tests/stats_test.hpp"
#include "tests/k_means_test.hpp"
#include "tests/lodePNG_test.hpp"
#include "tests/dataReader_test.hpp"
#include "tests/unsupervisedNormalization_test.hpp"
#include "tests/utilities_test.hpp"
#include "tests/clustering_test.hpp"
#include "tests/hierarchicalClustering_test.hpp"
#include "unsupervisedNormalization.hpp"
#include "distanceMatrix.hpp"
#include "normedVectorSpace.hpp"
#include "utilities.hpp"

using std::vector;
using std::string;

int main() {

	vector<string> cancers = { "BRCA" };
	int maxControl = 50;
	int maxTumor = 200;

	UnsupervisedNormalizationMethod method =
			UnsupervisedNormalizationMethod::BINARY_QUANTILE;
	UnsupervisedNormalizationParameters parameters;
	parameters.setBinaryQuantileParameters(0.4);
	parameters.setKMeansParameters(2, 1000);

	// clustering_KMeans_test(cancers, maxControl, maxTumor, method, parameters);

	clustering_Hierarchical_test(cancers, maxControl, maxTumor, method,
			parameters, DistanceMetric::EUCLIDEAN_DISTANCE,
			LinkageMethod::COMPLETE);

	return 0;
}
