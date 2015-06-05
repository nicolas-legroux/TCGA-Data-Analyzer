// Copyright Nicolas Legroux 2015

#include <iostream>
#include <vector>
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

#include "distanceMatrix.hpp"
#include "spectral_clustering.hpp"
#include "dataReader.hpp"

using std::vector;
using std::string;
using std::map;

int main() {

	vector<string> cancers = { "LUSC", "KIRC", "BRCA", "THCA"};
	int maxControl = 50;
	int maxTumor = 200;

	UnsupervisedNormalizationMethod method =
			UnsupervisedNormalizationMethod::BINARY_QUANTILE;
	UnsupervisedNormalizationParameters parameters;
	parameters.setBinaryQuantileParameters(0.45);
	parameters.setKMeansParameters(2, 1000);
	parameters.setBinaryIteratedKMeansParameters(6);

	clustering_KMeans_test(cancers, maxControl, maxTumor, method, parameters);

//	clustering_Hierarchical_test(cancers, maxControl, maxTumor, method,
//			parameters, DistanceMetric::PEARSON_CORRELATION,
//			LinkageMethod::COMPLETE);

//	clustering_Spectral_test(cancers, maxControl, maxTumor, method,
//		parameters, DistanceMetric::PEARSON_CORRELATION);



	return 0;
}
