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
#include "normedSpace.hpp"
#include "utilities.hpp"

int main() {
	meanVectorTest();
	return 0;
}
