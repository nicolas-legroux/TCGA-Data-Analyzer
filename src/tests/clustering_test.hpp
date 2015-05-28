#ifndef SRC_TESTS_CLUSTERING_TEST_HPP_
#define SRC_TESTS_CLUSTERING_TEST_HPP_

#include "../unsupervisedNormalization.hpp"
#include "../distanceMatrix.hpp"
#include "../hierarchical_clustering.hpp"

void clustering_KMeans_test();

void clustering_KMeans_test(const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters);

void clustering_Hierarchical_test(const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters,
		const DistanceMetric &distanceMetric,
		const LinkageMethod &linkageMethod);

#endif /* SRC_TESTS_CLUSTERING_TEST_HPP_ */
