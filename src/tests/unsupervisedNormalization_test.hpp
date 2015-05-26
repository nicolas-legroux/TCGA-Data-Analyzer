#ifndef SRC_TESTS_UNSUPERVISEDNORMALIZATION_TEST_HPP_
#define SRC_TESTS_UNSUPERVISEDNORMALIZATION_TEST_HPP_

#include "../unsupervisedNormalization.hpp"
#include "../distanceMatrix.hpp"

void unsupervisedNormalization_test(
		const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters,
		const DistanceMetric &distanceMetric = DistanceMetric::PEARSON_CORRELATION);

#endif /* SRC_TESTS_UNSUPERVISEDNORMALIZATION_TEST_HPP_ */
