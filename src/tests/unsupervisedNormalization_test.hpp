#ifndef SRC_TESTS_UNSUPERVISEDNORMALIZATION_TEST_HPP_
#define SRC_TESTS_UNSUPERVISEDNORMALIZATION_TEST_HPP_

void normalizationTest1KMeans(int K, int Nmax);
void normalizationTestIteratedKMeans(int Niter);
void normalizationTestQuantile(double cutPercentage);
void normalization_KMeans_Manhattan_test(int K, int Nmax);
void no_normalization_euclidean_test();

#endif /* SRC_TESTS_UNSUPERVISEDNORMALIZATION_TEST_HPP_ */
