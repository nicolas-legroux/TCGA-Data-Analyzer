#ifndef SRC_TESTS_K_MEANS_TEST_HPP_
#define SRC_TESTS_K_MEANS_TEST_HPP_

#include <string>

void kMeansTest1(int K, int Nmax, std::string cancerName, int patientId);
void iteratedBinaryKMeans_test(int N_iter, std::string cancerName,
		int patientId);

#endif /* SRC_TESTS_K_MEANS_TEST_HPP_ */
