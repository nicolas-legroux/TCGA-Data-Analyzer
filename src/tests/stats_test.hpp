/*
 * stats_test.hpp
 *
 *  Created on: May 4, 2015
 *      Author: nicolas
 */

#ifndef SRC_TESTS_STATS_TEST_HPP_
#define SRC_TESTS_STATS_TEST_HPP_

#include <iostream>
#include <vector>
#include "../stats.hpp"

//Test with ties
//R output is 0.136201
//Output is 0.136201
void spearmanCorrelationTest1() {
	std::vector<double> x { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	std::vector<double> y { 4, 0, 5, 0, 5, 3, -5, 1, 0 };
	std::cout << computeSpearmanCorrelation(x, y) << std::endl;
}

//Same output as R :
/*
 M = matrix(c(1,2,4,5,6,7,-2,9,10,4,0,5,0,5,3,-5,1,0,0,0,-5,9,-4,1,1,2,1), nrow=9, ncol=3)
 cor(M, method="spearman")
 */
void spearmanCorrelationTest2() {
	std::vector<double> x1 { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	std::vector<double> x2 { 4, 0, 5, 0, 5, 3, -5, 1, 0 };
	std::vector<double> x3 { 0, 0, -5, 9, -4, 1, 1, 2, 1 };
	std::vector<std::vector<double>> M { x1, x2, x3 };
	std::vector<double> correlationMatrix = computeSpearmanCorrelation(M);

	for (int l = 0; l < 3; ++l) {
		for (int c = 0; c < 3; ++c) {
			std::cout << correlationMatrix[3 * l + c] << "\t";
		}
		std::cout << std::endl;
	}
}

#endif /* SRC_TESTS_STATS_TEST_HPP_ */
