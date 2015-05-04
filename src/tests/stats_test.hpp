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
void spearmanCorrelationTest1(){
	std::vector<double> x{1,2,4,5,6,7,-2,9,10};
	std::vector<double> y{4,0,5,0,5,3,-5,1,0};
	std::cout << computeSpearmanCorrelation(x, y);
}


#endif /* SRC_TESTS_STATS_TEST_HPP_ */
