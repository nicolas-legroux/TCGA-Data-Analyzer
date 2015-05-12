/*
 * utilities_test.hpp
 *
 *  Created on: May 12, 2015
 *      Author: nicolas
 */

#ifndef SRC_TESTS_UTILITIES_TEST_HPP_
#define SRC_TESTS_UTILITIES_TEST_HPP_

#include <vector>
#include <string>
#include <iostream>

#include "../utilities.hpp"

void splitTest(){
	std::string test = "TCGA-A6-5661-01.genes.normalized.results";
	std::vector<std::string> strs = split(test, std::vector<char>{'.', '-'});
	for(const std::string &s : strs){
		std::cout << s << std::endl;
	}
}

#endif /* SRC_TESTS_UTILITIES_TEST_HPP_ */
