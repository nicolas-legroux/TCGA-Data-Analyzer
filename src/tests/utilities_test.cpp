#include <vector>
#include <string>
#include <iostream>
#include <iterator>

#include "../utilities.hpp"
#include "utilities_test.hpp"

void splitTest(){
	std::string test = "TCGA-A6-5661-01.genes.normalized.results";
	std::vector<std::string> strs = split(test, std::vector<char>{'.', '-'});
	for(const std::string &s : strs){
		std::cout << s << std::endl;
	}
}

void copy_if_two_ranges_test(){
	std::vector<double> data{1,2,3,4,5,6,7,8,9};
	std::vector<bool> mustCopy{true, true, false, false, true, true, true, false, false};
	std::vector<double> dataCopy;
	copy_if_two_ranges(data.cbegin(), data.cend(), mustCopy.cbegin(), back_inserter(dataCopy));
	for_each(dataCopy.cbegin(), dataCopy.cend(), [](double d){ std::cout << d << " "; });
}
