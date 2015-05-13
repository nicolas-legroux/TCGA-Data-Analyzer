#include <vector>
#include <string>
#include <iostream>

#include "../utilities.hpp"
#include "utilities_test.hpp"

void splitTest(){
	std::string test = "TCGA-A6-5661-01.genes.normalized.results";
	std::vector<std::string> strs = split(test, std::vector<char>{'.', '-'});
	for(const std::string &s : strs){
		std::cout << s << std::endl;
	}
}
