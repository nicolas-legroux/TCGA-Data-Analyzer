#ifndef SRC_UTILITIES_HPP_
#define SRC_UTILITIES_HPP_

#include <vector>

//Prints advancement of a task in %
void printAdvancement(unsigned int currentCount, unsigned int totalCount);

//Splits a string according to delimiters
std::vector<std::string> split(const std::string &s,
		const std::vector<char> &delimiters);

template<typename Iter>
std::string implode(Iter begin, Iter end, const std::string &delimiter){
	std::string result;
	for(Iter i = begin; i != end; ++i){
		if(i != begin){
			result += delimiter;
		}
		result += *i;
	}
	return result;
}

std::string removeTrailingZeros(std::string s);

double computeMean(const std::vector<double> &vec);
double computeStandardDeviation(const std::vector<double> &vec, bool correction = true);

#endif /* SRC_UTILITIES_HPP_ */
