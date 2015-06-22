#ifndef SRC_UTILITIES_HPP_
#define SRC_UTILITIES_HPP_

#include <vector>

//Prints advancement of a task in %
void printAdvancement(unsigned int currentCount, unsigned int totalCount);

//Splits a string according to delimiters
std::vector<std::string> split(const std::string &s,
		const std::vector<char> &delimiters);

double computeMean(const std::vector<double> &vec);
double computeStandardDeviation(const std::vector<double> &vec, bool correction = true);

#endif /* SRC_UTILITIES_HPP_ */
