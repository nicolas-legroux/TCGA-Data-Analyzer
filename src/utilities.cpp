#include <iostream>
#include <vector>
#include <algorithm>
#include "utilities.hpp"

void printAdvancement(unsigned int currentCount, unsigned int totalCount) {
	std::cout << (100.0 * (double) currentCount) / (double) (totalCount) << "% \r"
			<< std::flush;
}

std::vector<std::string> split(const std::string &s, const std::vector<char> &delimiters) {
	std::vector<std::string> strs;
	std::string currentString;

	for (const char &c : s) {
		if (find(delimiters.cbegin(), delimiters.cend(), c)
				== delimiters.cend()) {
			currentString.push_back(c);
		} else {
			strs.push_back(currentString);
			currentString.erase();
		}
	}

	strs.push_back(currentString);

	return strs;
}

double computeMean(const std::vector<double> &vec) {
	double sum = accumulate(vec.cbegin(), vec.cend(), 0.0);
	return sum / (double) vec.size();
}

double computeStandardDeviation(const std::vector<double> &vec, bool correction) {
	double m = computeMean(vec);
	double accum = 0.0;
	for_each(vec.cbegin(), vec.cend(), [&](const double d) {
		accum += (d-m)*(d-m);
	});
	if (correction) {
		return sqrt(accum / (vec.size() - 1));
	} else {
		return sqrt(accum / vec.size());
	}
}


