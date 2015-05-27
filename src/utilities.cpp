#include <iostream>
#include <vector>
#include <string>

#include "utilities.hpp"

using namespace std;

void printAdvancement(unsigned int currentCount, unsigned int totalCount) {
	cout << (100 * currentCount) / (totalCount) << "% \r" << flush;
}

vector<string> split(const string &s, const vector<char> &delimiters) {
	vector<string> strs;
	string currentString;

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

int numberOfPairs(int n) {
	return n * (n - 1) / 2;
}

/*
 *
 * DISTANCE MEASURES
 *
 */

double euclideanDistance(const std::vector<double> &a,
		const std::vector<double> &b) {
	double dist = 0.0;
	for_each_two_ranges(a.cbegin(), a.cend(), b.cbegin(),
			[&dist](double x, double y) {
				dist += (x-y)*(x-y);
			});
	return sqrt(dist);
}

double manhattanDistance(const std::vector<double> &a,
		const std::vector<double> &b) {
	double dist = 0.0;
	for_each_two_ranges(a.cbegin(), a.cend(), b.cbegin(),
			[&dist](double x, double y) {
				dist += abs(x-y);
			});
	return dist;
}

double cosineSimilarity(const std::vector<double> &a,
		const std::vector<double> &b) {
	double product = 0.0;
	double normA = 0.0;
	double normB = 0.0;
	for_each_two_ranges(a.cbegin(), a.cend(), b.cbegin(),
			[&product,&normA,&normB](double ai, double bi) {
				product += ai*bi;
				normA += ai*ai;
				normB += bi*bi;
			});
	return product / (sqrt(normA)*sqrt(normB));
}

/*
 *
 * UTILITIES FOR GENERIC K MEANS
 *
 */

double distanceDouble(const double &x, const double &y) {
	return abs(x - y);
}

void addToDouble(double &d, const double &x) {
	d += x;
}

void divideDoubleByConstant(double &d, double c) {
	d /= c;
}

void addToVector(std::vector<double> &v, const std::vector<double> &x) {
	assert(x.size() == v.size());
	for_each_two_ranges(v.begin(), v.end(), x.cbegin(),
			[](double &vd, double xd) {
				vd += xd;
			});
}

void divideVectorByConstant(std::vector<double> &v, double c) {
	for_each(v.begin(), v.end(), [c](double &vd) {
		vd /= c;
	});
}
