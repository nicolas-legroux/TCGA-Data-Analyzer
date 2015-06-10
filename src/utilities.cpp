#include <iostream>
#include <vector>
#include <string>

#include "utilities.hpp"

using namespace std;

void printAdvancement(unsigned int currentCount, unsigned int totalCount) {
	cout << (100.0 * (double) currentCount) / (double) (totalCount) << "% \r"
			<< flush;
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

double euclideanNorm(const std::vector<double> &x) {
	double norm = 0.0;
	for_each(x.cbegin(), x.cend(), [&](double xi) {
		norm += xi*xi;
	});
	return sqrt(norm);
}

double euclideanDistance(const std::vector<double> &a,
		const std::vector<double> &b) {
	double dist = 0.0;
	for_each_two_ranges(a.cbegin(), a.cend(), b.cbegin(),
			[&dist](double x, double y) {
				dist += (x-y)*(x-y);
			});
	return sqrt(dist);
}

double euclideanDistance(const VectorX &a, const VectorX &b) {
	return (a - b).norm();
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

double manhattanDistance(const VectorX &a, const VectorX &b) {
	return (a - b).lpNorm<1>();
}

double cosineSimilarity(const std::vector<double> &a,
		const std::vector<double> &b, double normA, double normB) {
	double product = 0.0;
	for_each_two_ranges(a.cbegin(), a.cend(), b.cbegin(),
			[&product](double ai, double bi) {
				product += ai*bi;});
	return product / (normA * normB);
}

double cosineSimilarity(const VectorX &a, const VectorX &b, double normA,
		double normB) {
	return a.dot(b) / (normA * normB);
}

double cosineSimilarity(const std::vector<double> &a,
		const std::vector<double> &b) {
	double product = 0.0;
	double normA = 0.0;
	double normB = 0.0;
	for_each_two_ranges(a.cbegin(), a.cend(), b.cbegin(),
			[&product, &normA, &normB](double ai, double bi) {
				product += ai*bi;
				normA += ai*ai;
				normB += bi*bi;
			});
	return product / (sqrt(normA) * sqrt(normB));
}

double cosineSimilarity(const VectorX &a, const VectorX &b) {
	return a.dot(b) / (a.norm() * b.norm());
}

