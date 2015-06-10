#include <vector>
#include <cmath>
#include <unordered_map>
#include <string>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cassert>

#include "stats.hpp"
#include "typedefs.hpp"
#include "utilities.hpp"

using namespace std;

//Let's assume it won't overflow
double computeMean(const vector<double> &vec) {
	double sum = accumulate(vec.cbegin(), vec.cend(), 0.0);
	return sum / (double) vec.size();
}

double computeMean(const VectorX &vec) {
	double sum = vec.sum();
	return sum / (double) vec.rows();
}

vector<double> computeMeanVector(const vector<vector<double>> &vecs) {
	unsigned int s = vecs.size();
	assert(s > 0);
	unsigned int dim = (*vecs.cbegin()).size();
	vector<double> mean(dim, 0.0);
	for_each(vecs.cbegin(), vecs.cend(),
			[dim, &mean](const vector<double> &vec) {
				assert(dim == vec.size());
				transform(mean.begin(), mean.end(), vec.cbegin(), mean.begin(), plus<double>());
			});

	transform(mean.begin(), mean.end(), mean.begin(),
			[s](double d) {return d/s;});
	return mean;
}

VectorX computeMeanVector(const MatrixX &matrix) {
	return matrix.colwise().mean();
}

double computeStandardDeviation(const vector<double> &vec, bool correction) {
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

double computeStandardDeviation(const VectorX &vec, bool correction) {
	double m = computeMean(vec);
	double accum = 0;
	for (unsigned int i = 0; i < vec.rows(); ++i) {
		double d = vec(i);
		accum += (d - m) * (d - m);
	}

	if (correction) {
		return std::sqrt(accum / (vec.rows() - 1));
	} else {
		return std::sqrt(accum / vec.rows());
	}
}

double computeZeroPercentage(const vector<double> &vec) {
	return (double) count(vec.cbegin(), vec.cend(), 0.0) / (double) vec.size();
}

double computePearsonCorrelation(const vector<double> &x,
		const vector<double> &y) {
	double accum = 0.0;
	double x_mean = computeMean(x);
	double y_mean = computeMean(y);
	for_each_two_ranges(x.cbegin(), x.cend(), y.cbegin(),
			[&accum, x_mean, y_mean](const double a, const double b) {
				accum += (a-x_mean)*(b-y_mean);
			});
	double x_stddev = computeStandardDeviation(x, false);
	double y_stddev = computeStandardDeviation(y, false);
	return accum / (x.size() * x_stddev * y_stddev);
}

double computePearsonCorrelation(const VectorX &x, const VectorX &y) {
	double accum = 0.0;
	double x_mean = computeMean(x);
	double y_mean = computeMean(y);
	assert(x.rows() == y.rows());
	for (unsigned int i = 0; i < x.rows(); ++i) {
		accum += (x(i) - x_mean) * (y(i) - y_mean);
	}
	double x_stddev = computeStandardDeviation(x, false);
	double y_stddev = computeStandardDeviation(y, false);
	return accum / (x.size() * x_stddev * y_stddev);
}

double computePearsonCorrelation(const vector<double> &x,
		const vector<double> &y, double x_mean, double x_stddev, double y_mean,
		double y_stddev) {
	double accum = 0.0;
	for_each_two_ranges(x.cbegin(), x.cend(), y.cbegin(),
			[&accum, x_mean, y_mean](const double a, const double b) {
				accum += (a-x_mean)*(b-y_mean);
			});
	return accum / (x.size() * x_stddev * y_stddev);
}

double computePearsonCorrelation(const VectorX &x, const VectorX &y,
		double x_mean, double x_stddev, double y_mean, double y_stddev) {
	double accum = 0.0;
	assert(x.rows() == y.rows());
	for (unsigned int i = 0; i < x.rows(); ++i) {
		accum += (x(i) - x_mean) * (y(i) - y_mean);
	}
	return accum / (x.size() * x_stddev * y_stddev);
}

//Assume x is sorted. Utility for Spearman correlation
void computeRankSorted(vector<double> &x) {
	unsigned int current = 0;
	while (current != x.size()) {
		double d = x[current];
		unsigned int next = current + 1;
		int sum = current;
		int count = 1;
		while (next < x.size() && d == x[next]) {
			++count;
			sum += next;
			++next;
		}

		for (unsigned int j = current; j < next; ++j) {
			x[j] = (double) sum / (double) count;
		}

		current = next;
	}
}

void computeRank(vector<double> &x) {
	vector<double> copyX(x);
	vector<size_t> sortedXIndexes = get_rank_increasing(copyX);
	sort(copyX.begin(), copyX.end());
	computeRankSorted(copyX);
	for (unsigned int i = 0; i != x.size(); ++i) {
		x[i] = copyX[sortedXIndexes[i]];
	}
}

double computeSpearmanCorrelation(const vector<double> &x,
		const vector<double> &y) {
	vector<double> copyX(x);
	vector<double> copyY(y);
	computeRank(copyX);
	computeRank(copyY);
	return computePearsonCorrelation(copyX, copyY);
}

double computeSpearmanCorrelation(const VectorX &x,
		const VectorX &y) {
	vector<double> copyX;
	vector<double> copyY;
	assert(x.rows() == y.rows());
	for (unsigned int i = 0; i != x.rows(); ++i) {
		double xi = x(i);
		double yi = y(i);
		copyX.push_back(xi);
		copyY.push_back(yi);
	}
	computeRank(copyX);
	computeRank(copyY);
	return computePearsonCorrelation(copyX, copyY);
}

unordered_map<string, vector<pair<double, double>>> computeControlDistribution(const RNASeqData &controlData) {

	unordered_map<string, vector<pair<double, double>>> controlDistributionParameters;

	cout << endl << "*********** Computing mean and standard deviation for each gene ************" << endl;

	for(const auto &mappedData : controlData) {
		string cancerName = mappedData.first;
		if(controlData.at(cancerName).at(0).empty()) {
			continue;
		}
		controlDistributionParameters.insert(make_pair(cancerName, vector<pair<double, double>>()));
		for(const vector<double> &vec : mappedData.second) {
			double mean = computeMean(vec);
			double stddev = computeStandardDeviation(vec);
			controlDistributionParameters[cancerName].push_back(make_pair(mean, stddev));
		}
	}

	return controlDistributionParameters;
}

void computeZScore(RNASeqData &tumorData,
		const unordered_map<string, vector<pair<double, double>>>&controlDistributionParameters) {

			cout << endl << "********* Computing Z Scores ********** " << endl;

			for(auto &mappedData : tumorData) {
				string cancerName = mappedData.first;
				if(controlDistributionParameters.find(cancerName) == controlDistributionParameters.end()) {
					cout << "Problem with Z Score for " << cancerName << endl;
					return;
				}

				int i=0;
				for(vector<double> &vec : mappedData.second) {
					for(double &d : vec) {
						double stddev = controlDistributionParameters.at(cancerName).at(i).second;
						if(stddev>0) {
							d = (d-controlDistributionParameters.at(cancerName).at(i).first)/stddev;
						}
						else {
							d = 0;
						}
					}
					++i;
				}
			}
		}
