#ifndef STATS_HPP_INCLUDED
#define STATS_HPP_INCLUDED

#include <vector>
#include <unordered_map>
#include <string>
#include <utility>
#include <algorithm>
#include "typedefs.hpp"

#include "dataReader.hpp"

//Compute the mean of a range of real values
double computeMean(const std::vector<double> &vec);
double computeMean(const VectorX &vec);

//Compute the mean of a range of vectors
std::vector<double> computeMeanVector(
		const std::vector<std::vector<double>> &vecs);
VectorX computeMeanVector(const MatrixX &matrix);

//Compute the standard deviation of a range of real values
double computeStandardDeviation(const std::vector<double> &vec,
		bool correction = true);
double computeStandardDeviation(const VectorX &vec, bool correction = true);

//Compute the percentage of null values in a vector
double computeZeroPercentage(const std::vector<double> &vec);

//Compute the Pearson correlation between two vectors
double computePearsonCorrelation(const std::vector<double> &x,
		const std::vector<double> &y);
double computePearsonCorrelation(const VectorX &x, const VectorX &y);

double computePearsonCorrelation(const std::vector<double> &x,
		const std::vector<double> &y, double x_mean, double x_stddev,
		double y_mean, double y_stddev);
double computePearsonCorrelation(const VectorX &x, const VectorX &y,
		double x_mean, double x_stddev, double y_mean, double y_stddev);

//utility for Spearman
void computeRank(std::vector<double> &x);

//Compute the Spearman correlation between two vectors
double computeSpearmanCorrelation(const std::vector<double> &x,
		const std::vector<double> &y);
double computeSpearmanCorrelation(const VectorX &x, const VectorX &y);

std::unordered_map<std::string, std::vector<std::pair<double, double>>>computeControlDistribution(const RNASeqData &controlData);
void computeZScore(RNASeqData &tumorData,
		const std::unordered_map<std::string,
				std::vector<std::pair<double, double>>>&controlDistributionParameters);

#endif // STATS_HPP_INCLUDEDs
