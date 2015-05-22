#ifndef STATS_HPP_INCLUDED
#define STATS_HPP_INCLUDED

#include <vector>
#include <unordered_map>
#include <string>
#include <utility>
#include <algorithm>

#include "dataReader.hpp"

//Compute the mean of a range of real values
double computeMean(const std::vector<double> &vec);

//Compute the mean of a range of vectors
std::vector<double> computeMeanVector(const std::vector<std::vector<double>> &vecs);

//Compute the standard deviation of a range of real values
double computeStandardDeviation(const std::vector<double> &vec, bool correction = true);

//Compute the percentage of null values in a vector
double computeZeroPercentage(const std::vector<double> &vec);

//Compute the pearson correlation between two vectors
double computePearsonCorrelation(const std::vector<double> &x, const std::vector<double> &y);

//Compute the pearson correlation matrix (output in 1D) between a range of vectors
std::vector<double> computePearsonCorrelation(const std::vector<std::vector<double>> &M);

//Compute the spearman correlation between two vectors
double computeSpearmanCorrelation(const std::vector<double> &x, const std::vector<double> &y);

//Compute the spearman correlation matrix (output in 1D) between a range of vectors
std::vector<double> computeSpearmanCorrelation(const std::vector<std::vector<double>> &M);

//Compute the euclidean distance matrix between a range of vectors
std::vector<double> computePairwiseEuclideanDistance(const std::vector<std::vector<double>> &M);

//Compute the manhattan distance matrix between a range of vectors
std::vector<double> computePairwiseManhattanDistance(const std::vector<std::vector<double>> &M);

std::unordered_map<std::string, std::vector<std::pair<double, double>>> computeControlDistribution(const RNASeqData &controlData);
void computeZScore(RNASeqData &tumorData, const std::unordered_map<std::string, std::vector<std::pair<double, double>>> &controlDistributionParameters);

#endif // STATS_HPP_INCLUDEDs
