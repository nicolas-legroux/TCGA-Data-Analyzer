#ifndef STATS_HPP_INCLUDED
#define STATS_HPP_INCLUDED

#include <vector>
#include <unordered_map>
#include <string>
#include <utility>
#include <algorithm>

#include "dataReader.hpp"

double computeMean(const std::vector<double> &vec);
double computeStandardDeviation(const std::vector<double> &vec, bool correction = true);
double computeZeroPercentage(const std::vector<double> &vec);
double computePearsonCorrelation(const std::vector<double> &x, const std::vector<double> &y);
std::vector<double> computePearsonCorrelation(const std::vector<std::vector<double>> &M);
double computeSpearmanCorrelation(const std::vector<double> &x, const std::vector<double> &y);
std::vector<double> computeSpearmanCorrelation(const std::vector<std::vector<double>> &M);

std::vector<double> computePairwiseEuclideanDistance(const std::vector<std::vector<double>> &M);
std::vector<double> computePairwiseManhattanDistance(const std::vector<std::vector<double>> &M);

std::unordered_map<std::string, std::vector<std::pair<double, double>>> computeControlDistribution(const RNASeqData &controlData);
void computeZScore(RNASeqData &tumorData, const std::unordered_map<std::string, std::vector<std::pair<double, double>>> &controlDistributionParameters);

#endif // STATS_HPP_INCLUDEDs
