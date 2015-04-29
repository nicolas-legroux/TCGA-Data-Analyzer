#ifndef STATS_HPP_INCLUDED
#define STATS_HPP_INCLUDED

#include <vector>
#include <unordered_map>
#include <string>
#include <utility>
#include <algorithm>

#include "dataReader.hpp"

double computeMean(const std::vector<double> &vec);
double computeStandardDeviation(const std::vector<double> &vec);
double computeZeroPercentage(const std::vector<double> &vec);
double computeCorrelation(const std::vector<double> &x, const std::vector<double> &y);

std::unordered_map<std::string, std::vector<std::pair<double, double>>> computeControlDistribution(const RNASeqData &controlData);
void computeZScore(RNASeqData &tumorData, const std::unordered_map<std::string, std::vector<std::pair<double, double>>> &controlDistributionParameters);

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

#endif // STATS_HPP_INCLUDEDs
