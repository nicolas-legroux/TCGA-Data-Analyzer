/*
 * utilities.hpp
 *
 *  Created on: May 11, 2015
 *      Author: nicolas
 */

#ifndef SRC_UTILITIES_HPP_
#define SRC_UTILITIES_HPP_

#include <vector>
#include <algorithm>

template <typename T>
std::vector<size_t> sort_indexes_decreasing(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

template <typename T>
std::vector<size_t> sort_indexes_increasing(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template <typename T>
std::vector<size_t> get_rank_increasing(const std::vector<T> &v){
	std::vector<size_t> sortedIdx{sort_indexes_increasing(v)};
	std::vector<size_t> ranks(sortedIdx.size());
	for(size_t i = 0; i != ranks.size(); ++i){
		ranks[sortedIdx[i]] = i;
	}
	return ranks;
}



#endif /* SRC_UTILITIES_HPP_ */
