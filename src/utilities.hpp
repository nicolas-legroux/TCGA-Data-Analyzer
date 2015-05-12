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
#include <iostream>

//Iterate on two ranges and call a binary-op function
template<typename InputIter1, typename InputIter2, typename Function>
Function for_each_two_ranges(InputIter1 first1, InputIter1 last1, InputIter2 first2, Function f) {
    for (; first1 != last1; ++first1, ++first2) {
        f(*first1, *first2);
    }
    return f;
}

//Sort the indexes of a vector in decreasing order
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

//Sort the indexes of an array in increasing order
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

//Get the rank of each value of an array, where the rank is taken in increasing order
template <typename T>
std::vector<size_t> get_rank_increasing(const std::vector<T> &v){
	std::vector<size_t> sortedIdx{sort_indexes_increasing(v)};
	std::vector<size_t> ranks(sortedIdx.size());
	for(size_t i = 0; i != ranks.size(); ++i){
		ranks[sortedIdx[i]] = i;
	}
	return ranks;
}

void printAdvancement(unsigned int currentCount, unsigned int totalCount);

#endif /* SRC_UTILITIES_HPP_ */
