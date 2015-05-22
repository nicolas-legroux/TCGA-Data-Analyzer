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
#include <assert.h>

//Iterate on two ranges and call a binary-op function
template<typename InputIter1, typename InputIter2, typename Function>
Function for_each_two_ranges(InputIter1 first1, InputIter1 last1, InputIter2 first2, Function f) {
    for (; first1 != last1; ++first1, ++first2) {
        f(*first1, *first2);
    }
    return f;
}

template<typename InputIter1, typename InputIter2, typename OutputIter>
//The second range should correspond to a container of boolean values saying whether the value should be added to the output
void copy_if_two_ranges(InputIter1 first1, InputIter1 last1, InputIter2 first2, OutputIter out) {
    for (; first1 != last1; ++first1, ++first2) {
        if(*first2){
        	*out = *first1;
        }
    }
}

//Print a vector
template <typename T>
void print_vector(const std::vector<T> &v) {
	std::cout << "{ ";
	for_each(v.cbegin(), v.cend(), [](const T &data){ std::cout << data << " ";});
	std::cout << "}" << std::endl;
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

//Prints advancement of a task in %
void printAdvancement(unsigned int currentCount, unsigned int totalCount);

//Splits a string according to delimiters
std::vector<std::string> split(const std::string &s, const std::vector<char> &delimiters);

//Returns the number of pairs in a set of cardinal n
int numberOfPairs(int n);

/*
 *
 * DISTANCE MEASURES
 *
 */

double euclideanDistance(const std::vector<double> &a,
		const std::vector<double> &b);
double manhattanDistance(const std::vector<double> &a, const std::vector<double> &b);

#endif /* SRC_UTILITIES_HPP_ */
