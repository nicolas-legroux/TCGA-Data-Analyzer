#ifndef SRC_NORMEDVECTORSPACE_HPP_
#define SRC_NORMEDVECTORSPACE_HPP_

#include <iostream>
#include <cmath>
#include <cassert>
#include "utilities.hpp"
#include <iostream>

template<typename T, typename K = double>
class NormedVectorSpace {
public:
	unsigned int getDimension() const;
	K distance(const T &, const T &) const;
	T zero() const;
	T & addTo(T &, const T &) const;
	T & substractFrom(T &, const T &) const;
	T & multiplyByConstant(T &, const K &) const;
};

template<typename K>
class NormedVectorSpace<double, K> {
public:
	unsigned int getDimension() const{
		return 1;
	}

	K distance(const double &left, const double &right) const {
		return std::fabs(right - left);
	}

	double zero() const {
		return 0;
	}

	double & addTo(double &left, const double &right) const {
		left += right;
		return left;
	}

	double & substractFrom(double &left, const double &right) const {
		left -= right;
		return left;
	}

	double & multiplyByConstant(double & left, const K &constant) const {
		left *= constant;
		return left;
	}
};

template<typename T, typename K>
class NormedVectorSpace<std::vector<T>, K> {
public:
	unsigned int dimension;
	NormedVectorSpace<std::vector<T>, K>(unsigned int n) : dimension(n) {
	}

	unsigned int getDimension() const {
		return dimension;
	}

	K distance(const std::vector<T> &left, const std::vector<T> &right) const {
		//std::cout << dimension << std::endl;
		assert(dimension == left.size() && dimension == right.size());
		K dist = 0;
		for_each_two_ranges(left.cbegin(), left.cend(), right.cbegin(),
				[&dist](const T &leftData, const T &rightData) {
					T temp = (rightData-leftData);
					dist += temp*temp;
				});
		return std::sqrt(dist);
	}

	std::vector<T> zero() const {
		return std::vector<T>(dimension, 0);
	}

	std::vector<T> & addTo(std::vector<T> &left,
			const std::vector<T> &right) const {
		assert(dimension == left.size() && dimension == right.size());
		transform(left.begin(), left.end(), right.cbegin(), left.begin(),
				std::plus<T>());
		return left;
	}

	std::vector<T> & substractFrom(std::vector<T> &left,
			const std::vector<T> &right) const {
		assert(dimension == left.size() && dimension == right.size());
		transform(left.begin(), left.end(), right.cbegin(), left.begin(),
				std::minus<T>());
		return left;
	}

	std::vector<T> & multiplyByConstant(std::vector<T> & left,
			const K &constant) const {
		assert(dimension == left.size());
		transform(left.begin(), left.end(), left.begin(),
				[&constant](const T& data) {
					return data*constant;
				});
		return left;
	}
};

template<typename T, typename K = double>
class EuclideanSpace: public NormedVectorSpace<std::vector<T>, K> {
public:
	EuclideanSpace(unsigned int n) :
			NormedVectorSpace<std::vector<T>, K>(n) {
	}
};

template<typename T, typename K = double>
class ManhattanSpace: public NormedVectorSpace<std::vector<T>, K> {
public:
	ManhattanSpace(unsigned int n) :
			NormedVectorSpace<std::vector<T>, K>(n) {
	}

	K distance(const std::vector<T> &left, const std::vector<T> &right) const {
		assert(
				this->dimension == left.size()
						&& this->dimension == right.size());
		K dist = 0;
		for_each_two_ranges(left.cbegin(), left.cend(), right.cbegin(),
				[&dist](const T &leftData, const T &rightData) {
					dist += std::fabs((rightData-leftData));
				});
		return dist;
	}
};

#endif /* SRC_NORMEDVECTORSPACE_HPP_ */
