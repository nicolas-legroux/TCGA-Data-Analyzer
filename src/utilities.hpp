#ifndef SRC_UTILITIES_HPP_
#define SRC_UTILITIES_HPP_

#include <vector>
#include <string>
#include <tuple>
#include <iostream>

//Prints advancement of a task in %
void printAdvancement(unsigned int currentCount, unsigned int totalCount);

//Splits a string according to delimiters
std::vector<std::string> split(const std::string &s,
		const std::vector<char> &delimiters);

template<typename Iter>
std::string implode(Iter begin, Iter end, const std::string &delimiter){
	std::string result;
	for(Iter i = begin; i != end; ++i){
		if(i != begin){
			result += delimiter;
		}
		result += *i;
	}
	return result;
}

template<std::size_t> struct int_{};

template <class Tuple, size_t Pos>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<Pos> ) {
  out << std::get< std::tuple_size<Tuple>::value-Pos >(t) << ',';
  return print_tuple(out, t, int_<Pos-1>());
}

template <class Tuple>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<1> ) {
  return out << std::get<std::tuple_size<Tuple>::value-1>(t);
}

template <class... Args>
std::ostream& operator<<(std::ostream& out, const std::tuple<Args...>& t) {
  out << '(';
  print_tuple(out, t, int_<sizeof...(Args)>());
  return out << ')';
}


std::string removeTrailingZeros(std::string s);

double computeMean(const std::vector<double> &vec);
double computeStandardDeviation(const std::vector<double> &vec, bool correction = true);

#endif /* SRC_UTILITIES_HPP_ */
