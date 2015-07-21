/*
 * typedefs.hpp
 *
 *  Created on: Jul 20, 2015
 *      Author: nicolas
 */

#ifndef SRC_HEINZ_ANALYZER_TYPEDEFS_HPP_
#define SRC_HEINZ_ANALYZER_TYPEDEFS_HPP_

#include <string>
#include <vector>
#include <tuple>
#include <map>

using WeightType = std::string;
using CancerNameType = std::string;
using IsTumorType = bool;

using HeinzClass = std::tuple<WeightType, CancerNameType, IsTumorType>;

using NegativeGeneCount = std::map<HeinzClass, std::map<std::string, unsigned int>>;
using DegreeStatistics = std::map<HeinzClass, std::vector<float>>;
using ClassCount = std::map<HeinzClass, unsigned int>;

#endif /* SRC_HEINZ_ANALYZER_TYPEDEFS_HPP_ */
