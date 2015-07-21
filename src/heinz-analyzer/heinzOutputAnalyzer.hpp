/*
 * heinzOutputAnalyzer.hpp
 *
 *  Created on: Jul 20, 2015
 *      Author: nicolas
 */

#ifndef SRC_HEINZ_ANALYZER_HEINZOUTPUTANALYZER_HPP_
#define SRC_HEINZ_ANALYZER_HEINZOUTPUTANALYZER_HPP_

#include "typedefs.hpp"

class HeinzOutputAnalyzer{
public:
	HeinzOutputAnalyzer(const std::string &weightsFilename, const std::string &patientIDsFilename);
	void analyze();
private:
	std::vector<WeightType> weights;
	std::vector<std::string> patientIDs;
	ClassCount classCount;
	NegativeGeneCount negativeGeneCount;
	DegreeStatistics degreeStatistics;
};



#endif /* SRC_HEINZ_ANALYZER_HEINZOUTPUTANALYZER_HPP_ */
