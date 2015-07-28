/*
 * heinzAnalyzer.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: nicolas
 */

#ifndef SRC_HEINZ_ANALYZER_HEINZMODULEANALYZER_HPP_
#define SRC_HEINZ_ANALYZER_HEINZMODULEANALYZER_HPP_

#include <memory>
#include <string>
#include "graph.hpp"
#include "typedefs.hpp"

class HeinzModuleAnalyzer {
public:
	HeinzModuleAnalyzer(const std::string &_fileBasename);
	void printModule();
	void analyze(ClassCount *classCount, NegativeGeneCount *negativeGeneCount, DegreeStatistics *degreeStatistics);
private:
	std::string fileBasename;
	HeinzClass heinzClass;
	std::string patientID;
	std::shared_ptr<PPIGraph> module;
	int positives = 0;
	int negatives = 0;
	float meanDegree = 0;
	void buildModule();
};

#endif /* SRC_HEINZ_ANALYZER_HEINZMODULEANALYZER_HPP_ */
