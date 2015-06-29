/*
 * heinzAnalyzer.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: nicolas
 */

#ifndef SRC_HEINZ_ANALYZER_HEINZANALYZER_HPP_
#define SRC_HEINZ_ANALYZER_HEINZANALYZER_HPP_

#include <memory>
#include <string>
#include "graph.hpp"

class HeinzAnalyzer {
public:
	HeinzAnalyzer(const std::string &_nodeFilename,
			const std::string &edgeFileCompleteGraph);
	void analyze();
private:
	std::string nodeFilename;
	std::shared_ptr<PPIGraph> completeGraph;
	std::shared_ptr<PPIGraph> module;
	void buildModule();
};

#endif /* SRC_HEINZ_ANALYZER_HEINZANALYZER_HPP_ */
