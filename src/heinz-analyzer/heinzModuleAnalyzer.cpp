/*
 * heinzAnalyzer.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: nicolas
 */

#include "heinzModuleAnalyzer.hpp"

#include <iostream>
#include <regex>
#include <fstream>
#include <algorithm>
#include "../config.hpp"
#include "../tcga-analyzer/TCGADataLoader.hpp"
#include "../utilities.hpp"

HeinzModuleAnalyzer::HeinzModuleAnalyzer(const std::string &_fileBasename) :
		fileBasename(_fileBasename) {
	std::vector<std::string> v = split(_fileBasename, { '_' });
	heinzClass = HeinzClass(v[0], v[1], (v[2] == "Tumor"));
	patientID = v[3];
	module = std::make_shared<PPIGraph>();
	buildModule();
}

void HeinzModuleAnalyzer::buildModule() {
	std::ifstream inputStream(
			HEINZ_RAW_OUTPUT_DIRECTORY + std::get<0>(heinzClass) + "/"
					+ fileBasename + ".txt");
	std::string line;
	std::smatch match;
	std::string regexNodePattern =
			"([0-9]+) \\[label=\"([A-Z0-9\\-]+)\\\\n([-.0-9]+)\\\\n";
	std::string regexEdgePattern = "([0-9]+) -- ([0-9]+)";
	std::regex regexNode(regexNodePattern);
	std::regex regexEdge(regexEdgePattern);
	std::map<std::string, std::string> heinz2Hgnc;

	bool readingNodes = false;
	bool readingEdges = false;

	std::map<std::string, std::string> heinzToHGNC;

	while (getline(inputStream, line)) {
		if (!readingEdges && std::regex_search(line, match, regexNode)) {
			readingNodes = true;
			std::string heinzName = match[1];
			std::string hgncName = match[2];
			heinz2Hgnc[heinzName] = hgncName;
			double value = std::stod(match[3]);
			module->addNode(hgncName, value);
		}

		else if (readingNodes && std::regex_search(line, match, regexEdge)) {
			readingEdges = true;
			std::string heinzName1 = match[1];
			std::string heinzName2 = match[2];
			module->addEdge(heinz2Hgnc.at(heinzName1),
					heinz2Hgnc.at(heinzName2));
		}
	}
	//std::cout << "Done building module, size : " << module->size() << std::endl;
}

void HeinzModuleAnalyzer::analyze(ClassCount *classCount,
		NegativeGeneCount *negativeGeneCount,
		DegreeStatistics *degreeStatistics) {
	++(*classCount)[heinzClass];
	std::for_each(module->getNodesHandler().cbegin(),
			module->getNodesHandler().cend(),
			[&](const std::pair<PPIGraph::NodeNameType, PPIGraph::NodeValueType> &pair) {
				if(pair.second<0) {
					++(*negativeGeneCount)[heinzClass][pair.first];
					++negatives;
				}
			});
	//Count negative nodes in module
	positives = module->size() - negatives;
	//Count mean degree
	meanDegree = (float) module->edgeCount() / module->size();
	(*degreeStatistics)[heinzClass].push_back(meanDegree);
}

void HeinzModuleAnalyzer::printModule(){
	module->printNodesToFile(HEINZ_OUTPUT_DIRECTORY + fileBasename + ".txt");
}

