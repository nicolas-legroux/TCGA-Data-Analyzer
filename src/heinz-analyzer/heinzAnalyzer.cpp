/*
 * heinzAnalyzer.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: nicolas
 */

#include <iostream>
#include <regex>
#include <fstream>
#include <algorithm>
#include "heinzAnalyzer.hpp"
#include "../config.hpp"
#include "../tcga-analyzer/TCGADataLoader.hpp"

HeinzAnalyzer::HeinzAnalyzer(const std::string &_nodeFilename,
		const std::string &edgeFileCompleteGraph) :
		nodeFilename(_nodeFilename) {
	completeGraph = PPIGraph::buildFromFile(
			HEINZ_INPUT_DIRECTORY + nodeFilename + ".txt",
			GRAPH_DATA_DIRECTORY + edgeFileCompleteGraph, true);
	module = std::make_shared<PPIGraph>();
	buildModule();
}

void HeinzAnalyzer::buildModule() {
	std::ifstream inputStream(
			HEINZ_RAW_OUTPUT_DIRECTORY + nodeFilename + ".txt");
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
		//std::cout << line << std::endl;
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

void HeinzAnalyzer::analyze() {
	//Count positive nodes in complete graph
	unsigned int countPositivesInCompleteGraph = 0;
	std::for_each(completeGraph->getNodesHandler().cbegin(),
			completeGraph->getNodesHandler().cend(),
			[&](const std::pair<PPIGraph::NodeNameType, PPIGraph::NodeValueType> &pair) {
				if(pair.second>0) {
					++countPositivesInCompleteGraph;
				}
			});

	//Count positive nodes in module
	unsigned int countPositivesInModule = 0;
	std::for_each(module->getNodesHandler().cbegin(),
			module->getNodesHandler().cend(),
			[&](const std::pair<PPIGraph::NodeNameType, PPIGraph::NodeValueType> &pair) {
				if(pair.second>0) {
					++countPositivesInModule;
				}
			});

	//Count negative nodes in module
	unsigned int countNegativesInModule = module->size()
			- countPositivesInModule;

	//Count number of edges
	unsigned int numberOfEdgesInModule = module->edgeCount();

	//Count mean degree
	float meanDegree = (float) numberOfEdgesInModule / module->size();

	std::ofstream moduleStatistics(
			HEINZ_OUTPUT_DIRECTORY + nodeFilename + ".stats");
	moduleStatistics << countPositivesInCompleteGraph << " "
			<< countPositivesInModule << " " << countNegativesInModule
			<< std::endl;
	moduleStatistics << numberOfEdgesInModule << " " <<  meanDegree;

	std::ofstream nodesOutput(
			HEINZ_OUTPUT_DIRECTORY + nodeFilename + ".nodes");

	//Export list of nodes
	auto mapping = TCGADataLoader::buildHgnc2IdMapping();
	std::for_each(module->getNodesHandler().cbegin(),
			module->getNodesHandler().cend(),
			[&nodesOutput, &mapping](const std::pair<PPIGraph::NodeNameType, PPIGraph::NodeValueType> &pair) {
			 	 nodesOutput << mapping.at(pair.first) << std::endl;
			});
}
