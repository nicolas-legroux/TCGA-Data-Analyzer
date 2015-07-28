/*
 * graph.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: nicolas
 */

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include "graph.hpp"
#include "../config.hpp"

PPIGraph::NodeIDType PPIGraph::getNodeId(const NodeNameType &name) {
	auto it = nameToId.find(name);
	if (it == nameToId.end()) {
		throw std::invalid_argument("The node is not in the graph.");
	}
	return it->second;
}

void PPIGraph::addNode(const NodeNameType &name, NodeValueType value) {
	if (nameToId.find(name) == nameToId.end()) {
		NodeIDType id = nodes.size();
		nodes.push_back( { name, value });
		nameToId.insert( { name, id });
		adjacencyList.push_back(std::set<NodeIDType>());
	}
}

void PPIGraph::addEdge(const NodeNameType &name1, const NodeNameType &name2) {
	NodeIDType id1 = getNodeId(name1);
	NodeIDType id2 = getNodeId(name2);
	adjacencyList[id1].insert(id2);
	adjacencyList[id2].insert(id1);
}

unsigned int PPIGraph::getOutDegree(const NodeNameType &name) {
	NodeIDType id = getNodeId(name);
	return adjacencyList[id].size();
}

PPIGraph::NodeValueType PPIGraph::getNodeValue(const NodeNameType &name) {
	NodeIDType id = getNodeId(name);
	return nodes[id].second;
}

unsigned int PPIGraph::size() {
	return nodes.size();
}

unsigned int PPIGraph::edgeCount() {
	unsigned int count = 0;
	std::for_each(adjacencyList.cbegin(), adjacencyList.cend(),
			[&](const std::set<NodeIDType> &s) {
				count += s.size();
			});
	return count / 2;
}

void PPIGraph::printDegreeStatistics() {
	std::map<int, int> degrees;
	std::for_each(adjacencyList.cbegin(), adjacencyList.cend(),
			[&](const std::set<NodeIDType> &s) {
				++degrees[s.size()];
			});
	std::for_each(degrees.cbegin(), degrees.cend(),
			[](const std::pair<int, int> &p) {
				std::cout << p.first << " " << p.second << std::endl;
			});
}

void PPIGraph::printNodes() {
	std::for_each(nodes.cbegin(), nodes.cend(),
			[](const std::pair<NodeNameType, NodeValueType> &pair) {
				std::cout << "{" << pair.first << ", " << pair.second << "}" << std::endl;
			});
}

std::shared_ptr<PPIGraph> PPIGraph::buildFromFile(const std::string &nodesFile,
		const std::string &edgesFile, bool labelledGraph) {
	std::shared_ptr<PPIGraph> graphPtr = std::make_shared<PPIGraph>();

	std::ifstream inputFileNodes(nodesFile);
	std::string line;

	while (std::getline(inputFileNodes, line)) {
		std::stringstream ss(line);
		NodeValueType value = NodeValueType();
		NodeNameType nodeName;
		if (labelledGraph) {
			ss >> nodeName >> value;
		} else {
			ss >> nodeName;
		}

		graphPtr->addNode(nodeName, value);
	}

	std::ifstream inputFileEdges(edgesFile);
	while (std::getline(inputFileEdges, line)) {
		std::stringstream ss(line);
		NodeNameType node1;
		NodeNameType node2;
		ss >> node1 >> node2;
		graphPtr->addEdge(node1, node2);
	}

	std::cout << "Built graph with " << graphPtr->size() << " nodes and "
			<< graphPtr->edgeCount() << " edges." << std::endl;

	return graphPtr;
}

const std::vector<std::pair<PPIGraph::NodeNameType, PPIGraph::NodeValueType>> &PPIGraph::getNodesHandler() const {
	return nodes;
}

std::vector<std::pair<PPIGraph::NodeNameType, PPIGraph::NodeValueType>> &PPIGraph::getNodesHandler() {
	return nodes;
}

void PPIGraph::printNodesToFile(const std::string &filename) {
	std::ofstream output(filename);
	std::for_each(nodes.cbegin(), nodes.cend(),
			[&](const std::pair<NodeNameType, NodeValueType> &pair) {
				output << pair.first << " " << pair.second << std::endl;
			});
}

