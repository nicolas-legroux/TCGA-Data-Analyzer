/*
 * graph.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: nicolas
 */

#ifndef SRC_HEINZ_ANALYZER_GRAPH_HPP_
#define SRC_HEINZ_ANALYZER_GRAPH_HPP_

#include <vector>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <memory>

class PPIGraph {
public:
	typedef unsigned int NodeIDType;
	typedef double NodeValueType;
	typedef std::string NodeNameType;

	PPIGraph() = default;

	void addNode(const NodeNameType &name, NodeValueType value = NodeValueType());
	void addEdge(const NodeNameType &node1, const NodeNameType &node2);

	unsigned int size();
	unsigned int edgeCount();
	unsigned int getOutDegree(const NodeNameType &name);
	NodeValueType getNodeValue(const NodeNameType &name);

	static std::shared_ptr<PPIGraph> buildFromFile(const std::string &nodesFile,
			const std::string &edgesFile, bool labelledGraph = false);

	void printDegreeStatistics();
	void printNodes();
	void printNodesToFile(const std::string &filename);

	const std::vector<std::pair<NodeNameType, NodeValueType>> &getNodesHandler() const;
	std::vector<std::pair<NodeNameType, NodeValueType>> &getNodesHandler();

private:
	std::vector<std::pair<NodeNameType, NodeValueType>> nodes;
	std::map<NodeNameType, NodeIDType> nameToId;
	std::vector<std::set<NodeIDType>> adjacencyList;
	NodeIDType getNodeId(const NodeNameType &name);
};

#endif /* SRC_HEINZ_ANALYZER_GRAPH_HPP_ */
