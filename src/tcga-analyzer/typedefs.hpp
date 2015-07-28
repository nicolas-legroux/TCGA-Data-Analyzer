/*
 * typedefs.hpp
 *
 *  Created on: May 11, 2015
 *      Author: nicolas
 */

#ifndef SRC_TYPEDEFS_HPP_
#define SRC_TYPEDEFS_HPP_

#include <unordered_map>
#include <vector>

// cancerName -> geneID -> patientID -> RNASeq value
typedef std::vector<std::vector<double>> RNASeqData;
// geneID -> (HNSC Symbol, Entrez ID)
typedef std::vector<std::pair<std::string, int>> GeneList;
// class -> list of Patient IDs
typedef std::map<std::string, std::vector<int>> ClassMap;

#endif /* SRC_TYPEDEFS_HPP_ */
