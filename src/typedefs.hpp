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
typedef std::unordered_map<std::string, std::vector<std::vector<double>>> RNASeqDataCancerMap;
// cancerName -> patientID -> patientName
typedef std::unordered_map<std::string, std::vector<std::string>> PatientNamesCancerMap;
// geneID -> (HNSC Symbol, Entrez ID)
typedef std::vector<std::pair<std::string, int>> GeneList;
// cancerName -> list of Patient IDs
typedef std::unordered_map<std::string, std::vector<int>> PatientIDsCancerMap;

#endif /* SRC_TYPEDEFS_HPP_ */
