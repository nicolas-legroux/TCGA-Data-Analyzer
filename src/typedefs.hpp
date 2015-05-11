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
#include <string>

// cancerName -> geneID -> patientID -> RNASeq
typedef std::unordered_map<std::string, std::vector<std::vector<double>>> RNASeqData;
// cancerName -> patientID -> patientName
typedef std::unordered_map<std::string, std::vector<std::string>> PatientList;
// geneID -> (HNSC Symbol, Entrez ID)
typedef std::vector<std::pair<std::string, int>> GeneList;
// cancerName -> list of Patient IDs
typedef std::unordered_map<std::string, std::vector<int>> DataTypeMapping;

#endif /* SRC_TYPEDEFS_HPP_ */
