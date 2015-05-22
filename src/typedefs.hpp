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

// cancerName -> geneID -> patientID -> RNASeq value
typedef std::unordered_map<std::string, std::vector<std::vector<double>>> RNASeqData;
// cancerName -> patientID -> patientName
typedef std::unordered_map<std::string, std::vector<std::string>> PatientList;
// geneID -> (HNSC Symbol, Entrez ID)
typedef std::vector<std::pair<std::string, int>> GeneList;
// cancerName -> list of Patient IDs
typedef std::unordered_map<std::string, std::vector<int>> CancerPatientIDList;

struct SampleIdentifier {
	std::string cancerName;
	bool isTumor;
	std::string patientId;

	SampleIdentifier(std::string _cancerName, bool _isTumor,
			std::string _patientId) :
			cancerName(_cancerName), isTumor(_isTumor), patientId(_patientId) {
	}

	std::string toString() const {
		return cancerName + "-" + ((isTumor) ? "Tumor" : "Control") + " ("
				+ patientId + ")";
	}
};

struct Data {
	PatientList controlPatientList;
	PatientList tumorPatientList;
	GeneList geneList;
	RNASeqData controlRNASeqData;
	RNASeqData tumorRNASeqData;
};

#endif /* SRC_TYPEDEFS_HPP_ */
