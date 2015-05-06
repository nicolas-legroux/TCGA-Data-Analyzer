/*
 * correlationMatrix.hpp
 *
 *  Created on: May 5, 2015
 *      Author: nicolas
 */

#ifndef SRC_CORRELATIONMATRIX_HPP_
#define SRC_CORRELATIONMATRIX_HPP_

#include <string>
#include <unordered_map>
#include <vector>

struct DataIdentifier {
	std::string cancerName;
	bool isTumor;
	std::string patientId;

	DataIdentifier(std::string _cancerName, bool _isTumor,
			std::string _patientId) :
			cancerName(_cancerName), isTumor(_isTumor), patientId(_patientId) {
	}

	std::string toString() const {
		return cancerName + "-" + ((isTumor) ? "Tumor" : "Control") + " ("
				+ patientId + ")";
	}
};

typedef std::unordered_map<std::string, std::vector<int>> DataTypeMapping;
typedef std::unordered_map<std::string, std::vector<std::vector<double>>>RNASeqData;
typedef std::unordered_map<std::string, std::vector<std::string>> PatientList;

void prepareData(std::vector<std::vector<double>> &data,
		std::vector<DataIdentifier> &dataIdentifiers,
		DataTypeMapping &dataTypeMapping, const PatientList &patientControlList,
		const PatientList &patientTumorList, const RNASeqData &controlData,
		const RNASeqData &tumorData);

std::vector<double> pearson(std::vector<std::vector<double>> &data);
std::vector<double> spearman(std::vector<std::vector<double>> &data);

void exportCorrelationMatrix(const std::vector<double> &correlationMatrix,
		const std::vector<DataIdentifier> &dataIdentifiers,
		const std::string &filemaneMatrix, const std::string &patientsIds);

void exportGeneralStats(const std::vector<double> &correlationMatrix,
		const DataTypeMapping &dataTypeMapping, const std::string &filemane);

#endif /* SRC_CORRELATIONMATRIX_HPP_ */
