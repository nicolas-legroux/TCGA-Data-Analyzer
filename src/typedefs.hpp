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
#include <iostream>

#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorX;

// cancerName -> geneID -> patientID -> RNASeq value
typedef std::unordered_map<std::string, std::vector<std::vector<double>>>RNASeqData;
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

	std::vector<SampleIdentifier> sampleIdentifiers;
	// Column = Patient; Row = Gene
	Eigen::MatrixXd transposedData;
	CancerPatientIDList cancerPatientIDList;

	unsigned int getNumberOfProteins() {
		return geneList.size();
	}

	unsigned int getNumberOfSamples() {
		unsigned int numberOfSamples = 0;
		for (const auto &kv : controlPatientList) {
			numberOfSamples += kv.second.size();
		}
		for (const auto &kv : tumorPatientList) {
			numberOfSamples += kv.second.size();
		}
		return numberOfSamples;
	}

	std::vector<double> getPatientTumorData(std::string &cancer,
			int patientId) {
		int numberOfProteins = getNumberOfProteins();
		std::vector<double> data(numberOfProteins);
		for (int i = 0; i < numberOfProteins; ++i) {
			data[i] = tumorRNASeqData.at(cancer).at(i).at(patientId);
		}
		return data;
	}

	void transposeData() {

		std::cout << std::endl << "Transposing data... " << std::flush;
		unsigned int numberOfGenes = getNumberOfProteins();
		unsigned int numberOfSamples = getNumberOfSamples();
		transposedData.resize(numberOfGenes, numberOfSamples);

		int countPatients = 0;

		for (const auto &kv : tumorRNASeqData) {
			std::string cancerName = kv.first;
			cancerPatientIDList.insert(
					make_pair(cancerName + "-" + "Tumor", std::vector<int>()));
			cancerPatientIDList.insert(
					make_pair(cancerName + "-" + "Control",
							std::vector<int>()));

			for (unsigned int j = 0;
					j < controlRNASeqData.at(cancerName).at(0).size(); ++j) {

				cancerPatientIDList[cancerName + "-" + "Control"].push_back(
						countPatients);
				sampleIdentifiers.push_back(
						SampleIdentifier(cancerName, false,
								controlPatientList.at(cancerName).at(j)));
				for (unsigned int k = 0; k < numberOfGenes; ++k) {
					transposedData(k, countPatients) = controlRNASeqData.at(
							cancerName).at(k).at(j);
				}
				countPatients++;
			}

			for (unsigned int j = 0;
					j < tumorRNASeqData.at(cancerName).at(0).size(); ++j) {

				cancerPatientIDList[cancerName + "-" + "Tumor"].push_back(
						countPatients);
				sampleIdentifiers.push_back(
						SampleIdentifier(cancerName, true,
								tumorPatientList.at(cancerName).at(j)));
				for (unsigned int k = 0; k < numberOfGenes; ++k) {
					transposedData(k, countPatients) = tumorRNASeqData.at(
							cancerName).at(k).at(j);
				}
				countPatients++;
			}
		}

		std::cout << "Done." << std::endl;
	}
};

#endif /* SRC_TYPEDEFS_HPP_ */
