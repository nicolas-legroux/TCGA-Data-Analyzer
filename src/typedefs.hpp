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

	int getNumberOfProteins() {
		return geneList.size();
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

	void transposeData(std::vector<std::vector<double>> &data,
			std::vector<SampleIdentifier> &sampleIdentifiers,
			CancerPatientIDList &cancerPatientIDList) {

		std::cout << std::endl << "Transposing data... " << std::flush;

		int countPatients = 0;

		for (const auto &kv : tumorRNASeqData) {
			std::string cancerName = kv.first;
			cancerPatientIDList.insert(
					make_pair(cancerName + "-" + "Tumor", std::vector<int>()));
			cancerPatientIDList.insert(
					make_pair(cancerName + "-" + "Control",
							std::vector<int>()));

			unsigned int numberOfGenes = getNumberOfProteins();

			for (unsigned int j = 0;
					j < controlRNASeqData.at(cancerName).at(0).size(); ++j) {

				cancerPatientIDList[cancerName + "-" + "Control"].push_back(
						countPatients);
				sampleIdentifiers.push_back(
						SampleIdentifier(cancerName, false,
								controlPatientList.at(cancerName).at(j)));
				std::vector<double> patientData(numberOfGenes);
				for (unsigned int k = 0; k < numberOfGenes; ++k) {
					patientData[k] = controlRNASeqData.at(cancerName).at(k).at(j);
				}
				data.push_back(patientData);
				countPatients++;
			}

			for (unsigned int j = 0; j < tumorRNASeqData.at(cancerName).at(0).size();
					++j) {

				cancerPatientIDList[cancerName + "-" + "Tumor"].push_back(
						countPatients);
				sampleIdentifiers.push_back(
						SampleIdentifier(cancerName, true,
								tumorPatientList.at(cancerName).at(j)));
				std::vector<double> patientData(numberOfGenes);
				for (unsigned int k = 0; k < numberOfGenes; ++k) {
					patientData[k] = tumorRNASeqData.at(cancerName).at(k).at(j);
				}
				data.push_back(patientData);
				countPatients++;
			}
		}

		std::cout << "Done." << std::endl;
	}
};

#endif /* SRC_TYPEDEFS_HPP_ */
