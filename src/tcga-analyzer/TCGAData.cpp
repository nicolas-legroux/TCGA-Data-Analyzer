/*
 * TCGAData.cpp
 *
 *  Created on: Jun 17, 2015
 *      Author: nicolas
 */

#include "../tcga-analyzer/TCGAData.hpp"

#include <fstream>
#include <iostream>
#include <set>
#include "../utilities.hpp"
#include "../config.hpp"

GeneList & TCGAData::getGeneListHandler() {
	return geneList;
}

const GeneList & TCGAData::getGeneListHandler() const {
	return geneList;
}

PatientNamesCancerMap & TCGAData::getPatientListHandler(bool getTumorData) {
	if (getTumorData) {
		return tumorPatientList;
	} else {
		return controlPatientList;
	}
}

const PatientNamesCancerMap & TCGAData::getPatientListHandler(
		bool getTumorData) const {
	if (getTumorData) {
		return tumorPatientList;
	} else {
		return controlPatientList;
	}
}

RNASeqDataCancerMap & TCGAData::getRNASeqDataHandler(bool getTumorData) {
	if (getTumorData) {
		return tumorRNASeqData;
	} else {
		return controlRNASeqData;
	}
}

const RNASeqDataCancerMap & TCGAData::getRNASeqDataHandler(
		bool getTumorData) const {
	if (getTumorData) {
		return tumorRNASeqData;
	} else {
		return controlRNASeqData;
	}
}

Eigen::MatrixXd & TCGAData::getDataMatrixHandler() {
	return dataMatrix;
}

const Eigen::MatrixXd & TCGAData::getDataMatrixHandler() const {
	return dataMatrix;
}

std::vector<Sample> & TCGAData::getSamplesHandler() {
	return samples;
}

const std::vector<Sample> & TCGAData::getSamplesHandler() const {
	return samples;
}

PatientIDsCancerMap & TCGAData::getPatientsIDsHandler() {
	return patientIDs;
}

const PatientIDsCancerMap & TCGAData::getPatientsIDsHandler() const {
	return patientIDs;
}

unsigned int TCGAData::getNumberOfGenes() const {
	return geneList.size();
}

unsigned int TCGAData::getNumberOfSamples() const {
	unsigned int numberOfSamples = 0;
	for (const auto &kv : controlRNASeqData) {
		numberOfSamples += kv.second.at(0).size();
	}
	for (const auto &kv : tumorRNASeqData) {
		numberOfSamples += kv.second.at(0).size();
	}
	return numberOfSamples;
}

std::vector<double> TCGAData::getPatientTumorData(const std::string &cancer,
		int patientIndex, bool getTumor) const {
	unsigned int numberOfGenes = getNumberOfGenes();
	std::vector<double> data(numberOfGenes);
	const std::vector<std::vector<double>> *dataHandler = 0;
	if (getTumor) {
		dataHandler = &tumorRNASeqData.at(cancer);
	} else {
		dataHandler = &controlRNASeqData.at(cancer);
	}

	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		data[i] = dataHandler->at(i).at(patientIndex);
	}

	return data;
}

void TCGAData::transposeData(bool verbose) {
	unsigned int numberOfGenes = getNumberOfGenes();
	unsigned int numberOfSamples = getNumberOfSamples();

	if (verbose) {
		std::cout << "Transposing data (number of samples : "
				<< numberOfSamples << ")... " << std::flush;
	}

	dataMatrix.resize(numberOfGenes, numberOfSamples);
	patientIDs.clear();
	samples.clear();

	int countPatients = 0;

	for (const auto &kv : tumorRNASeqData) {
		std::string cancerName = kv.first;
		patientIDs.insert(
				make_pair(cancerName + "_" + "Tumor", std::vector<int>()));
		patientIDs.insert(
				make_pair(cancerName + "_" + "Control", std::vector<int>()));

		for (unsigned int j = 0;
				j < controlRNASeqData.at(cancerName).at(0).size(); ++j) {

			patientIDs[cancerName + "_" + "Control"].push_back(countPatients);
			samples.push_back(
					Sample(cancerName, false,
							controlPatientList.at(cancerName).at(j)));
			for (unsigned int k = 0; k < numberOfGenes; ++k) {
				dataMatrix(k, countPatients) =
						controlRNASeqData.at(cancerName).at(k).at(j);
			}
			countPatients++;
		}

		for (unsigned int j = 0;
				j < tumorRNASeqData.at(cancerName).at(0).size(); ++j) {

			patientIDs[cancerName + "_" + "Tumor"].push_back(countPatients);
			samples.push_back(
					Sample(cancerName, true,
							tumorPatientList.at(cancerName).at(j)));
			for (unsigned int k = 0; k < numberOfGenes; ++k) {
				dataMatrix(k, countPatients) =
						tumorRNASeqData.at(cancerName).at(k).at(j);
			}
			countPatients++;
		}
	}

	if (verbose) {
		std::cout << "Done." << std::endl;
	}
}

void TCGAData::exportToMatrix(const std::string &matrixFilenamePath,
		const std::string &patientListFilenamePath, bool verbose) const {

	std::ofstream matrixOutputStream(EXPORT_DIRECTORY + matrixFilenamePath);
	std::ofstream patientsOutputStream(
			EXPORT_DIRECTORY + patientListFilenamePath);

	if (verbose) {
		std::cout << std::endl << "****** Exporting raw matrix ******"
				<< std::endl;
	}

	unsigned int numberOfGenes = getNumberOfGenes();

	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		if (verbose && i % 1000 == 0) {
			printAdvancement(i, numberOfGenes);
		}

		for (const auto &kv : tumorRNASeqData) {
			std::string cancerName = kv.first;
			for (unsigned int j = 0;
					j < controlRNASeqData.at(cancerName).at(0).size(); ++j) {
				if (i == 0) {
					patientsOutputStream << cancerName << "-Control ("
							<< controlPatientList.at(cancerName).at(j) << ")"
							<< std::endl;
				}
				matrixOutputStream
						<< controlRNASeqData.at(cancerName).at(i).at(j) << "\t";
			}

			for (unsigned int j = 0;
					j < tumorRNASeqData.at(cancerName).at(0).size(); ++j) {
				if (i == 0) {
					patientsOutputStream << cancerName << "-Tumor ("
							<< tumorPatientList.at(cancerName).at(j) << ")"
							<< std::endl;
				}
				matrixOutputStream << tumorRNASeqData.at(cancerName).at(i).at(j)
						<< "\t";
			}
		}
		matrixOutputStream << std::endl;
	}

	if (verbose) {
		std::cout << "********* Done exporting *********" << std::endl;
	}
}

void TCGAData::keepOnlyGenesInGraph(const std::string &filenameNodes) {
	std::cout << "The number of genes is currently " << getNumberOfGenes()
			<< ". ";
	std::cout << "Keeping only genes which are in the graph..." << std::endl;

	std::string pathToFile = GRAPH_DATA_DIRECTORY + filenameNodes;

	std::ifstream inputStream(pathToFile);
	GeneList newGeneList;
	RNASeqDataCancerMap newControlData;
	RNASeqDataCancerMap newTumorData;
	std::set<std::string> genesInGraph;
	std::string gene;

	while (inputStream >> gene) {
		genesInGraph.insert(gene);
	}

	for (unsigned int i = 0; i < geneList.size(); ++i) {
		std::string HGNCSymbol = geneList[i].first;
		if (genesInGraph.find(HGNCSymbol) != genesInGraph.end()) {
			newGeneList.push_back(geneList[i]);
		}
	}

	for (const auto &kv : controlRNASeqData) {
		std::string cancer = kv.first;
		newControlData.insert( { cancer, std::vector<std::vector<double>>() });
		for (unsigned int i = 0; i < geneList.size(); ++i) {
			std::string HGNCSymbol = geneList[i].first;
			if (genesInGraph.find(HGNCSymbol) != genesInGraph.end()) {
				newControlData[cancer].push_back(controlRNASeqData[cancer][i]);
			}
		}
	}

	controlRNASeqData = std::move(newControlData);

	for (const auto &kv : tumorRNASeqData) {
		std::string cancer = kv.first;
		newTumorData.insert( { cancer, std::vector<std::vector<double>>() });
		for (unsigned int i = 0; i < geneList.size(); ++i) {
			std::string HGNCSymbol = geneList[i].first;
			if (genesInGraph.find(HGNCSymbol) != genesInGraph.end()) {
				newTumorData[cancer].push_back(tumorRNASeqData[cancer][i]);
			}
		}
	}

	tumorRNASeqData = std::move(newTumorData);
	geneList = std::move(newGeneList);

	std::cout << "Done. Gene count is now " << getNumberOfGenes() << "."
			<< std::endl;
}

