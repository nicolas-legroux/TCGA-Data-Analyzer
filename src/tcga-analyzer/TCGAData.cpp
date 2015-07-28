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

RNASeqData & TCGAData::getDataHandler() {
	return data;
}

const RNASeqData & TCGAData::getDataHandler() const {
	return data;
}

Eigen::MatrixXd & TCGAData::getDataMatrixHandler() {
	return dataMatrix;
}

const Eigen::MatrixXd & TCGAData::getDataMatrixHandler() const {
	return dataMatrix;
}

std::vector<TCGAPatientData> & TCGAData::getPatientsHandler() {
	return patients;
}

const std::vector<TCGAPatientData> & TCGAData::getPatientsHandler() const {
	return patients;
}

ClassMap & TCGAData::getClassMapHandler() {
	return classMap;
}

const ClassMap & TCGAData::getClassMapHandler() const {
	return classMap;
}

unsigned int TCGAData::getNumberOfGenes() const {
	return geneList.size();
}

unsigned int TCGAData::getNumberOfSamples() const {
	return patients.size();
}

std::vector<double> TCGAData::getPatientRNASeqData(int patientIndex) const {
	unsigned int numberOfGenes = getNumberOfGenes();
	std::vector<double> patientData(numberOfGenes);
	for (unsigned int i = 0; i < numberOfGenes; ++i) {
		patientData[i] = data[i][patientIndex];
	}

	return patientData;
}

std::vector<std::string> TCGAData::getPatientLabels() const {
	std::vector<std::string> v;
	for (const auto &sample : patients) {
		v.push_back(sample.toString());
	}
	return v;
}

void TCGAData::buildDataMatrix(bool verbose) {

	if (!dataMatrixIsComputed) {
		unsigned int numberOfGenes = getNumberOfGenes();
		unsigned int numberOfSamples = getNumberOfSamples();

		if (verbose) {
			std::cout << "Building data matrix (number of samples : "
					<< numberOfSamples << ")... " << std::flush;
		}

		dataMatrix.resize(numberOfGenes, numberOfSamples);
		classMap.clear();

		//Deal with gene data
		for (unsigned int i = 0; i < numberOfGenes; ++i) {
			for (unsigned int j = 0; j < numberOfSamples; ++j) {
				dataMatrix(i, j) = data[i][j];
			}
		}

		//Deal with patient data
		for (unsigned int j = 0; j < numberOfSamples; ++j) {
			classMap[patients[j].toClassString(clinicalKeys)].push_back(j);
		}

		if (verbose) {
			std::cout << "Done." << std::endl;
		}

		dataMatrixIsComputed = true;
	}
}

void TCGAData::keepOnlyGenesInGraph(const std::string &filenameNodes) {
	std::cout << "The number of genes is currently " << getNumberOfGenes()
			<< ". ";
	std::cout << "Keeping only genes which are in the graph..." << std::endl;

	std::string pathToFile = GRAPH_DATA_DIRECTORY + filenameNodes;

	std::ifstream inputStream(pathToFile);
	GeneList newGeneList;
	RNASeqData newData;
	std::set<std::string> genesInGraph;
	std::string gene;

	while (inputStream >> gene) {
		genesInGraph.insert(gene);
	}

	for (unsigned int i = 0; i < geneList.size(); ++i) {
		std::string HGNCSymbol = geneList[i].first;
		if (genesInGraph.find(HGNCSymbol) != genesInGraph.end()) {
			newGeneList.push_back(geneList[i]);
			newData.push_back(data[i]);
		}
	}

	data = std::move(newData);
	geneList = std::move(newGeneList);

	dataMatrixIsComputed = false;

	std::cout << "Done. Gene count is now " << getNumberOfGenes() << "."
			<< std::endl;
}

void TCGAData::reorderSamples() {
	RNASeqData newData;
	std::vector<TCGAPatientData> newPatients;
	newData.resize(getNumberOfGenes());

	dataMatrixIsComputed = false;

	for (const auto &kv : classMap) {
		for (auto i : kv.second) {
			newPatients.push_back(patients[i]);
			for (unsigned int j = 0; j < getNumberOfGenes(); ++j) {
				newData[j].push_back(data[j][i]);
			}
		}
	}

	data = std::move(newData);
	patients = std::move(newPatients);
	buildDataMatrix();
}
