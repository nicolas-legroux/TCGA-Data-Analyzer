/*
 * TCGA-Data.hpp
 *
 *  Created on: Jun 17, 2015
 *      Author: nicolas
 */

#ifndef SRC_TCGADATA_HPP_
#define SRC_TCGADATA_HPP_

#include <Eigen/Dense>
#include <map>
#include "TCGAPatientData.hpp"
#include "../tcga-analyzer/typedefs.hpp"

class TCGAData {
public:
	TCGAData() : dataMatrixIsComputed(false) {}
	//Handlers
	GeneList &getGeneListHandler();
	const GeneList &getGeneListHandler() const;
	std::vector<TCGAPatientData> &getPatientsHandler();
	const std::vector<TCGAPatientData> &getPatientsHandler() const;
	RNASeqData &getDataHandler();
	const RNASeqData &getDataHandler() const;
	Eigen::MatrixXd &getDataMatrixHandler();
	const Eigen::MatrixXd &getDataMatrixHandler() const;
	ClassMap &getClassMapHandler();
	const ClassMap &getClassMapHandler() const;

	//Utilities
	unsigned int getNumberOfGenes() const;
	unsigned int getNumberOfSamples() const;

	std::vector<double> getPatientRNASeqData(int patientIndex) const;
	std::vector<std::string> getPatientLabels() const;

	void setClinicalAttributes(const std::set<std::string> &attributes);

	void buildDataMatrix(bool verbose = false);

	void keepOnlyGenesInGraph(const std::string &filenameNodes);

	void reorderSamples();

private:
	GeneList geneList;
	std::vector<TCGAPatientData> patients;
	RNASeqData data;

	// Column = Patient; Row = Gene
	bool dataMatrixIsComputed;
	std::set<std::string> clinicalAttributes;
	Eigen::MatrixXd dataMatrix;
	ClassMap classMap;
};

#endif /* SRC_TCGADATA_HPP_ */
