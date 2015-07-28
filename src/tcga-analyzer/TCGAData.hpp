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

/*Identifies uniquely a sample by knowing :
 *  - the Cancer set from which it comes from
 *  - whether it is a Control or Tumor Sample
 *  - the TCGA ID of the patient
 */

class TCGAData {
public:
	TCGAData() = default;
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

	void buildDataMatrix(const std::set<std::string> &keys = {}, bool verbose = true);

	void keepOnlyGenesInGraph(const std::string &filenameNodes);

	std::vector<std::string> getPatientLabels(){
		std::vector<std::string> v;
		for(const auto &sample : patients){
			v.push_back(sample.toString());
		}
		return v;
	}

private:
	GeneList geneList;
	std::vector<TCGAPatientData> patients;
	RNASeqData data;

	// Column = Patient; Row = Gene
	Eigen::MatrixXd dataMatrix;
	ClassMap classMap;
};

#endif /* SRC_TCGADATA_HPP_ */
