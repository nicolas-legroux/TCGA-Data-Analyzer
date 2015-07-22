/*
 * TCGA-Data.hpp
 *
 *  Created on: Jun 17, 2015
 *      Author: nicolas
 */

#ifndef SRC_TCGADATA_HPP_
#define SRC_TCGADATA_HPP_

#include <Eigen/Dense>
#include "../tcga-analyzer/typedefs.hpp"
#include <map>

/*Identifies uniquely a sample by knowing :
 *  - the Cancer set from which it comes from
 *  - whether it is a Control or Tumor Sample
 *  - the TCGA ID of the patient
 */

struct Sample {
	std::string cancerName;
	bool isTumor;
	std::string patientId;

	Sample(std::string _cancerName, bool _isTumor, std::string _patientId) :
			cancerName(_cancerName), isTumor(_isTumor), patientId(_patientId) {
	}

	std::string toFullString() const {
		return cancerName + "_" + ((isTumor) ? "Tumor" : "Control") + "_"
				+ patientId;
	}

	std::string toClassString() const {
		return cancerName + "_" + ((isTumor) ? "Tumor" : "Control");
	}
};

class TCGAData {
public:
	TCGAData() = default;

	//Handlers
	GeneList &getGeneListHandler();
	const GeneList &getGeneListHandler() const;
	PatientNamesCancerMap &getPatientListHandler(
			bool getTumorData = true);
	const PatientNamesCancerMap &getPatientListHandler(
			bool getTumorData = true) const;
	RNASeqDataCancerMap &getRNASeqDataHandler(bool getTumorData = true);
	const RNASeqDataCancerMap &getRNASeqDataHandler(bool getTumorData =
			true) const;
	Eigen::MatrixXd &getDataMatrixHandler();
	const Eigen::MatrixXd &getDataMatrixHandler() const;
	std::vector<Sample> &getSamplesHandler();
	const std::vector<Sample> &getSamplesHandler() const;
	PatientIDsCancerMap &getPatientsIDsHandler();
	const PatientIDsCancerMap &getPatientsIDsHandler() const;

	//Utilities
	unsigned int getNumberOfGenes() const;
	unsigned int getNumberOfSamples() const;

	//Handler to copy data and get specific sample info
	std::vector<double> getPatientTumorData(const std::string &cancer,
			int patientIndex, bool getTumor = true) const;

	void transposeData(bool verbose = true);

	void exportToMatrix(const std::string &matrixFilename,
			const std::string &patientListFilename, bool verbose = true) const;

	void keepOnlyGenesInGraph(const std::string &filenameNodes);

	std::vector<std::string> getPatientLabels(){
		std::vector<std::string> v;
		for(const auto &sample : samples){
			v.push_back(sample.toFullString());
		}
		return v;
	}

	//To export data in TSV file

private:
	PatientNamesCancerMap controlPatientList;
	PatientNamesCancerMap tumorPatientList;
	GeneList geneList;
	RNASeqDataCancerMap controlRNASeqData;
	RNASeqDataCancerMap tumorRNASeqData;

	std::vector<Sample> samples;
	// Column = Patient; Row = Gene
	Eigen::MatrixXd dataMatrix;
	PatientIDsCancerMap patientIDs;
};

#endif /* SRC_TCGADATA_HPP_ */
