#include "../dataReader.hpp"
#include "dataReader_test.hpp"

void exportDataTest1(){
	//STEP1 : read the data
	std::string filenameCancers = "cancer.list";
	PatientList patientControlList;
	PatientList patientTumorList;
	RNASeqData controlData;
	RNASeqData tumorData;
	GeneList geneMapping(
			makeGeneMapping(
					"data/BRCA-normalized/TCGA-A1-A0SJ-01.genes.normalized.results"));
	readPatientData(filenameCancers, patientControlList, patientTumorList);
	readRNASeqData(patientControlList, patientTumorList, geneMapping,
			controlData, tumorData, 5);

	exportToMatrix(patientControlList, patientTumorList, controlData, tumorData, "data_matrix.out", "patients.out", geneMapping.size());
}
