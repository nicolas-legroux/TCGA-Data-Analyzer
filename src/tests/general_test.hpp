/*
 * general_tests.hpp
 *
 *  Created on: May 5, 2015
 *      Author: nicolas
 */

#ifndef SRC_TESTS_GENERAL_TEST_HPP_
#define SRC_TESTS_GENERAL_TEST_HPP_

#include "../dataReader.hpp"

//Use this to export a Matrix to compute correlations in R
void general_test1() {
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
			controlData, tumorData, 50);

	exportToMatrix(patientControlList, patientTumorList, controlData, tumorData,
			"matrix.out", "patients.out", geneMapping.size());

}

#endif /* SRC_TESTS_GENERAL_TEST_HPP_ */
