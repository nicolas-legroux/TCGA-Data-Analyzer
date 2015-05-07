/*
 * correlationMatrix_test.hpp
 *
 *  Created on: May 5, 2015
 *      Author: nicolas
 */

#ifndef SRC_TESTS_CORRELATIONMATRIX_TEST_HPP_
#define SRC_TESTS_CORRELATIONMATRIX_TEST_HPP_

#include "../correlationMatrix.hpp"
#include "../dataReader.hpp"

void correlationMatrixTest1() {
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

	std::vector<std::vector<double>> data;
	std::vector<DataIdentifier> dataIdentifiers;
	DataTypeMapping dataTypeMapping;

	prepareData(data, dataIdentifiers, dataTypeMapping, patientControlList,
			patientTumorList, controlData, tumorData);

	std::vector<double> correlationMatrix = pearson(data);
	exportCorrelationMatrix(correlationMatrix, dataIdentifiers, "matrix.out.test", "patients.out.test");
}

void correlationMatrixTest2() {
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
			controlData, tumorData, 500);

	std::vector<std::vector<double>> data;
	std::vector<DataIdentifier> dataIdentifiers;
	DataTypeMapping dataTypeMapping;

	prepareData(data, dataIdentifiers, dataTypeMapping, patientControlList,
			patientTumorList, controlData, tumorData);

	std::vector<double> correlationMatrixPearson = pearson(data);
	exportCorrelationMatrix(correlationMatrixPearson, dataIdentifiers, "matrix.pearson", "patients.pearson");
	exportGeneralStats(correlationMatrixPearson, dataTypeMapping, "general_stats_pearson.out");

	std::vector<double> correlationMatrixSpearman = spearman(data);
	exportCorrelationMatrix(correlationMatrixSpearman, dataIdentifiers, "matrix.spearman", "patients.spearman");
	exportGeneralStats(correlationMatrixSpearman, dataTypeMapping, "general_stats_spearman.out");
}

#endif /* SRC_TESTS_CORRELATIONMATRIX_TEST_HPP_ */
