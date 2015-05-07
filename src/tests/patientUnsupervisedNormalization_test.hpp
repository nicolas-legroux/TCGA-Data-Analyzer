/*
 * patientUnsupervisedNormalization_test.hpp
 *
 *  Created on: May 7, 2015
 *      Author: nicolas
 */

#ifndef SRC_TESTS_PATIENTUNSUPERVISEDNORMALIZATION_TEST_HPP_
#define SRC_TESTS_PATIENTUNSUPERVISEDNORMALIZATION_TEST_HPP_

#include "../correlationMatrix.hpp"
#include "../dataReader.hpp"
#include "../heatMap.hpp"
#include "../patientUnsupervisedNormalization.hpp"

void normalizationTest1(int K, int Nmax) {
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

		//NORMALIZE
		normalizeKMeans(controlData, tumorData, K, Nmax);

		std::vector<std::vector<double>> data;
		std::vector<DataIdentifier> dataIdentifiers;
		DataTypeMapping dataTypeMapping;

		prepareData(data, dataIdentifiers, dataTypeMapping, patientControlList,
				patientTumorList, controlData, tumorData);

		std::vector<double> correlationMatrixPearson = pearson(data);
		exportCorrelationMatrix(correlationMatrixPearson, dataIdentifiers, "matrix.pearson", "patients.pearson", "labels.pearson");
		exportGeneralStats(correlationMatrixPearson, dataTypeMapping, "classes_correlation_pearson.out", "classes_size_pearson.out");
		makeHeatMap(correlationMatrixPearson,"export/heat_map_pearson.png");

		std::vector<double> correlationMatrixSpearman = spearman(data);
		exportCorrelationMatrix(correlationMatrixSpearman, dataIdentifiers, "matrix.spearman", "patients.spearman", "labels.spearman");
		exportGeneralStats(correlationMatrixSpearman, dataTypeMapping, "classes_correlatio_spearman.out", "classes_size_spearman.out");
		makeHeatMap(correlationMatrixSpearman,"export/heat_map_spearman.png");
}

#endif /* SRC_TESTS_PATIENTUNSUPERVISEDNORMALIZATION_TEST_HPP_ */
