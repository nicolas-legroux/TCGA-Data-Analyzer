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
#include "../unsupervisedNormalization.hpp"

void normalizationTest1KMeans(int K, int Nmax) {
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
			controlData, tumorData, 500);

	//STEP2 : NORMALIZE
	normalizeKMeans(controlData, tumorData, K, Nmax);

	//STEP3 : COMPUTE CORRELATION
	std::vector<std::vector<double>> data;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList dataTypeMapping;

	prepareData(data, sampleIdentifiers, dataTypeMapping, patientControlList,
			patientTumorList, controlData, tumorData);

	std::vector<double> correlationMatrixPearson = pearson(data);
	exportCorrelationMatrix(correlationMatrixPearson, sampleIdentifiers,
			"matrix.pearson", "patients.pearson", "labels.pearson");
	exportClassStats(correlationMatrixPearson, dataTypeMapping,
			"classes_correlation_pearson.out");
	makeHeatMap(correlationMatrixPearson, "heat_map_pearson.png");

	std::vector<double> correlationMatrixSpearman = spearman(data);
	exportCorrelationMatrix(correlationMatrixSpearman, sampleIdentifiers,
			"matrix.spearman", "patients.spearman", "labels.spearman");
	exportClassStats(correlationMatrixSpearman, dataTypeMapping,
			"classes_correlatio_spearman.out");
	makeHeatMap(correlationMatrixSpearman, "heat_map_spearman.png");
}

void normalizationTestQuantile(double cutPercentage) {
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
			controlData, tumorData, 500);

	//STEP2 : NORMALIZE
	normalizeQuantile(controlData, tumorData, cutPercentage);

	//STEP3 : COMPUTE CORRELATION
	std::vector<std::vector<double>> data;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList dataTypeMapping;

	prepareData(data, sampleIdentifiers, dataTypeMapping, patientControlList,
			patientTumorList, controlData, tumorData);

	std::vector<double> correlationMatrixPearson = pearson(data);
	exportCorrelationMatrix(correlationMatrixPearson, sampleIdentifiers,
			"matrix.pearson", "patients.pearson", "labels.pearson");
	exportClassStats(correlationMatrixPearson, dataTypeMapping,
			"classes_correlation_pearson.out");
	makeHeatMap(correlationMatrixPearson, "heat_map_pearson.png");

	std::vector<double> correlationMatrixSpearman = spearman(data);
	exportCorrelationMatrix(correlationMatrixSpearman, sampleIdentifiers,
			"matrix.spearman", "patients.spearman", "labels.spearman");
	exportClassStats(correlationMatrixSpearman, dataTypeMapping,
			"classes_correlatio_spearman.out");
	makeHeatMap(correlationMatrixSpearman, "heat_map_spearman.png");
}

void normalizeAndPrintMaxGenes(int K, int Nmax) {
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
			controlData, tumorData, 500);

	//STEP2 : NORMALIZE
	normalizeKMeans(controlData, tumorData, K, Nmax);

	//STEP 3 : PRINT MAX EXPRESSED GENES
	printMaxExpressedGenes(controlData, tumorData, geneMapping);
}

#endif /* SRC_TESTS_PATIENTUNSUPERVISEDNORMALIZATION_TEST_HPP_ */
