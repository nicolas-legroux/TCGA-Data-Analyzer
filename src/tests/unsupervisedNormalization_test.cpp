#include "../correlationMatrix.hpp"
#include "../dataReader.hpp"
#include "../heatMap.hpp"
#include "../unsupervisedNormalization.hpp"
#include "unsupervisedNormalization_test.hpp"

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
	CancerPatientIDList cancerPatientIDList;

	prepareData(data, sampleIdentifiers, cancerPatientIDList,
			patientControlList, patientTumorList, controlData, tumorData);

	std::vector<double> correlationMatrixPearson = pearson(data);
	exportCorrelationMatrix(correlationMatrixPearson, sampleIdentifiers,
			"matrix.pearson", "patients.pearson", "labels.pearson");
	exportClassStats(correlationMatrixPearson, cancerPatientIDList,
			sampleIdentifiers, "classes_correlation_pearson.tsv");
	makeHeatMap(correlationMatrixPearson, "heat_map_pearson.png",
			buildClassDivision(sampleIdentifiers), 20);

	if (K == 2) {
		printMaxExpressedGenes(controlData, tumorData, geneMapping, 15,
				"most_expressed_genes.out");
	}
}

void normalizationTestIteratedKMeans(int Niter) {
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
			controlData, tumorData, 50);

	//STEP2 : NORMALIZE
	normalizeIteratedBinaryKMeans(controlData, tumorData, Niter);

	//STEP3 : COMPUTE CORRELATION
	std::vector<std::vector<double>> data;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;

	prepareData(data, sampleIdentifiers, cancerPatientIDList,
			patientControlList, patientTumorList, controlData, tumorData);

	std::vector<double> correlationMatrixPearson = pearson(data);
	exportCorrelationMatrix(correlationMatrixPearson, sampleIdentifiers,
			"matrix.pearson", "patients.pearson", "labels.pearson");
	exportClassStats(correlationMatrixPearson, cancerPatientIDList,
			sampleIdentifiers, "classes_correlation_pearson.tsv");
	makeHeatMap(correlationMatrixPearson, "heat_map_pearson.png",
			buildClassDivision(sampleIdentifiers), 3);
}

void normalization_KMeans_Manhattan_test(int K, int Nmax) {
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
			controlData, tumorData, 50);

	//STEP2 : NORMALIZE
	normalizeKMeans(controlData, tumorData, K, Nmax);

	//STEP3 : COMPUTE CORRELATION
	std::vector<std::vector<double>> data;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;

	prepareData(data, sampleIdentifiers, cancerPatientIDList,
			patientControlList, patientTumorList, controlData, tumorData);

	std::vector<double> manhattanDistanceMatrix = manhattan(data);
	exportCorrelationMatrix(manhattanDistanceMatrix, sampleIdentifiers,
			"matrix.manhattan", "patients.manhattan", "labels.manhattan");
	exportClassStats(manhattanDistanceMatrix, cancerPatientIDList,
			sampleIdentifiers, "classes_manhattan.tsv");
	makeHeatMap(manhattanDistanceMatrix, "heat_map_manhattan.png",
			buildClassDivision(sampleIdentifiers), 3);
}

void no_normalization_euclidean_test() {
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
			controlData, tumorData, 50);

	//STEP3 : COMPUTE CORRELATION
	std::vector<std::vector<double>> data;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;

	prepareData(data, sampleIdentifiers, cancerPatientIDList,
			patientControlList, patientTumorList, controlData, tumorData);

	std::vector<double> euclideanDistanceMatrix = euclidean(data);
	exportCorrelationMatrix(euclideanDistanceMatrix, sampleIdentifiers,
			"matrix.manhattan", "patients.manhattan", "labels.manhattan");
	exportClassStats(euclideanDistanceMatrix, cancerPatientIDList,
			sampleIdentifiers, "classes_manhattan.tsv");
	makeHeatMap(euclideanDistanceMatrix, "heat_map_manhattan.png",
			buildClassDivision(sampleIdentifiers), 3);
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
			controlData, tumorData, 100);

	//STEP2 : NORMALIZE
	normalizeQuantile(controlData, tumorData, cutPercentage);

	//STEP3 : COMPUTE CORRELATION
	std::vector<std::vector<double>> data;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;

	prepareData(data, sampleIdentifiers, cancerPatientIDList,
			patientControlList, patientTumorList, controlData, tumorData);

	std::vector<double> correlationMatrixPearson = pearson(data);
	exportCorrelationMatrix(correlationMatrixPearson, sampleIdentifiers,
			"matrix.pearson", "patients.pearson", "labels.pearson");
	exportClassStats(correlationMatrixPearson, cancerPatientIDList,
			sampleIdentifiers, "classes_correlation_pearson.tsv");
	makeHeatMap(correlationMatrixPearson, "heat_map_pearson.png",
			buildClassDivision(sampleIdentifiers), 20);

	printMaxExpressedGenes(controlData, tumorData, geneMapping, 15,
			"most_expressed_genes.out");
}
