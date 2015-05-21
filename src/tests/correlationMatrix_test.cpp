#include "correlationMatrix_test.hpp"
#include "../correlationMatrix.hpp"
#include "../dataReader.hpp"
#include "../heatMap.hpp"

using namespace std;

void correlationMatrixTest1() {
	std::string filenameCancers = "cancer.list";
	vector<string> cancers{"BRCA","THCA"};
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
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;

	prepareData(data, sampleIdentifiers, cancerPatientIDList, patientControlList,
			patientTumorList, controlData, tumorData);

	std::vector<double> correlationMatrix = pearson(data);
	exportCorrelationMatrix(correlationMatrix, sampleIdentifiers, "matrix.pearson.out", "patients.out", "labels.pearson.out");
	exportClassStats(correlationMatrix, cancerPatientIDList, sampleIdentifiers, "classes_correlation_pearson.tsv");
	makeHeatMap(correlationMatrix,"heatmap_pearson.png", buildClassDivision(sampleIdentifiers), 20);
}
