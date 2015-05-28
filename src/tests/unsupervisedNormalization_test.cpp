#include "unsupervisedNormalization_test.hpp"

#include "../dataReader.hpp"
#include "../heatMap.hpp"
#include "../unsupervisedNormalization.hpp"
#include "../distanceMatrix.hpp"
#include "../typedefs.hpp"

using namespace std;

void unsupervisedNormalization_test(
		const UnsupervisedNormalizationMethod &method,
		const UnsupervisedNormalizationParameters &parameters,
		const DistanceMetric &distanceMetric) {

	//STEP1 : READ DATA
	string filenameCancers = "cancer.list";
	Data data;
	readData(filenameCancers, data, 100, 1000);

	//STEP2 : NORMALIZE
	unsupervisedNormalization(data, method, parameters);

	//STEP3 : COMPUTE DISTANCE MATRIX
	std::vector<std::vector<double>> rawData;
	std::vector<SampleIdentifier> sampleIdentifiers;
	CancerPatientIDList cancerPatientIDList;
	data.transposeData(rawData, sampleIdentifiers, cancerPatientIDList);

	std::vector<double> distanceMatrix = computeDistanceMatrix(rawData,
			distanceMetric);

	exportClassStats(distanceMatrix, cancerPatientIDList, sampleIdentifiers,
			"classes_" + distanceMetricName(distanceMetric) + ".tsv");

	string heatMapFilename = "heatmap_" + distanceMetricName(distanceMetric)
			+ ".png";

	int divisionThickness = (int) (0.004 * (double) rawData.size()) + 1;

	makeHeatMap(distanceMatrix, heatMapFilename.c_str(),
			buildClassDivision(sampleIdentifiers), divisionThickness);

	if ((method == UnsupervisedNormalizationMethod::KMEANS && parameters.K == 2)
			|| method == UnsupervisedNormalizationMethod::BINARY_ITERATED_KMEANS
			|| method == UnsupervisedNormalizationMethod::BINARY_QUANTILE) {
		printMaxExpressedGenes(data, 15, "mostExpressedGenes.out");
	}
}
