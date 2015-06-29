#include "../tcga-analyzer/TCGADataDistanceMatrixAnalyzer.hpp"

#include <fstream>
#include <algorithm>
#include <ClusterXX/metrics/metrics.hpp>
#include <ClusterXX/utils/heatMapBuilder.hpp>
#include "../config.hpp"
#include "../utilities.hpp"

TCGADataDistanceMatrixAnalyser::TCGADataDistanceMatrixAnalyser(
		TCGAData *_ptrToData, const std::shared_ptr<ClusterXX::Metric> &_metric,
		bool _verbose) :
		ptrToData(_ptrToData), metric(_metric), verbose(_verbose), matrixIsComputed(false) {
}

void TCGADataDistanceMatrixAnalyser::computeDistanceMatrix() {
	if (!matrixIsComputed) {
		ptrToData->transposeData(verbose);
		metric->setVerbose(verbose);
		distanceMatrix = metric->computeMatrix(ptrToData->getDataMatrixHandler());
		//std::cout << "Distance metric : " << ptrToData->getDataMatrixHandler() << std::endl << distanceMatrix;
		matrixIsComputed = true;
	}
}

void TCGADataDistanceMatrixAnalyser::exportDistanceMatrix() {
	if (verbose) {
		std::cout << std::endl << "Exporting Correlation matrix..."
				<< std::flush;
	}

	std::ofstream matrixOutputStream(
			EXPORT_DIRECTORY + "matrix-" + metric->toString() + ".txt");
	std::ofstream patientsOutputStream(
			EXPORT_DIRECTORY + "patients-" + metric->toString() + ".txt");

	unsigned int numberOfSamples = ptrToData->getNumberOfSamples();

	for (unsigned int i = 0; i < numberOfSamples; ++i) {
		for (unsigned int j = 0; j < numberOfSamples; ++j) {
			matrixOutputStream << distanceMatrix(i, j) << "\t";
		}
		matrixOutputStream << std::endl;
	}

	for (const Sample &sample : ptrToData->getSamplesHandler()) {
		patientsOutputStream << sample.toFullString() << std::endl;
	}

	std::cout << " Done." << std::endl;
}

void TCGADataDistanceMatrixAnalyser::exportHeatMap(bool withClassDivision,
		std::array<unsigned char, 3> separatorColor) {
	std::ofstream outputStreamLabels(
			EXPORT_DIRECTORY + "class-sizes-" + metric->toString()
					+ ".txt");
	if (verbose) {
		std::cout << "Making heat map... " << std::flush;
	}
	std::string current = "";
	int countCurrent = 0;

	for (const Sample &sample : ptrToData->getSamplesHandler()) {
		std::string newCurrent = sample.toClassString();
		if (newCurrent != current) {
			if (countCurrent != 0) {
				outputStreamLabels << current << " " << countCurrent
						<< std::endl;
			}
			current = newCurrent;
			countCurrent = 1;
		} else {
			countCurrent++;
		}
	}

	outputStreamLabels << current << " " << countCurrent << std::endl;

	ClusterXX::Utilities::HeatMapBuilder heatMapBuilder;
	std::vector<unsigned int> classDivision;
	unsigned int lineThickness = 0;

	if (withClassDivision) {
		classDivision = buildClassDivisionForHeatmap();
		lineThickness = (unsigned int) (0.005
				* (double) ptrToData->getNumberOfSamples()) + 1;
	}

	heatMapBuilder.build(distanceMatrix,
			EXPORT_DIRECTORY + "heatmap-" + metric->toString() + ".png",
			classDivision, lineThickness, separatorColor);

	if (verbose) {
		std::cout << "Done." << std::endl;
	}
}

void TCGADataDistanceMatrixAnalyser::exportClassStats() {

	if (verbose) {
		std::cout << std::endl << "Exporting class stats... " << std::flush;
	}

	std::vector<std::string> classes;
	unsigned int numberOfSamples = distanceMatrix.cols();

	for (const auto &sample : ptrToData->getSamplesHandler()) {
		classes.push_back(sample.toClassString());
	}

	auto end_unique = unique(classes.begin(), classes.end());
	classes.erase(end_unique, classes.end());

	unsigned int numberOfClasses = classes.size();
	std::vector<double> mean_correlation(numberOfClasses * numberOfClasses);
	std::vector<double> standard_dev_correlation(
			numberOfClasses * numberOfClasses);

	for (unsigned int i = 0; i < numberOfClasses; ++i) {
		for (unsigned int j = i; j < numberOfClasses; ++j) {
			std::vector<double> data;
			for (int I : ptrToData->getPatientsIDsHandler().at(classes[i])) {
				for (int J : ptrToData->getPatientsIDsHandler().at(classes[j])) {
					// When I = J : we are comparing the same patients,
					//  we know the distance is null
					if (I != J) {
						double d = distanceMatrix(I, J);
						data.push_back(d);
					}
				}
			}

			double mean = computeMean(data);
			double standard_dev = computeStandardDeviation(data, true);
			mean_correlation[numberOfClasses * i + j] = mean;
			mean_correlation[numberOfClasses * j + i] = mean;
			standard_dev_correlation[numberOfClasses * i + j] = standard_dev;
			standard_dev_correlation[numberOfClasses * j + i] = standard_dev;
		}
	}

	std::ofstream outputStream(
			EXPORT_DIRECTORY + "class-statistics" + metric->toString()
					+ ".tsv");
	outputStream << "CLASSES";
	for (const std::string &s : classes) {
		outputStream << "\t" << s << " (" << ptrToData->getPatientsIDsHandler().at(s).size()
				<< ")";
	}
	outputStream << std::endl;

	for (unsigned int i = 0; i < numberOfClasses; ++i) {
		outputStream << classes[i] << " ("
				<< ptrToData->getPatientsIDsHandler().at(classes[i]).size() << ")";
		for (unsigned int j = 0; j < numberOfClasses; ++j) {
			outputStream << "\t" << mean_correlation[numberOfClasses * i + j]
					<< " (" << standard_dev_correlation[numberOfClasses * i + j]
					<< ")";
		}
		outputStream << std::endl;
	}

	if(verbose){
		std::cout << "Done." << std::endl;
	}
}

std::vector<unsigned int> TCGADataDistanceMatrixAnalyser::buildClassDivisionForHeatmap() {
	std::vector<unsigned int> classDivision;
	Sample prev = ptrToData->getSamplesHandler()[0];
	unsigned int count = 1;
	for (auto it = ptrToData->getSamplesHandler().begin() + 1; it != ptrToData->getSamplesHandler().end();
			++it) {
		Sample next = *it;
		if (prev.cancerName != next.cancerName
				|| prev.isTumor != next.isTumor) {
			classDivision.push_back(count);
		}
		++count;
		prev = next;
	}
	return classDivision;
}
