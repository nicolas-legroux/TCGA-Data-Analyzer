/*
 * CommandLineProcessor.cpp
 *
 *  Created on: Jun 29, 2015
 *      Author: nicolas
 */

#include <iostream>
#include <set>
#include <ClusterXX/metrics/metrics.hpp>
#include "command_line_processor.hpp"
#include "config.hpp"
#include "parameters.hpp"
#include "utilities.hpp"
#include "tcga-analyzer/TCGA-Analyzer.hpp"
#include "heinz-analyzer/heinzAnalyzer.hpp"

CommandLineProcessor::CommandLineProcessor(int argc, char *argv[]) {
	if (argc % 2 != 1) {
		throw wrong_usage_exception(
				"Invalid command line : number of arguments should be an odd number.");
	}

	for (int i = 1; i < argc; i += 2) {
		process(argv[i], argv[i + 1]);
	}
}

void CommandLineProcessor::process(const std::string &optionName,
		const std::string &optionValue) {
	if (optionName == "-mode") {
		if (std::isdigit(optionValue[0])) {
			int i = std::atoi(optionValue.c_str());
			if (i >= 0 && i <= 2) {
				PROGRAM_MODE = i;
			} else {
				throw wrong_usage_exception(
						"-mode option value should be a digit between 0 and 2.");

			}
		} else {
			throw wrong_usage_exception(
					"-mode option value should be a digit between 0 and 2.");
		}
	}

	else if (optionName == "-cancers") {
		CANCERS.clear();
		std::vector<std::string> cancers = split(optionValue, { ',' });
		for (const auto &s : cancers) {
			if (ALLOWED_CANCERS.find(s) != ALLOWED_CANCERS.end()) {
				CANCERS.insert(s);
			} else {
				throw wrong_usage_exception(
						"Error when trying to process cancer with name '" + s
								+ "'. -cancers must be a subset of "
								+ implode(ALLOWED_CANCERS.begin(),
										ALLOWED_CANCERS.end(), ","));
			}
		}
	}

	else if (optionName == "-weights") {
		WEIGHTS.clear();
		std::vector<std::string> weights = split(optionValue, { ',' });
		for (const auto &s : weights) {
			double d = std::atof(s.c_str());
			if (d <= 0) {
				throw wrong_usage_exception(
						"Error when trying to process weight '" + s
								+ "'. -weights should be a comma-separated list of weights (without the minus sign).");
			}
			WEIGHTS.push_back(d);
		}
	}

	else if (optionName == "-maxcontrol") {
		MAX_CONTROL_SAMPLES = std::atoi(optionValue.c_str());
	}

	else if (optionName == "-maxtumor") {
		MAX_TUMOR_SAMPLES = std::atoi(optionValue.c_str());
	}

	else if (optionName == "-verbose") {
		VERBOSE = std::atoi(optionValue.c_str());
	}

	else if (optionName == "-k") {
		K_CLUSTER = std::atoi(optionValue.c_str());
	}

	else if (optionName == "-cutpercentage") {
		BINARY_QUANTILE_NORMALIZATION_PARAM = std::atof(optionValue.c_str());
		if (BINARY_QUANTILE_NORMALIZATION_PARAM < 0
				|| BINARY_QUANTILE_NORMALIZATION_PARAM > 1) {
			throw wrong_usage_exception(
					"Error while processing " + optionName + " " + optionValue
							+ ".\n The cut percentage should be between 0 and 1");
		}
	}

	else if (optionName == "-metric") {
		if (ALLOWED_METRICS.find(optionValue) == ALLOWED_METRICS.end()) {
			throw wrong_usage_exception(
					"Error while processing " + optionName + " " + optionValue
							+ ".\n The metric name should be in the following set: \n{ "
							+ implode(ALLOWED_METRICS.begin(),
									ALLOWED_METRICS.end(), ", ") + " }");
		}
		METRIC = ClusterXX::buildMetric(optionValue);
	}

	else if(optionName == "-f"){
		workingFile = optionValue;
	}

	else {
		throw wrong_usage_exception(
				"Unknown command line option '" + optionName + "'.");
	}
}

void CommandLineProcessor::runProgram() {
	std::cout << "--------------------------------------" << std::endl;
	std::cout << "|            TCGA-ANALYZER           |" << std::endl;
	std::cout << "--------------------------------------" << std::endl;

	if (PROGRAM_MODE == 0) {
		std::cout << std::endl << "Program mode : 0 (Clustering mode)"
				<< std::endl << std::endl;
	}

	else if (PROGRAM_MODE == 1) {
		std::cout << std::endl << "Program mode : 1 (Entry of the Heinz pipeline)"
				<< std::endl << std::endl;
	}

	else if (PROGRAM_MODE == 2) {
		std::cout << std::endl << "Program mode : 2 (Heinz raw output analyzer)"
				<< std::endl << std::endl;
	}

	if (PROGRAM_MODE == 0 || PROGRAM_MODE == 1) {

		std::cout << "------------------- Data Parameters --------------------"
				<< std::endl;
		std::cout << "* Cancers : "
				<< implode(CANCERS.begin(), CANCERS.end(), ", ") << std::endl;
		std::cout << "* Max control samples : " << MAX_CONTROL_SAMPLES
				<< std::endl;
		std::cout << "* Max tumor samples : " << MAX_TUMOR_SAMPLES << std::endl;
		std::cout << "--------------------------------------------------------"
				<< std::endl << std::endl;

		/* Read Data */
		std::cout << "-------------------- Loading data ----------------------"
				<< std::endl;
		TCGAData data;
		TCGADataLoader loader(&data, CANCERS, MAX_CONTROL_SAMPLES,
				MAX_TUMOR_SAMPLES, VERBOSE);
		loader.loadData();
		data.keepOnlyGenesInGraph(GRAPH_NODE_FILE);
		std::cout << "--------------------------------------------------------"
				<< std::endl << std::endl;

		/* Normalize */
		std::cout << "-------------- Normalization parameters ----------------"
				<< std::endl;
		std::cout << "* Normalization method : Binary quantile" << std::endl;
		std::cout << "* Binary quantile cut percentage : "
				<< BINARY_QUANTILE_NORMALIZATION_PARAM << std::endl;
		std::cout << "--------------------------------------------------------"
				<< std::endl << std::endl;

		std::cout << "------------------ Normalizing data --------------------"
				<< std::endl;
		std::shared_ptr<Normalizer> normalizer = std::make_shared<
				BinaryQuantileNormalizer>(BINARY_QUANTILE_NORMALIZATION_PARAM);
		TCGADataNormalizer tcgaNormalizer(&data, normalizer, VERBOSE);
		tcgaNormalizer.normalize();
		std::cout << "--------------------------------------------------------"
				<< std::endl << std::endl;

		if (PROGRAM_MODE == 0) {
			/* Output distance matrix */
			std::cout
					<< "------------------ Distance matrix ---------------------"
					<< std::endl;
			std::cout << "* Metric : " << METRIC->toString() << std::endl;
			TCGADataDistanceMatrixAnalyser distanceMetricAnalyzer(&data, METRIC,
					VERBOSE);
			distanceMetricAnalyzer.computeDistanceMatrix();
			distanceMetricAnalyzer.exportClassStats();
			distanceMetricAnalyzer.exportHeatMap();
			std::cout
					<< "--------------------------------------------------------"
					<< std::endl << std::endl;

			std::cout
					<< "---------------- Clustering parameters -----------------"
					<< std::endl;
			if (K_CLUSTER == 0) {
				std::cout
						<< "Number of clusters to find : automatic (= number of real classes in the data)"
						<< std::endl;
			} else {
				std::cout << "Number of clusters to find : " << K_CLUSTER
						<< std::endl;
			}
			std::cout
					<< "--------------------------------------------------------"
					<< std::endl << std::endl;

			std::cout
					<< "------------------ KMeans Clustering -------------------"
					<< std::endl;

			TCGADataKMeansClusterer kMeansClusterer(&data, K_CLUSTER,
					K_MEANS_MAX_ITERATIONS, VERBOSE);
			kMeansClusterer.computeClustering();
			kMeansClusterer.printClusteringInfo();

			std::cout
					<< "--------------------------------------------------------"
					<< std::endl << std::endl;

			std::cout
					<< "--------------- Hierarchical Clustering ----------------"
					<< std::endl;

			TCGADataHierarchicalClusterer hierarchicalClusterer(&data,
					distanceMetricAnalyzer.getDistanceMatrixHandler(), METRIC,
					K_CLUSTER, DEFAULT_LINKAGE_METHOD, true);
			hierarchicalClusterer.computeClustering();
			hierarchicalClusterer.printClusteringInfo();

			std::cout
					<< "--------------------------------------------------------"
					<< std::endl << std::endl;

			std::cout
					<< "----------------- Spectral Clustering ------------------"
					<< std::endl;

			TCGADataSpectralClusterer spectralClusterer(&data,
					distanceMetricAnalyzer.getDistanceMatrixHandler(), METRIC,
					K_CLUSTER, DEFAULT_GRAPH_TRANSFORMATION, VERBOSE);
			spectralClusterer.computeClustering();
			spectralClusterer.printClusteringInfo();

			std::cout
					<< "--------------------------------------------------------"
					<< std::endl << std::endl;
		}

		else {
			std::ofstream negativeWeightsOutput(HEINZ_NEGATIVEWEIGHT_LIST);
			std::cout
					<< "----------------- Writing Heinz input ------------------"
					<< std::endl;
			std::vector<std::string> weights_string;
			for(double d : WEIGHTS){
				weights_string.push_back(removeTrailingZeros(std::to_string(d)));
			}
			std::cout << "* Weights : " << implode(weights_string.begin(), weights_string.end(), ", ") << std::endl;
			for (double d : WEIGHTS) {
				negativeWeightsOutput << removeTrailingZeros(std::to_string(d)) << std::endl;
				std::cout << "Writing files for d=-" << d << "... " << std::endl;
				tcgaNormalizer.exportToFile(1, -d);
			}

			std::cout
					<< "--------------------------------------------------------"
					<< std::endl << std::endl;
		}
	}

	else if (PROGRAM_MODE == 2) {
		if(workingFile.empty()){
			std::cout << "Missing -f option.";
		}
		else{
			HeinzAnalyzer heinzAnalyzer(workingFile, "biogrid-edges.txt");
			heinzAnalyzer.analyze();
		}
	}
}
