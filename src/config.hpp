/*
 * config.hpp
 *
 *  Created on: Jun 17, 2015
 *      Author: nicolas
 */

#ifndef SRC_CONFIG_HPP_
#define SRC_CONFIG_HPP_

#include <string>
#include <set>

const std::set<std::string> ALLOWED_CANCERS =
		{ "BRCA", "COAD", "GBM", "HNSC", "KIRC", "KIRP", "LGG", "LUAD", "LUSC", "OV",
				"PRAD", "SARC", "THCA", "UCEC", "SARCB" };
const std::set<std::string> ALLOWED_METRICS =
		{ "absolute-cosine", "cosine-absolute-similarity", "cosine-distance",
				"cosine", "cosine-smilarity", "euclidean", "euclidean-distance",
				"squared-euclidean", "squared-euclidean-distance", "manhattan",
				"manhattan-distance", "absolute-pearson",
				"pearson-absolute-correlation", "pearson",
				"pearson-correlation", "pearson-distance", "absolute-spearman",
				"spearman-absolute-correlation", "spearman",
				"spearman-correlation", "spearman-distance", "jaccard",
				"jaccard-similarity", "jaccard-distance" };

enum UnsupervisedNormalizationMethod {
	KMEANS_NORMALIZATION, BINARY_QUANTILE_NORMALIZATION, NO_NORMALIZATION
};

const std::string TCGA_DATA_DIRECTORY = "data/tcga/";
const std::string GRAPH_DATA_DIRECTORY = "data/graph/";
const std::string HEINZ_DIRECTORY = "data/heinz/";
const std::string HEINZ_INPUT_DIRECTORY = HEINZ_DIRECTORY + "input/";
const std::string HEINZ_RAW_OUTPUT_DIRECTORY = HEINZ_DIRECTORY + "raw_output/";
const std::string HEINZ_OUTPUT_DIRECTORY = HEINZ_DIRECTORY + "output/";
const std::string HEINZ_SAMPLES_LIST = HEINZ_DIRECTORY + "samples.txt";
const std::string HEINZ_NEGATIVEWEIGHT_LIST = HEINZ_DIRECTORY
		+ "negative-weights.txt";

const std::string EXPORT_DIRECTORY = "export/";
const std::string SAMPLE_TCGA_FILE = TCGA_DATA_DIRECTORY
		+ "BRCA-normalized/TCGA-3C-AAAU-01.genes.normalized.results";
const std::string SAMPLE_BERGONIE_FILE = TCGA_DATA_DIRECTORY
		+ "SARCB-normalized/M976-01.genes.normalized.results";

#endif /* SRC_CONFIG_HPP_ */
