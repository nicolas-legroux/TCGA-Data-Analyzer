/*
 * heinzOutputAnalyzer.cpp
 *
 *  Created on: Jul 20, 2015
 *      Author: nicolas
 */

#include "heinzOutputAnalyzer.hpp"
#include "heinzModuleAnalyzer.hpp"
#include "../config.hpp"
#include "../utilities.hpp"
#include <fstream>
#include <iostream>

HeinzOutputAnalyzer::HeinzOutputAnalyzer(const std::string &weightsFilename,
		const std::string &patientIDsFilename) {
	std::string line;
	std::ifstream weightsFile(HEINZ_DIRECTORY + weightsFilename);
	while (weightsFile >> line) {
		weights.push_back(line);
	}
	std::ifstream samplesFile(HEINZ_DIRECTORY + patientIDsFilename);
	while (samplesFile >> line) {
		patientIDs.push_back(line);
	}
}

void HeinzOutputAnalyzer::analyze() {
	for (unsigned int i = 0; i < 5; ++i) {
		for (unsigned int j = 0; j < patientIDs.size(); ++j) {
			std::cout << patientIDs[j] << std::endl;
			HeinzModuleAnalyzer hma(weights[i] + "_" + patientIDs[j]);
			hma.analyze(&classCount, &negativeGeneCount, &degreeStatistics);
		}
	}

	for (const auto &kv : negativeGeneCount) {
		unsigned int count = classCount[kv.first];
		std::cout << kv.first << " (" << count << " samples)" << std::endl;
		for (const auto &kv2 : kv.second) {
			float p = (float) kv2.second / count;
			if (p >= 0.05) {
				std::cout << "\t" << kv2.first << " " << (100*p) << "%"
						<< std::endl;
			}
		}
	}
}

