#include "../tcga-analyzer/TCGADataLoader.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "../utilities.hpp"
#include "../config.hpp"

TCGADataLoader::TCGADataLoader(TCGAData *_ptrToData,
		const std::set<std::string> &_cancers,
		unsigned int _maxControlSamples, unsigned int _maxTumorSamples,
		bool _verbose) :
		cancers(_cancers), ptrToData(_ptrToData), verbose(_verbose), maxControlSamples(
				_maxControlSamples), maxTumorSamples(_maxTumorSamples) {
	//Nothing to do
}

TCGADataLoader::TCGADataLoader(TCGAData *_ptrToData,
		const std::string &_filenameWithCancerList,
		unsigned int _maxControlSamples, unsigned int _maxTumorSamples,
		bool _verbose) :
		ptrToData(_ptrToData), verbose(_verbose), maxControlSamples(
				_maxControlSamples), maxTumorSamples(_maxTumorSamples) {
	std::ifstream input(TCGA_DATA_DIRECTORY + _filenameWithCancerList);
	std::string cancerName;
	while (input >> cancerName) {
		cancers.insert(cancerName);
	}
}

std::map<std::string, int> TCGADataLoader::buildHgnc2IdMapping() {
	std::map<std::string, int> mapping;
	std::fstream input(SAMPLE_TCGA_FILE);
	std::string firstLine;
	std::getline(input, firstLine);
	std::string geneId;
	double score;
	int count = 0;
	while (input >> geneId >> score) {
		std::vector<std::string> strs = split(geneId,
				std::vector<char> { '|' });
		std::string hgncSymbol = boost::to_upper_copy<std::string>(strs[0]);
		mapping.insert( { hgncSymbol, count });
		++count;
	}
	return mapping;
}

void TCGADataLoader::loadGeneData() {
	std::fstream input(SAMPLE_TCGA_FILE);
	std::string firstLine;
	std::getline(input, firstLine);
	std::string geneId;
	double score;
	while (input >> geneId >> score) {
		std::vector<std::string> strs = split(geneId,
				std::vector<char> { '|' });
		std::string hgncSymbol = boost::to_upper_copy<std::string>(strs[0]);
		int entrezId = atoi(strs[1].c_str());
		ptrToData->getGeneListHandler().push_back(
				make_pair(hgncSymbol, entrezId));
	}
}

void TCGADataLoader::loadPatientDataByCancer(const std::string &cancer) {

	std::string patientListFilename = TCGA_DATA_DIRECTORY + cancer
			+ "-normalized/patient.list";
	std::ifstream input(patientListFilename);
	std::string patientId;
	int countControl = 0;
	int countTumor = 0;

//	if (verbose) {
//		std::cout << "**** Processing " << cancer << " Patients ****"
//				<< std::endl;
//	}

	while (input >> patientId) {
		std::vector<std::string> strs = split(patientId, std::vector<char> {
				'-', '.' });
		if (strs.size() >= 4) {
			std::string patientName = strs[0] + "-" + strs[1] + "-" + strs[2];
			if (strs[3] == "01") {
				ptrToData->getPatientListHandler(true)[cancer].push_back(
						patientName);
				++countTumor;
			} else if (strs[3] == "11") {
				ptrToData->getPatientListHandler(false)[cancer].push_back(
						patientName);
				++countControl;
			}
//				else if (verbose) {
//				std::cout << "Did not process " << patientId << std::endl;
//			}
		}
	}

//	if (verbose) {
//		std::cout << "Data types : " << std::endl;
//		std::cout << "01 (Tumor) --> " << countTumor << std::endl;
//		std::cout << "11 (Control) --> " << countControl << std::endl;
//
//		//Checking that for normal samples we have the corresponding cancer sample
//		for (const std::string &s : ptrToData->getPatientListHandler(false)[cancer]) {
//			if (find(ptrToData->getPatientListHandler(true)[cancer].begin(),
//					ptrToData->getPatientListHandler(true)[cancer].end(), s)
//					== ptrToData->getPatientListHandler(true)[cancer].end()) {
//				std::cout
//						<< "Did not find a matching tumor sample for the control sample "
//						<< s << std::endl;
//			}
//		}
//
//		std::cout << "*************************" << std::endl << std::endl;
//	}
}

void TCGADataLoader::loadPatientData() {
	for (const std::string &cancer : cancers) {
		ptrToData->getPatientListHandler(false).insert(
				std::make_pair(cancer, std::vector<std::string>()));
		ptrToData->getPatientListHandler(true).insert(
				std::make_pair(cancer, std::vector<std::string>()));
		loadPatientDataByCancer(cancer);
	}
}

void TCGADataLoader::initializeRNASeqData(bool isTumorData) {
	unsigned int numberOfGenes = ptrToData->getNumberOfGenes();
	for (auto itr = ptrToData->getPatientListHandler(isTumorData).begin();
			itr != ptrToData->getPatientListHandler(isTumorData).end(); ++itr) {
		ptrToData->getRNASeqDataHandler(isTumorData).insert(
				std::make_pair(itr->first,
						std::vector<std::vector<double>>(numberOfGenes)));

	}
}

void TCGADataLoader::loadRNASeqSample(const std::string &cancer,
		bool isTumorData, const std::string &patient) {
	std::string filename = patient + "-" + ((isTumorData) ? "01" : "11")
			+ ".genes.normalized.results";
	std::string filePath = TCGA_DATA_DIRECTORY + cancer + "-normalized/"
			+ filename;
	std::string geneId;
	std::string firstLine;
	double score;
	std::fstream input(filePath);
	std::getline(input, firstLine);

	int i = 0;
	while (input >> geneId >> score) {
		ptrToData->getRNASeqDataHandler(isTumorData)[cancer][i].push_back(
				score);
		i++;
	}
}

void TCGADataLoader::loadRNASeqData(bool tumorData) {
	initializeRNASeqData(tumorData);

	unsigned int maxSamples = (tumorData) ? maxTumorSamples : maxControlSamples;

	for (const auto &pairedData : ptrToData->getPatientListHandler(tumorData)) {
		const std::string &cancer = pairedData.first;
		if (verbose) {
			std::string type;
			(tumorData) ? type = "Tumor" : type = "Control";
			std::cout << "Reading " << cancer << " " << type
					<< " patient files... " << std::flush;
		}

		unsigned int i = 0;
		for (const std::string &patient : pairedData.second) {
			if (i >= maxSamples) {
				break;
			}
			loadRNASeqSample(cancer, tumorData, patient);
			i++;
		}

		if (verbose) {
			std::cout << "Found "
					<< ptrToData->getRNASeqDataHandler(tumorData)[cancer][0].size()
					<< " patients." << std::endl;
		}
	}
}

void TCGADataLoader::loadRNASeqData() {
	//Load Control Data
	loadRNASeqData(false);
	//Load tumor Data
	loadRNASeqData(true);
}

void TCGADataLoader::loadData() {
	loadGeneData();
	loadPatientData();
	loadRNASeqData();
}
