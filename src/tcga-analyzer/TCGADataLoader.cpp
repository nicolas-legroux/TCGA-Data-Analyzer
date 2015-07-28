#include "../tcga-analyzer/TCGADataLoader.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "../utilities.hpp"
#include "../config.hpp"

TCGADataLoader::TCGADataLoader(TCGAData *_ptrToData,
		const std::set<std::string> &_cancers, unsigned int _maxControlSamples,
		unsigned int _maxTumorSamples, bool _verbose) :
		cancers(_cancers), ptrToData(_ptrToData), verbose(_verbose), maxControlSamples(
				_maxControlSamples), maxTumorSamples(_maxTumorSamples) {
	//Nothing to do
}

std::map<std::string, int> TCGADataLoader::buildHgnc2IdMapping(
		const std::string &file) {
	std::map<std::string, int> mapping;
	std::fstream input(file);
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

void TCGADataLoader::loadGeneData(const std::string &file) {
	std::fstream input(file);
	std::string firstLine;
	std::getline(input, firstLine);
	std::string geneId;
	double score;
	while (input >> geneId >> score) {
		std::vector<std::string> strs = split(geneId, { '|' });
		std::string hgncSymbol = boost::to_upper_copy<std::string>(strs[0]);
		int entrezId = -1;
		if (strs.size() > 1) {
			entrezId = atoi(strs[1].c_str());
		}
		ptrToData->getGeneListHandler().push_back(
				make_pair(hgncSymbol, entrezId));
	}
}

void TCGADataLoader::initializeRNASeqData() {
	unsigned int numberOfGenes = ptrToData->getNumberOfGenes();
	ptrToData->getDataHandler().resize(numberOfGenes);
}

void TCGADataLoader::loadRNASeqData(const std::string &cancer,
		const std::string &patientId) {
	std::string filePath = TCGA_DATA_DIRECTORY + cancer + "-normalized/"
			+ patientId + ".genes.normalized.results";
	std::string geneId;
	std::string firstLine;
	double score;
	std::fstream input(filePath);
	std::getline(input, firstLine);
	int i = 0;
	while (input >> geneId >> score) {
		ptrToData->getDataHandler()[i].push_back(score);
		i++;
	}
}

void TCGADataLoader::loadDataByCancer(const std::string &cancer) {
	std::string patientListFilename = TCGA_DATA_DIRECTORY + cancer
			+ "-normalized/patient.list";
	std::ifstream input(patientListFilename);
	std::string patientId;
	unsigned int countControl = 0;
	unsigned int countTumor = 0;
	unsigned int countLoadedControl = 0;
	unsigned int countLoadedTumor = 0;

	if (verbose) {
		std::cout << "*Processing " << cancer << "..." << std::endl;
	}

	while (input >> patientId) {
		std::vector<std::string> strs = split(patientId, { '-' });
		unsigned int isTumorInfoPosition = strs.size() - 1;
		std::string patientName = strs[0];
		for (unsigned int i = 1; i < isTumorInfoPosition; ++i) {
			patientName += "-" + strs[i];
		}

		if (strs[isTumorInfoPosition] == "01"
				&& ++countTumor <= maxTumorSamples) {
			loadRNASeqData(cancer, patientId);
			ptrToData->getPatientsHandler().push_back(
					TCGAPatientData(patientName, cancer, true));
			++countLoadedTumor;
		} else if (strs[isTumorInfoPosition] == "11"
				&& ++countControl <= maxControlSamples) {
			loadRNASeqData(cancer, patientId);
			ptrToData->getPatientsHandler().push_back(
					TCGAPatientData(patientName, cancer, false));
			++countLoadedControl;
		}
	}

	if (verbose) {
		std::cout << "\tLoaded " << countLoadedControl
				<< " control samples (out of " << countControl
				<< " samples) and " << countLoadedTumor
				<< " tumor samples (out of " << countTumor << " samples)."
				<< std::endl;
	}
}

void TCGADataLoader::loadData(const std::string &sampleFilePath) {
	loadGeneData(sampleFilePath);
	initializeRNASeqData();
	for (const auto &cancer : cancers) {
		loadDataByCancer(cancer);
	}
}
