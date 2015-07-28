#ifndef DATAREADER_H_INCLUDED
#define DATAREADER_H_INCLUDED

#include <vector>
#include <map>
#include <set>

#include "../tcga-analyzer/TCGAData.hpp"

class TCGADataLoader {
public:
	TCGADataLoader() : ptrToData (nullptr), verbose(false), maxControlSamples(0), maxTumorSamples(0) { };
	TCGADataLoader(TCGAData *_ptrToData,
			const std::set<std::string> &_cancers,
			unsigned int _maxControlSamples,
			unsigned int _maxTumorSamples, bool verbose);
	void loadGeneExpressionData(const std::string &sampleFilePath);
	void loadClinicalData(const std::set<std::string> &clinicalAttributes);

	static std::map<std::string, int> buildHgnc2IdMapping(const std::string &file);
private:
	std::set<std::string> cancers;
	TCGAData *ptrToData;
	bool verbose;
	unsigned int maxControlSamples;
	unsigned int maxTumorSamples;

	std::vector<std::string> clinicalKeys;
	std::map<std::string, std::map<std::string, std::string>> clinicalData;

	void loadGeneData(const std::string &file);
	void initializeRNASeqData();
	void loadRNASeqData(const std::string &cancer, const std::string &patientId);
	void loadDataByCancer(const std::string &cancer);

	void loadClinicalDataByCancer(const std::set<std::string> &clinicalAttributes, const std::string &cancer);
};

#endif // DATAREADER_H_INCLUDED
