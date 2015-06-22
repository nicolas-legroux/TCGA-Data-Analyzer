#ifndef DATAREADER_H_INCLUDED
#define DATAREADER_H_INCLUDED

#include <vector>
#include "TCGAData.hpp"

class TCGADataLoader {
public:
	TCGADataLoader() : ptrToData (nullptr), verbose(false), maxControlSamples(0), maxTumorSamples(0) { };
	TCGADataLoader(TCGAData *_ptrToData,
			const std::vector<std::string> &_cancers,
			unsigned int _maxControlSamples = 100,
			unsigned int _maxTumorSamples = 1000, bool verbose = true);
	TCGADataLoader(TCGAData *_ptrToData,
			const std::string &_filenameWithCancerList,
			unsigned int _maxControlSamples = 100,
			unsigned int _maxTumorSamples = 1000, bool verbose = true);
	void loadData();
private:
	std::vector<std::string> cancers;
	TCGAData *ptrToData;
	bool verbose;
	unsigned int maxControlSamples;
	unsigned int maxTumorSamples;

	void loadGeneData();
	void loadPatientDataByCancer(const std::string &cancer);
	void loadPatientData();
	void initializeRNASeqData(bool isTumorData);
	void loadRNASeqSample(const std::string &cancer, bool isTumorData,
			const std::string &patient);
	void loadRNASeqData(bool tumorData);
	void loadRNASeqData();
};

#endif // DATAREADER_H_INCLUDED
