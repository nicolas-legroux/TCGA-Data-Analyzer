#ifndef DATAREADER_H_INCLUDED
#define DATAREADER_H_INCLUDED

#include <vector>
#include <string>
#include <unordered_map>

typedef std::unordered_map<std::string, std::vector<std::vector<double>>>RNASeqData;
typedef std::unordered_map<std::string, std::vector<std::string>> PatientList;
typedef std::vector<std::pair<std::string, int>> GeneList;

std::vector<std::string> getCancerNames(const std::string &filename);
void buildPatientIDs(const std::string &cancerName,
		PatientList &patientControlData, PatientList &patientTumorData);
void readPatientData(const std::string &filename,
		PatientList &patientControlData, PatientList &patientTumorData);
void readPatientData(const std::vector<std::string> &cancers,
		PatientList &patientControlData, PatientList &patientTumorData);
GeneList makeGeneMapping(const std::string &filename);

void readRNASeq(const std::string &cancerName, const std::string& dataType,
		const std::string &patientName, const GeneList &geneMapping,
		RNASeqData &rnaSeqData);

void initializeRNASeqData(const PatientList &patientData,
		unsigned int numberOfGenes, RNASeqData &rnaSeqData);

void readRNASeqData(const PatientList &patientData, const GeneList &geneMapping,
		RNASeqData &rnaSeqData, bool isTumorData, int maxPatients);
void readData(const PatientList &patientControlData,
		const PatientList &patientTumorData, const GeneList &geneMapping,
		RNASeqData &controlData, RNASeqData &tumorData, int maxPatients = 2000);

void exportDataToFile(const PatientList &patientControlData,
		const PatientList &patientTumorData, int numberOfProteins,
		const RNASeqData &controlSeqData, const RNASeqData &cancerSeqData,
		const std::string &fileName);
void importDataFromFile(PatientList &patientControlData,
		PatientList &patientTumorData, RNASeqData &controlSeqData,
		RNASeqData &cancerSeqData, const std::string &fileName);

#endif // DATAREADER_H_INCLUDED
