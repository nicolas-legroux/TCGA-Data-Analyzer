#ifndef DATAREADER_H_INCLUDED
#define DATAREADER_H_INCLUDED

#include <vector>
#include <string>
#include <unordered_map>

typedef std::unordered_map<std::string, std::vector<std::vector<double>>> RNASeqData;
typedef std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> PatientData;
typedef std::vector<std::pair<std::string, int>> GeneMapping;

std::vector<std::string> getCancerNames(const std::string &filename);
std::unordered_map<std::string, std::vector<std::string>> getPatientIDs(const std::string &cancerName);
void readPatientData(const std::string &filename, PatientData &patientData);
void readPatientData(const std::vector<std::string> &cancers, PatientData &patientData);
GeneMapping makeGeneMapping(const std::string &filename);

void readRNASeq(const std::string &cancerName, const std::string& dataType, const std::string &patientName, const GeneMapping &geneMapping , RNASeqData &rnaSeqData);

void initializeRNASeqData(const PatientData &patientData, unsigned int numberOfGenes, RNASeqData &rnaSeqData);

void readRNASeqData(const PatientData &patientData, const std::string &dataType, const GeneMapping &geneMapping, RNASeqData &cancerData, int maxPatients);
void readControlData(const PatientData &patientData, const GeneMapping &geneMapping, RNASeqData &controlData, int maxPatients=2000);
void readCancerData(const PatientData &patientData, const GeneMapping &geneMapping, RNASeqData &controlData, int maxPatients=2000);

void exportDataToFile(const PatientData &patientData, int numberOfProteins, const RNASeqData &controlSeqData, const RNASeqData &cancerSeqData, const std::string &fileName);
void importDataFromFile(PatientData &patientData, RNASeqData &controlSeqData, RNASeqData &cancerSeqData, const std::string &fileName);



#endif // DATAREADER_H_INCLUDED
