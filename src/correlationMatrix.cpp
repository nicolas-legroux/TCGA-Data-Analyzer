#include "correlationMatrix.hpp"
#include "stats.hpp"
#include <iostream>
#include <fstream>

using namespace std;

typedef std::unordered_map<std::string, std::vector<int>> DataTypeMapping;
typedef std::unordered_map<std::string, std::vector<std::vector<double>>>RNASeqData;
typedef std::unordered_map<std::string, std::vector<std::string>> PatientList;

void prepareData(std::vector<std::vector<double>> &data,
		std::vector<DataIdentifier> &dataIdentifiers,
		DataTypeMapping &dataTypeMapping, const PatientList &patientControlList,
		const PatientList &patientTumorList, const RNASeqData &controlData,
		const RNASeqData &tumorData) {

	cout << endl << "Preparing data for correlation computation... ";

	int countPatients = 0;

	for (const auto &kv : tumorData) {
		string cancerName = kv.first;
		dataTypeMapping.insert(
				make_pair(cancerName + "-" + "Tumor", vector<int>()));
		dataTypeMapping.insert(
				make_pair(cancerName + "-" + "Control", vector<int>()));

		int numberOfGenes = kv.second.size();

		for (unsigned int j = 0; j < controlData.at(cancerName).at(0).size();
				++j) {

			dataTypeMapping[cancerName + "-" + "Control"].push_back(
					countPatients);
			dataIdentifiers.push_back(
					DataIdentifier(cancerName, false,
							patientControlList.at(cancerName).at(j)));
			vector<double> patientData(numberOfGenes);
			for(int k=0; k<numberOfGenes; ++k){
				patientData[k] = controlData.at(cancerName).at(k).at(j);
			}
			data.push_back(patientData);
		}

		for (unsigned int j = 0; j < tumorData.at(cancerName).at(0).size();
				++j) {

			dataTypeMapping[cancerName + "-" + "Tumor"].push_back(
					countPatients);
			dataIdentifiers.push_back(
					DataIdentifier(cancerName, true,
							patientTumorList.at(cancerName).at(j)));
			vector<double> patientData(numberOfGenes);
			for(int k=0; k<numberOfGenes; ++k){
				patientData[k] = tumorData.at(cancerName).at(k).at(j);
			}
			data.push_back(patientData);
		}
	}

	cout << "Done." << endl;
}

std::vector<double> pearson(std::vector<std::vector<double>> &data) {
	cout << endl << "Computing Pearson Correlation... " << flush;
	return computePearsonCorrelation(data);
	cout << " Done." << endl;
}

std::vector<double> spearman(std::vector<std::vector<double>> &data) {
	cout << endl << "Computing Spearman Correlation... " << flush;
	return computeSpearmanCorrelation(data);
	cout << " Done." << endl;
}

void exportCorrelationMatrix(const std::vector<double> &correlationMatrix,
		const std::vector<DataIdentifier> &dataIdentifiers,
		const std::string &filemaneMatrix, const std::string &patientsIds){

	cout << endl << "Exporting Correlation matrix...";

	ofstream matrixOutputStream("export/" + filemaneMatrix);
	ofstream patientsOutputStream("export/" + patientsIds);

	int N = dataIdentifiers.size();

	for(int i=0; i<N; ++i){
		for(int j=0; j<N; ++j){
			matrixOutputStream << correlationMatrix[i*N+j] << "\t";
		}
		matrixOutputStream << endl;
	}

	for(const DataIdentifier &dataIdentifier : dataIdentifiers){
		patientsOutputStream << dataIdentifier.toString() << endl;
	}

	cout << " Done." << endl;
}
