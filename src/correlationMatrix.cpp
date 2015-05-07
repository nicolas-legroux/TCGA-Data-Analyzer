#include "correlationMatrix.hpp"
#include "stats.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

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
			for (int k = 0; k < numberOfGenes; ++k) {
				patientData[k] = controlData.at(cancerName).at(k).at(j);
			}
			data.push_back(patientData);
			countPatients++;
		}

		for (unsigned int j = 0; j < tumorData.at(cancerName).at(0).size();
				++j) {

			dataTypeMapping[cancerName + "-" + "Tumor"].push_back(
					countPatients);
			dataIdentifiers.push_back(
					DataIdentifier(cancerName, true,
							patientTumorList.at(cancerName).at(j)));
			vector<double> patientData(numberOfGenes);
			for (int k = 0; k < numberOfGenes; ++k) {
				patientData[k] = tumorData.at(cancerName).at(k).at(j);
			}
			data.push_back(patientData);
			countPatients++;
		}
	}

	cout << "Done." << endl;
}

std::vector<double> pearson(std::vector<std::vector<double>> &data) {
	cout << endl << "Computing Pearson Correlation... " << flush;
	return computePearsonCorrelation(data);
}

std::vector<double> spearman(std::vector<std::vector<double>> &data) {
	cout << endl << "Computing Spearman Correlation... " << flush;
	return computeSpearmanCorrelation(data);
}

void exportCorrelationMatrix(const std::vector<double> &correlationMatrix,
		const std::vector<DataIdentifier> &dataIdentifiers,
		const std::string &filemaneMatrix, const std::string &patientsIds,
		const std::string &filenameHeatMapLabels) {

	cout << endl << "Exporting Correlation matrix...";

	ofstream matrixOutputStream("export/" + filemaneMatrix);
	ofstream patientsOutputStream("export/" + patientsIds);

	int N = dataIdentifiers.size();

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			matrixOutputStream << correlationMatrix[i * N + j] << "\t";
		}
		matrixOutputStream << endl;
	}

	string current = "";
	int countCurrent = 0;
	ofstream outputStreamLabels("export/" + filenameHeatMapLabels);

	for (const DataIdentifier &dataIdentifier : dataIdentifiers) {
		string newCurrent= dataIdentifier.cancerName + "-" + ((dataIdentifier.isTumor)? "Tumor" : "Control");
		if(newCurrent != current){
			if(countCurrent != 0){
				outputStreamLabels << current << " " << countCurrent << endl;
			}
			current = newCurrent;
			countCurrent  = 1;
		}
		else{
			countCurrent++;
		}
		patientsOutputStream << dataIdentifier.toString() << endl;
	}

	cout << " Done." << endl;
}

void exportGeneralStats(const std::vector<double> &correlationMatrix,
		const DataTypeMapping &dataTypeMapping,
		const string &filemaneCorrelationMeans, const string &filemaneClasses) {

	cout << endl;

	vector<string> classes;
	int N = (int) sqrt(correlationMatrix.size());

	for (const auto &kv : dataTypeMapping) {
		string dataType = kv.first;
		if (kv.second.size() > 0) {
			classes.push_back(dataType);
		}
	}

	sort(classes.begin(), classes.end());

	ofstream outputStreamClasses("export/" + filemaneClasses);

	for (const string &s : classes) {
		outputStreamClasses << s << " " << dataTypeMapping.at(s).size() << endl;
	}

	outputStreamClasses.close();

	int n = classes.size();
	vector<double> mean_correlation(n * n);

	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			double sum = 0;
			for (double I : dataTypeMapping.at(classes[i])) {
				for (double J : dataTypeMapping.at(classes[j])) {
					if (I != J) {
						sum += correlationMatrix[N * I + J];
					}
				}
			}
			int mi = dataTypeMapping.at(classes[i]).size();
			int mj = dataTypeMapping.at(classes[j]).size();
			if (i == j) {
				sum /= (double) (mi * (mi - 1));
			} else {
				sum /= (double) (mi * mj);
			}
			mean_correlation[n * i + j] = sum;
			mean_correlation[n * j + i] = sum;
		}
	}

	ofstream outputStream("export/" + filemaneCorrelationMeans);
	outputStream << "CLASSES";
	for (const string &s : classes) {
		outputStream << "\t" << s;
	}
	outputStream << endl;

	for (int i = 0; i < n; ++i) {
		outputStream << classes[i];
		for (int j = 0; j < n; ++j) {
			outputStream << "\t" << mean_correlation[n * i + j];
		}
		outputStream << endl;
	}
}
