#include "dataReader.hpp"

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <unordered_map>
#include <algorithm>

using namespace std;

vector<string> getCancerNames(const string &filename) {
	ifstream input(filename);
	vector<string> cancerNames;
	string cancerName;
	while (input >> cancerName) {
		cancerNames.push_back(cancerName);
	}
	return cancerNames;
}

void buildPatientIDsFromFile(const std::string &cancerName,
		PatientList &patientControlData, PatientList &patientTumorData) {
	cout << "**** Processing " << cancerName << " Patients ****" << endl;
	string patientListFilename = "data/" + cancerName
			+ "-normalized/patient.list";
	ifstream input(patientListFilename);

	string patientId;
	int countControl = 0;
	int countTumor = 0;

	while (input >> patientId) {
		vector<string> strs;
		boost::split(strs, patientId, boost::is_any_of("-."));
		if (strs.size() >= 4) {

			string patientName = strs[0] + "-" + strs[1] + "-" + strs[2];

			if (strs[3] == "01") {
				patientTumorData[cancerName].push_back(patientName);
				++countTumor;
			} else if (strs[3] == "11") {
				patientControlData[cancerName].push_back(patientName);
				++countControl;
			} else {
				cout << "Did not process " << patientId << endl;
			}
		}
	}

	cout << "Data types : " << endl;
	cout << "01 (Tumor) --> " << countTumor << endl;
	cout << "11 (Control) --> " << countControl << endl;

	//Checking that for normal samples we have the corresponding cancer sample

	for (const string &s : patientControlData[cancerName]) {
		if (find(patientTumorData[cancerName].begin(),
				patientTumorData[cancerName].end(), s)
				== patientTumorData[cancerName].end()) {
			cout
					<< "Did not find a matching tumor sample for the control sample "
					<< s << endl;
		}
	}

	cout << "*************************" << endl << endl;
}

void readPatientData(const std::string &filename,
		PatientList &patientControlData, PatientList &patientTumorData) {
	vector<string> cancerNames(getCancerNames("data/" + filename));
	for (const string &cancer : cancerNames) {
		patientControlData.insert(make_pair(cancer, vector<string>()));
		patientTumorData.insert(make_pair(cancer, vector<string>()));
		buildPatientIDsFromFile(cancer, patientControlData, patientTumorData);
	}
}

void readPatientData(const vector<string> &cancerNames,
		PatientList &patientControlData, PatientList &patientTumorData) {
	for (const string &cancer : cancerNames) {
		patientControlData.insert(make_pair(cancer, vector<string>()));
		patientTumorData.insert(make_pair(cancer, vector<string>()));
		buildPatientIDsFromFile(cancer, patientControlData, patientTumorData);
	}
}

GeneList makeGeneMapping(const string &pathToFile) {
	GeneList geneMapping;
	fstream input(pathToFile);
	string firstLine;
	getline(input, firstLine);
	string geneId;
	double score;
	while (input >> geneId >> score) {
		vector<string> strs;
		boost::split(strs, geneId, boost::is_any_of("|"));
		string hgncSymbol = strs[0];
		int entrezId = atoi(strs[1].c_str());
		geneMapping.push_back(make_pair(hgncSymbol, entrezId));
	}

	return geneMapping;
}

void readRNASeqFromFile(const string &cancerName, const string& dataType,
		const string &patientName, const GeneList &geneMapping,
		RNASeqData &rnaSeqData) {
	string filename = patientName + "-" + dataType
			+ ".genes.normalized.results";

	string filePath = "data/" + cancerName + "-normalized/" + filename;
	string geneId;
	string firstLine;
	double score;
	fstream input(filePath);
	getline(input, firstLine);

	int i = 0;

	while (input >> geneId >> score) {
		rnaSeqData[cancerName][i].push_back(score);
		i++;
	}
}

void initializeRNASeqData(const PatientList &patientData,
		unsigned int numberOfGenes, RNASeqData &rnaSeqData) {
	for (auto itr = patientData.begin(); itr != patientData.end(); ++itr) {
		rnaSeqData.insert(make_pair(itr->first, vector<vector<double>>()));
		rnaSeqData[itr->first].resize(numberOfGenes);
	}
}

void readRNASeqData(const PatientList &patientData, const GeneList &geneMapping,
		RNASeqData &rnaSeqData, bool isTumorData, int maxPatients) {
	initializeRNASeqData(patientData, geneMapping.size(), rnaSeqData);

	for (const auto &pairedData : patientData) {
		const string &cancerName = pairedData.first;
		string type;
		(isTumorData) ? type = "Tumor" : type = "Control";
		cout << "Reading " << cancerName << " " << type << " patient files... "
				<< flush;
		int i = 0;

		for (const string &patientName : pairedData.second) {
			if (isTumorData) {
				readRNASeqFromFile(cancerName, "01", patientName, geneMapping,
						rnaSeqData);
			} else {
				readRNASeqFromFile(cancerName, "11", patientName, geneMapping,
						rnaSeqData);
			}

			i++;

			if (i >= maxPatients) {
				break;
			}
		}

		cout << "Found " << rnaSeqData[cancerName][0].size() << " patients."
				<< endl;
	}
}

void readRNASeqData(const PatientList &patientControlData,
		const PatientList &patientTumorData, const GeneList &geneMapping,
		RNASeqData &controlData, RNASeqData &tumorData, int maxPatients) {
	readRNASeqData(patientControlData, geneMapping, controlData, false,
			maxPatients);
	readRNASeqData(patientTumorData, geneMapping, tumorData, true, maxPatients);
}

void exportDataToFile(const PatientList &patientControlData,
		const PatientList &patientTumorData, int numberOfProteins,
		const RNASeqData &controlSeqData, const RNASeqData &cancerSeqData,
		const string &fileName) {

	cout << endl << "-------- Exporting data to file --------" << endl;

	ofstream outputStream("export/" + fileName);

	//Print number of cancers
	outputStream << patientTumorData.size() << endl;

	for (const auto &mappedPatientData : patientTumorData) {
		string cancerName = mappedPatientData.first;
		//print cancer, number of control samples, number of tumor samples
		outputStream << cancerName << " "
				<< patientControlData.at(cancerName).size() << " "
				<< patientTumorData.at(cancerName).size() << endl;

		//Start by writing the control sample names, then the tumor sample names
		for (const string &s : patientControlData.at(cancerName)) {
			outputStream << s << endl;
		}
		for (const string &s : patientTumorData.at(cancerName)) {
			outputStream << s << endl;
		}
	}

	// Number of proteins
	outputStream << numberOfProteins << endl;

	//Write the control Data
	for (const auto &mappedData : controlSeqData) {
		const string &cancerName = mappedData.first;
		const vector<vector<double>> &data = mappedData.second;
		cout << "Exporting " << cancerName << " control data..." << endl;
		outputStream << cancerName << " " << data[0].size() << endl;
		if (data[0].size() > 0) {
			for (const vector<double> &seqData : data) {
				for (double d : seqData) {
					outputStream << d << " ";
				}
				outputStream << endl;
			}
		}
	}

	//Write the tumor data
	for (const auto &mappedData : cancerSeqData) {
		const string &cancerName = mappedData.first;
		const vector<vector<double>> &data = mappedData.second;
		cout << "Exporting " << cancerName << " tumor data..." << endl;
		outputStream << cancerName << " " << data[0].size() << endl;
		if (data[0].size() > 0) {
			for (const vector<double> &seqData : data) {
				for (double d : seqData) {
					outputStream << d << " ";
				}
				outputStream << endl;
			}
		}
	}

	cout << "-------- Done exporting data --------" << endl;

	outputStream.close();
}

void importSampleNames(ifstream &inputStream, PatientList &patientData,
		const string &cancerName, int numberOfSamples) {

	string sampleName;
	for (int i = 0; i < numberOfSamples; ++i) {
		inputStream >> sampleName;
		patientData[cancerName].push_back(sampleName);
	}
}

void importSeqData(ifstream &inputStream, const string& cancerName,
		RNASeqData &seqData, int numberOfProteins, int numberOfSamples) {
	double d;
	for (int i = 0; i < numberOfProteins; ++i) {
		for (int j = 0; j < numberOfSamples; ++j) {
			inputStream >> d;
			seqData[cancerName][i].push_back(d);
		}
	}
}

void importDataFromFile(PatientList &patientControlData,
		PatientList &patientTumorData, RNASeqData &controlSeqData,
		RNASeqData &cancerSeqData, const string &fileName) {

	cout << endl << "****** IMPORTING DATA FROM FILE ******" << endl;
	int numberCancers;
	ifstream inputStream("export/" + fileName);

	inputStream >> numberCancers;
	for (int i = 0; i < numberCancers; ++i) {
		string cancerName;
		int numberOfControlSamples;
		int numberOfTumorSamples;

		inputStream >> cancerName >> numberOfControlSamples
				>> numberOfTumorSamples;

		patientControlData.insert(make_pair(cancerName, vector<string>()));
		patientTumorData.insert(make_pair(cancerName, vector<string>()));

		cout << "Processing " << cancerName << " Sample Names..." << endl;
		importSampleNames(inputStream, patientControlData, cancerName,
				numberOfControlSamples);
		importSampleNames(inputStream, patientTumorData, cancerName,
				numberOfTumorSamples);
	}

	int numberOfProteins;
	inputStream >> numberOfProteins;
	cout << "Number of proteins : " << numberOfProteins << endl;

	///Import control data
	for (int i = 0; i < numberCancers; ++i) {
		string cancerName;
		int numberOfSamples;
		inputStream >> cancerName >> numberOfSamples;
		cout << "Importing " << cancerName << " control samples ("
				<< numberOfSamples << ")" << endl;
		controlSeqData.insert(make_pair(cancerName, vector<vector<double>>()));
		if (numberOfSamples > 0) {
			controlSeqData[cancerName].resize(numberOfProteins);
			importSeqData(inputStream, cancerName, controlSeqData,
					numberOfProteins, numberOfSamples);
		}
	}

	for (int i = 0; i < numberCancers; ++i) {
		string cancerName;
		int numberOfSamples;
		inputStream >> cancerName >> numberOfSamples;
		cout << "Importing " << cancerName << " tumor samples ("
				<< numberOfSamples << ")" << endl;
		cancerSeqData.insert(make_pair(cancerName, vector<vector<double>>()));
		if (numberOfSamples > 0) {
			cancerSeqData[cancerName].resize(numberOfProteins);
			importSeqData(inputStream, cancerName, cancerSeqData,
					numberOfProteins, numberOfSamples);
		}
	}

	cout << "**************************************" << endl;
}

void exportToMatrix(const PatientList &patientControlData,
		const PatientList &patientTumorData, const RNASeqData &controlSeqData,
		const RNASeqData &tumorSeqData, const string &matrixFileName,
		const string &patientListFileName, int numberOfGenes) {

	cout << endl << "****** Exporting raw matrix ******" << endl;

	ofstream matrixOutputStream("export/" + matrixFileName);
	ofstream patientsOutputStream("export/" + patientListFileName);

	for (int i = 0; i < numberOfGenes; ++i) {

		if (i%1000 == 0){
			cout << (100*i)/numberOfGenes << "%..." << endl;
		}

		for (const auto &kv : tumorSeqData) {
			string cancerName = kv.first;

			for (unsigned int j = 0; j < controlSeqData.at(cancerName).at(0).size();
					++j) {
				if(i == 0){
					patientsOutputStream << cancerName << "-Control ("
							<< patientControlData.at(cancerName).at(j) << ")" << endl;
				}

				matrixOutputStream << controlSeqData.at(cancerName).at(i).at(j) << "\t";
			}

			for (unsigned int j = 0; j < tumorSeqData.at(cancerName).at(0).size();
					++j) {
				if(i == 0){
					patientsOutputStream << cancerName << "-Tumor ("
							<< patientTumorData.at(cancerName).at(j) << ")" << endl;
				}

				matrixOutputStream << tumorSeqData.at(cancerName).at(i).at(j) << "\t";
			}
		}

		matrixOutputStream << endl;
	}

	cout << "********* Done exporting *********" << endl;
}
