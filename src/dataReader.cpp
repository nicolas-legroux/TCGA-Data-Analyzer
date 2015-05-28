#include <fstream>
#include <iostream>
#include <unordered_map>
#include <algorithm>

#include "dataReader.hpp"
#include "utilities.hpp"

using namespace std;

/*
 *
 * FUNCTIONS TO READ DATA FROM TCGA FILES
 *
 */

vector<string> getCancerNames(const string &filename) {
	ifstream input("data/" + filename);
	vector<string> cancerNames;
	string cancerName;
	while (input >> cancerName) {
		cancerNames.push_back(cancerName);
	}
	return cancerNames;
}

void readPatientNameByCancer(const std::string &cancerName,
		PatientList &controlPatientList, PatientList &tumorPatientList) {
	cout << "**** Processing " << cancerName << " Patients ****" << endl;
	string patientListFilename = "data/" + cancerName
			+ "-normalized/patient.list";
	ifstream input(patientListFilename);

	string patientId;
	int countControl = 0;
	int countTumor = 0;

	while (input >> patientId) {
		vector<string> strs = split(patientId, vector<char> { '-', '.' });
		if (strs.size() >= 4) {
			string patientName = strs[0] + "-" + strs[1] + "-" + strs[2];
			if (strs[3] == "01") {
				tumorPatientList[cancerName].push_back(patientName);
				++countTumor;
			} else if (strs[3] == "11") {
				controlPatientList[cancerName].push_back(patientName);
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
	for (const string &s : controlPatientList[cancerName]) {
		if (find(tumorPatientList[cancerName].begin(),
				tumorPatientList[cancerName].end(), s)
				== tumorPatientList[cancerName].end()) {
			cout
					<< "Did not find a matching tumor sample for the control sample "
					<< s << endl;
		}
	}

	cout << "*************************" << endl << endl;
}

void readPatientList(const vector<string> &cancers, Data &data) {
	for (const string &cancer : cancers) {
		data.controlPatientList.insert(make_pair(cancer, vector<string>()));
		data.tumorPatientList.insert(make_pair(cancer, vector<string>()));
		readPatientNameByCancer(cancer, data.controlPatientList,
				data.tumorPatientList);
	}
}

void readGeneList(const string &cancer, const string &patient, Data &data) {
	string pathToFile = "data/" + cancer + "-normalized/" + patient
			+ "-01.genes.normalized.results";
	fstream input(pathToFile);
	string firstLine;
	getline(input, firstLine);
	string geneId;
	double score;
	while (input >> geneId >> score) {
		vector<string> strs = split(geneId, vector<char> { '|' });
		string hgncSymbol = strs[0];
		int entrezId = atoi(strs[1].c_str());
		data.geneList.push_back(make_pair(hgncSymbol, entrezId));
	}
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
			if (i >= maxPatients) {
				break;
			}

			if (isTumorData) {
				readRNASeqFromFile(cancerName, "01", patientName, geneMapping,
						rnaSeqData);
			} else {
				readRNASeqFromFile(cancerName, "11", patientName, geneMapping,
						rnaSeqData);
			}

			i++;
		}

		cout << "Found " << rnaSeqData[cancerName][0].size() << " patients."
				<< endl;
	}
}

void readRNASeqData(Data &data, int maxControl, int maxTumor) {
	readRNASeqData(data.controlPatientList, data.geneList,
			data.controlRNASeqData, false, maxControl);
	readRNASeqData(data.tumorPatientList, data.geneList, data.tumorRNASeqData,
			true, maxTumor);
}

void readData(const vector<string> &cancers, Data &data, int maxControl,
		int maxTumor) {
	//Load Patient Data
	readPatientList(cancers, data);
	//Load GeneList
	string firstCancer = (*data.tumorPatientList.begin()).first;
	string firstPatient = *(*data.tumorPatientList.begin()).second.begin();
	readGeneList(firstCancer, firstPatient, data);
	//Load RNASeq Data
	readRNASeqData(data, maxControl, maxTumor);
}

void readData(const string &filename, Data &data, int maxControl,
		int maxTumor) {
	vector<string> cancers = getCancerNames(filename);
	readData(cancers, data, maxControl, maxTumor);
}

/*
 *
 * FUNCTION TO EXPORT TCGA DATA IN ONE BIG FILE
 *
 */

void exportDataToFile(const Data &data, const string &filename) {

	cout << endl << "-------- Exporting data to file --------" << endl;

	ofstream outputStream("export/" + filename);

	//Print number of cancers
	outputStream << data.tumorPatientList.size() << endl;

	for (const auto &mappedPatientData : data.tumorPatientList) {
		string cancerName = mappedPatientData.first;
		//print cancer, number of control samples, number of tumor samples
		outputStream << cancerName << " "
				<< data.controlPatientList.at(cancerName).size() << " "
				<< data.tumorPatientList.at(cancerName).size() << endl;

		//Start by writing the control sample names, then the tumor sample names
		for (const string &s : data.controlPatientList.at(cancerName)) {
			outputStream << s << endl;
		}
		for (const string &s : data.tumorPatientList.at(cancerName)) {
			outputStream << s << endl;
		}
	}

	// Number of genes
	int numberOfGenes = (*data.tumorRNASeqData.begin()).second.size();
	outputStream << numberOfGenes << endl;

	//Write the control Data
	for (const auto &kv : data.controlRNASeqData) {
		const string &cancerName = kv.first;
		const vector<vector<double>> &rnaSeqData = kv.second;
		cout << "Exporting " << cancerName << " control data..." << endl;
		outputStream << cancerName << " " << rnaSeqData[0].size() << endl;
		if (rnaSeqData[0].size() > 0) {
			for (const vector<double> &geneData : rnaSeqData) {
				for (double d : geneData) {
					outputStream << d << " ";
				}
				outputStream << endl;
			}
		}
	}

	//Write the tumor data
	for (const auto &kv : data.tumorRNASeqData) {
		const string &cancerName = kv.first;
		const vector<vector<double>> &rnaSeqData = kv.second;
		cout << "Exporting " << cancerName << " tumor data..." << endl;
		outputStream << cancerName << " " << rnaSeqData[0].size() << endl;
		if (rnaSeqData[0].size() > 0) {
			for (const vector<double> &geneData : rnaSeqData) {
				for (double d : geneData) {
					outputStream << d << " ";
				}
				outputStream << endl;
			}
		}
	}

	cout << "-------- Done exporting data --------" << endl;

	outputStream.close();
}

/*
 *
 * FUNCTIONS TO IMPORT DATA FROM FILE (QUICKER THAN READING TCGA FILES)
 *
 */

void importSampleNames(ifstream &inputStream, PatientList &patientData,
		const string &cancerName, int numberOfSamples) {

	string sampleName;
	for (int i = 0; i < numberOfSamples; ++i) {
		inputStream >> sampleName;
		patientData[cancerName].push_back(sampleName);
	}
}

void importRNASeqData(ifstream &inputStream, const string& cancerName,
		RNASeqData &seqData, int numberOfGenes, int numberOfSamples) {
	double d;
	for (int i = 0; i < numberOfGenes; ++i) {
		for (int j = 0; j < numberOfSamples; ++j) {
			inputStream >> d;
			seqData[cancerName][i].push_back(d);
		}
	}
}

void importDataFromFile(Data &data, const string &filename) {

	cout << endl << "****** IMPORTING DATA FROM FILE ******" << endl;
	int numberCancers;
	ifstream inputStream("export/" + filename);

	//Import patient list data
	inputStream >> numberCancers;
	for (int i = 0; i < numberCancers; ++i) {
		string cancerName;
		int numberOfControlSamples;
		int numberOfTumorSamples;

		inputStream >> cancerName >> numberOfControlSamples
				>> numberOfTumorSamples;

		data.controlPatientList.insert(make_pair(cancerName, vector<string>()));
		data.tumorPatientList.insert(make_pair(cancerName, vector<string>()));

		cout << "Processing " << cancerName << " Sample Names..." << endl;
		importSampleNames(inputStream, data.controlPatientList, cancerName,
				numberOfControlSamples);
		importSampleNames(inputStream, data.tumorPatientList, cancerName,
				numberOfTumorSamples);
	}

	int numberOfGenes;
	inputStream >> numberOfGenes;
	cout << "Number of proteins : " << numberOfGenes << endl;

	//Import control RNASeq data
	for (int i = 0; i < numberCancers; ++i) {
		string cancerName;
		int numberOfSamples;
		inputStream >> cancerName >> numberOfSamples;
		cout << "Importing " << cancerName << " control samples ("
				<< numberOfSamples << ")" << endl;
		data.controlRNASeqData.insert(
				make_pair(cancerName, vector<vector<double>>()));
		if (numberOfSamples > 0) {
			data.controlRNASeqData[cancerName].resize(numberOfGenes);
			importRNASeqData(inputStream, cancerName, data.controlRNASeqData,
					numberOfGenes, numberOfSamples);
		}
	}

	//Import tumor RNASeq Data
	for (int i = 0; i < numberCancers; ++i) {
		string cancerName;
		int numberOfSamples;
		inputStream >> cancerName >> numberOfSamples;
		cout << "Importing " << cancerName << " tumor samples ("
				<< numberOfSamples << ")" << endl;
		data.tumorRNASeqData.insert(
				make_pair(cancerName, vector<vector<double>>()));
		if (numberOfSamples > 0) {
			data.tumorRNASeqData[cancerName].resize(numberOfGenes);
			importRNASeqData(inputStream, cancerName, data.tumorRNASeqData,
					numberOfGenes, numberOfSamples);
		}
	}

	//Build gene list
	string firstCancer = (*data.tumorPatientList.begin()).first;
	string firstPatient = *(*data.tumorPatientList.begin()).second.begin();
	readGeneList(firstCancer, firstPatient, data);

	cout << "**************************************" << endl;
}

/*
 *
 * FUNCTION TO EXPORT DATA TO ONE BIG TSV MATRIX
 * LINES = GENES
 * COLUMNS = PATIENTS
 *
 */

void exportToMatrix(const Data &data, const string &matrixFilename,
		const string &patientListFilename) {

	cout << endl << "****** Exporting raw matrix ******" << endl;
	ofstream matrixOutputStream("export/" + matrixFilename);
	ofstream patientsOutputStream("export/" + patientListFilename);

	int numberOfGenes = (*data.tumorRNASeqData.begin()).second.size();

	for (int i = 0; i < numberOfGenes; ++i) {
		if (i % 1000 == 0) {
			printAdvancement(i, numberOfGenes);
		}
		for (const auto &kv : data.tumorRNASeqData) {
			string cancerName = kv.first;
			for (unsigned int j = 0;
					j < data.controlRNASeqData.at(cancerName).at(0).size();
					++j) {
				if (i == 0) {
					patientsOutputStream << cancerName << "-Control ("
							<< data.controlPatientList.at(cancerName).at(j)
							<< ")" << endl;
				}
				matrixOutputStream
						<< data.controlRNASeqData.at(cancerName).at(i).at(j)
						<< "\t";
			}

			for (unsigned int j = 0;
					j < data.tumorRNASeqData.at(cancerName).at(0).size(); ++j) {
				if (i == 0) {
					patientsOutputStream << cancerName << "-Tumor ("
							<< data.tumorPatientList.at(cancerName).at(j) << ")"
							<< endl;
				}
				matrixOutputStream
						<< data.tumorRNASeqData.at(cancerName).at(i).at(j)
						<< "\t";
			}
		}
		matrixOutputStream << endl;
	}
	cout << "********* Done exporting *********" << endl;
}
