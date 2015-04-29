#include "dataReader.hpp"

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

using namespace std;

typedef unordered_map<string, vector<vector<double>>> RNASeqData;
typedef std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> PatientData;
typedef std::vector<std::pair<std::string, int>> GeneMapping;

vector<string> getCancerNames(const string &filename){
    ifstream input(filename);
    vector<string> cancerNames;
    string cancerName;
    while(input >> cancerName){
        cancerNames.push_back(cancerName);
    }
    return cancerNames;
}

unordered_map<string, vector<string>> getPatientIDs(const string &cancerName){
    cout << "**** Processing " << cancerName << " ****" << endl;
    string patientListFilename = "data/" + cancerName + "-normalized/patient.list";
    ifstream input(patientListFilename);

    string patientId;
    unordered_map<string, vector<string>> patientIdsByType;

    while (input >> patientId){
        vector<string> strs;
        boost::split(strs, patientId, boost::is_any_of("-."));
        if(strs.size()>=4){

            if(patientIdsByType.find(strs[3]) == patientIdsByType.end()){
                patientIdsByType.insert(make_pair(strs[3], vector<string>()));
            }

            patientIdsByType[strs[3]].push_back(strs[0] + "-" + strs[1] + "-" + strs[2]);
        }
    }

    cout << "Data types : " << endl;

    for(auto itr = patientIdsByType.begin(); itr != patientIdsByType.end(); ++itr){
        cout << itr->first << "-->" << itr->second.size() << endl;
    }

    //Checking that for normal samples we have the corresponding cancer sample
    for(auto itr = patientIdsByType.begin(); itr != patientIdsByType.end(); ++itr){
        if(itr->first != "01"){
            for(const string &s : itr->second){
                if(find(patientIdsByType["01"].begin(), patientIdsByType["01"].end(), s) == patientIdsByType["01"].end()){
                    cout << "Did not find a matching tumor sample for the control sample " << s << "-" << itr->first << endl;
                }
            }
        }
    }

    cout << "*************************" << endl << endl;

    return patientIdsByType;
}

void readPatientData(const std::string &filename, PatientData &patientData){
    vector<string> cancerNames(getCancerNames("data/" + filename));
    for(const string &cancer : cancerNames){
        unordered_map<string, vector<string>> patientIds(getPatientIDs(cancer));
        patientData.insert(std::make_pair(cancer, patientIds));
    }
}

void readPatientData(const std::vector<std::string> &cancerNames, PatientData &patientData){
    for(const string &cancer : cancerNames){
        unordered_map<string, vector<string>> patientIds(getPatientIDs(cancer));
        patientData.insert(std::make_pair(cancer, patientIds));
    }
}

GeneMapping makeGeneMapping(const string &pathToFile){
    GeneMapping geneMapping;
    fstream input(pathToFile);
    string firstLine;
    getline(input, firstLine);
    string geneId;
    double score;
    while(input >> geneId >> score){
        vector<string> strs;
        boost::split(strs, geneId, boost::is_any_of("|"));
        string hgncSymbol = strs[0];
        int entrezId = atoi(strs[1].c_str());
        geneMapping.push_back(make_pair(hgncSymbol, entrezId));
    }

    return geneMapping;
}

void readRNASeq(const string &cancerName, const string& dataType, const string &patientName, const GeneMapping &geneMapping , RNASeqData &rnaSeqData){
    string filename = patientName + "-" + dataType + ".genes.normalized.results";
    string geneId;
    string firstLine;
    double score;
    fstream input("data/" + cancerName + "-normalized/" + filename);
    getline(input, firstLine);

    int i=0;

    while(input >> geneId >> score){
        //It seems that this is actually unnecessary and costs a lot (try profiling !) !
        /*
        vector<string> strs;
        boost::split(strs, geneId, boost::is_any_of("|"));
        string hgncSymbol = strs[0];
        int entrezId = atoi(strs[1].c_str());
        if(geneMapping[i] == make_pair(hgncSymbol, entrezId)){
        }
        else{
            cout << "Genes don't match !" << endl;
        }
        */

        rnaSeqData[cancerName][i].push_back(score);
        i++;
    }
}

void initializeRNASeqData(const PatientData &patientData, unsigned int numberOfGenes, RNASeqData &rnaSeqData){
    for(auto itr = patientData.begin(); itr != patientData.end(); ++itr){
        rnaSeqData.insert(make_pair(itr->first, vector<vector<double>>()));
        rnaSeqData[itr->first].resize(numberOfGenes);
    }
}

void readRNASeqData(const PatientData &patientData, const string &dataType, const GeneMapping &geneMapping, RNASeqData &rnaSeqData, int maxPatients){
    initializeRNASeqData(patientData, geneMapping.size(), rnaSeqData);

    for(const auto &pairedData : patientData){
        const string &cancerName = pairedData.first;
        cout << "---------- Reading " << cancerName << " patient files -----------" <<endl;
        int i=0;
        if (pairedData.second.find(dataType) != pairedData.second.end()){
            for(const string &patientName : (pairedData.second).at(dataType)){

                readRNASeq(cancerName, dataType, patientName, geneMapping, rnaSeqData);
                //cout << dataType << " " << i << endl;
                i++;

                if(i>=maxPatients){
                    break;
                }
            }
        }
    }
}

void readControlData(const PatientData &patientData, const GeneMapping &geneMapping, RNASeqData &controlData, int maxPatients){
    cout << "************ READING CONTROL FILES **************" << endl;
    readRNASeqData(patientData, "11", geneMapping, controlData, maxPatients);
    cout << "*************************************************" << endl;
}

void readCancerData(const PatientData &patientData, const GeneMapping &geneMapping, RNASeqData &controlData, int maxPatients){
    cout << endl << "************ READING TUMOR FILES **************" << endl;
    readRNASeqData(patientData, "01", geneMapping, controlData, maxPatients);
    cout << "*************************************************" << endl;
}

void exportDataToFile(const PatientData &patientData, int numberOfProteins, const RNASeqData &controlSeqData, const RNASeqData &cancerSeqData, const string &fileName){

    cout << endl << "-------- Exporting data to file --------" << endl;

    ofstream outputStream("export/" + fileName);
    outputStream << patientData.size() << endl;

    for(const auto &mappedPatientData : patientData){
        outputStream << mappedPatientData.first << endl;

        //Start by writing the patientData

        for(const auto &mappedPatientDataCancerType :  mappedPatientData.second){
            bool writeData = false;
            if(mappedPatientDataCancerType.first == "01"){
                outputStream << "T " << mappedPatientDataCancerType.second.size() << endl;
                writeData = true;
            }
            else if(mappedPatientDataCancerType.first == "11"){
                outputStream << "N " << mappedPatientDataCancerType.second.size() << endl;
                writeData = true;
            }

            if(writeData){
                for (const string &patientName : mappedPatientDataCancerType.second){
                    outputStream << patientName << endl;
                }
            }
        }

        if(mappedPatientData.second.find("01") == mappedPatientData.second.end()){
            outputStream << "T 0" << endl;
        }
        if(mappedPatientData.second.find("11") == mappedPatientData.second.end()){
            outputStream << "N 0" << endl;
        }
    }

    // Number of proteins
    outputStream << numberOfProteins << endl;

    //Write the control Data
    for(const auto &mappedData : controlSeqData){
        const string &cancerName = mappedData.first;
        const vector<vector<double>> &data = mappedData.second;
        cout << "Exporting " << cancerName << " control data..." << endl;
        outputStream << cancerName << " " << data[0].size() << endl;
        if(data[0].size()>0){
            for(const vector<double> &seqData : data){
                for(double d : seqData){
                    outputStream << d << " ";
                }
                outputStream << endl;
            }
        }
    }

    //Write the tumor data
    for(const auto &mappedData : cancerSeqData){
        const string &cancerName = mappedData.first;
        const vector<vector<double>> &data = mappedData.second;
        cout << "Exporting " << cancerName << " tumor data..." << endl;
        outputStream << cancerName << " " << data[0].size() << endl;
        if(data[0].size()>0){
            for(const vector<double> &seqData : data){
                for(double d : seqData){
                    outputStream << d << " ";
                }
                outputStream << endl;
            }
        }
    }

    cout << "-------- Done exporting data --------" << endl;

    outputStream.close();
}

void importSampleNames(ifstream &inputStream, PatientData &patientData, const string &cancerName, const string &dataType, int numberOfSamples){

    patientData[cancerName].insert(make_pair(dataType, vector<string>()));

    string sampleName;
    for(int i=0; i<numberOfSamples; ++i){
        inputStream >> sampleName;
        patientData[cancerName][dataType].push_back(sampleName);
    }
}

void importSeqData(ifstream &inputStream, const string& cancerName, RNASeqData &seqData, int numberOfProteins, int numberOfSamples){
    double d;
    for(int i=0; i<numberOfProteins; ++i){
        for(int j=0; j<numberOfSamples; ++j){
            inputStream >> d;
            seqData[cancerName][i].push_back(d);
        }
    }
}

void importDataFromFile(PatientData &patientData, RNASeqData &controlSeqData, RNASeqData &cancerSeqData, const string &fileName){

    cout << endl << "****** IMPORTING DATA FROM FILE ******" << endl;
    int numberCancers;
    ifstream inputStream("export/" + fileName);

    inputStream >> numberCancers;
    for(int i=0; i<numberCancers; ++i){
        string cancerName;
        string dataType;
        int numberOfSamples;
        inputStream >> cancerName;
        patientData.insert(make_pair(cancerName, unordered_map<string, vector<string>>()));
        cout << "-- Processing " << cancerName << " Sample Names --" << endl;
        inputStream >> dataType >> numberOfSamples;
        (dataType == "T") ? dataType = "01" : dataType = "11";
        importSampleNames(inputStream, patientData, cancerName, dataType, numberOfSamples);
        cout << dataType << " --> " << numberOfSamples << endl;
        inputStream >> dataType >> numberOfSamples;
        (dataType == "T") ? dataType = "01" : dataType = "11";
        importSampleNames(inputStream, patientData, cancerName, dataType, numberOfSamples);
        cout << dataType << " --> " << numberOfSamples << endl;
        cout << "--------------------------" << endl;
    }

    int numberOfProteins;
    inputStream >> numberOfProteins;
    cout << endl << "Number of proteins : "  << numberOfProteins << endl;

    ///Import control data
    for(int i=0; i<numberCancers; ++i){
        string cancerName;
        int numberOfSamples;
        inputStream >> cancerName >> numberOfSamples;
        cout << "Importing " << cancerName << " control samples (" << numberOfSamples << ")" << endl;
        controlSeqData.insert(make_pair(cancerName, vector<vector<double>>()));
        if(numberOfSamples>0){
            controlSeqData[cancerName].resize(numberOfProteins);
            importSeqData(inputStream, cancerName, controlSeqData, numberOfProteins, numberOfSamples);
        }
    }

    for(int i=0; i<numberCancers; ++i){
        string cancerName;
        int numberOfSamples;
        inputStream >> cancerName >> numberOfSamples;
        cout << "Importing " << cancerName << " tumor samples (" << numberOfSamples << ")" << endl;
        cancerSeqData.insert(make_pair(cancerName, vector<vector<double>>()));
        if(numberOfSamples>0){
            cancerSeqData[cancerName].resize(numberOfProteins);
            importSeqData(inputStream, cancerName, cancerSeqData, numberOfProteins, numberOfSamples);
        }
    }
}
