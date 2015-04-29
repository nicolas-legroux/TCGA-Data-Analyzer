#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include "dataReader.hpp"
#include "stats.hpp"

using namespace std;

typedef std::unordered_map<std::string, std::vector<std::vector<double>>> RNASeqData;
typedef unordered_map<string, vector<string>> PatientList;
typedef std::vector<std::pair<std::string, int>> GeneList;


int main()
{

    vector<string> cancers{"BRCA"};
    PatientList patientControlData;
    PatientList patientTumorData;
    RNASeqData controlData;
    RNASeqData tumorData;
    GeneList geneMapping(makeGeneMapping("data/BRCA-normalized/TCGA-A1-A0SJ-01.genes.normalized.results"));

    /*
    readPatientData(cancers, patientControlData, patientTumorData);
    readData(patientControlData, patientTumorData, geneMapping, controlData, tumorData, 5);
	*/
    importDataFromFile(patientControlData, patientTumorData, controlData, tumorData, "brca.export");
    exportDataToFile(patientControlData, patientTumorData, geneMapping.size(), controlData, tumorData, "brca.export2");

    /*


    unordered_map<string, vector<pair<double, double>>> test = computeControlDistribution(controlData);
    computeZScore(tumorData, test);

    vector<double> sumZScores;
    for(unsigned int i=0; i<geneMapping.size(); ++i){
        sumZScores.push_back(accumulate(tumorData["BRCA"][i].begin(), tumorData["BRCA"][i].end(), 0.0)/(double)tumorData["BRCA"][0].size());
    }

    cout << geneMapping[5678].first << " --> " << sumZScores[5678]<< endl;



    vector<size_t> sortedIndexes = sort_indexes(sumZScores);

    for(int i=0; i<100; ++i){
        cout << sortedIndexes[i] << " " << geneMapping[sortedIndexes[i]].first << " " << sumZScores[sortedIndexes[i]] << endl;
    }
    */


    /*
    PatientData patientData;
    GeneMapping geneMapping(makeGeneMapping("Data/BRCA-normalized/TCGA-A1-A0SJ-01.genes.normalized.results"));
    RNASeqData controlData;
    RNASeqData tumorData;
    importDataFromFile(patientData, controlData, tumorData, "data.export");
    */

    /*
    unsigned int gene = 13065;

   // cout << endl << "---------- Showing stats for gene " << geneMapping[gene].first << " ------------" << endl;
   // cout << "Cancer\t\t#Samples\t\tMean\t\tStandard deviation\t\t%0" << endl;

    for(const auto &pairedData : controlData){
        const string &cancerName = pairedData.first;
        int count0 = 0;
        cout << "Cancer " << cancerName << ": ";
        for(int i=0; i<geneMapping.size(); i++){
            if(computeZeroPercentage(pairedData.second.at(i))>0.80){
                count0++;
            }
        }

        cout << count0 << " genes with more than 80% zero values" << endl;
    }
    */

    /*
    vector<double> brca1;
    vector<double> brca2;
    vector<double> thca1;
    vector<double> thca2;

    for(int i=0; i < geneMapping.size(); i++){
        brca1.push_back(tumorData["KIRC"][i][0]);
        brca2.push_back(tumorData["KIRC"][i][1]);
        thca1.push_back(tumorData["THCA"][i][1]);
        thca2.push_back(tumorData["THCA"][i][4]);
    }

    cout << computeCorrelation(brca1, brca2) << endl;
    cout << computeCorrelation(brca1, thca2) << endl;
    cout << computeCorrelation(brca1, thca1) << endl;
    cout << computeCorrelation(brca2, thca1) << endl;
    cout << computeCorrelation(brca2, thca2) << endl;
    cout << computeCorrelation(thca1, thca2) << endl;
    */
    return 0;
}
