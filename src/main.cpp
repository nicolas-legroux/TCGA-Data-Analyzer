#include "tests/stats_test.hpp"
#include "tests/k_means_test.hpp"
#include "tests/correlationMatrix_test.hpp"
#include "tests/general_test.hpp"
#include "tests/lodePNG_test.hpp"
#include "tests/patientUnsupervisedNormalization_test.hpp"

using namespace std;

int main() {

	kMeansTest1(2, 1000);

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

	return 0;
}
