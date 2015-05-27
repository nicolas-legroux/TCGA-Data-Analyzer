#include <vector>
#include <string>
#include <iostream>
#include <iterator>

#include "../utilities.hpp"
#include "utilities_test.hpp"

using namespace std;

void splitTest(){
	string test = "TCGA-A6-5661-01.genes.normalized.results";
	vector<string> strs = split(test, vector<char>{'.', '-'});
	for(const string &s : strs){
		cout << s << endl;
	}
}

void copy_if_two_ranges_test(){
	vector<double> data{1,2,3,4,5,6,7,8,9};
	vector<bool> mustCopy{true, true, false, false, true, true, true, false, false};
	vector<double> dataCopy;
	copy_if_two_ranges(data.cbegin(), data.cend(), mustCopy.cbegin(), back_inserter(dataCopy));
	for_each(dataCopy.cbegin(), dataCopy.cend(), [](double d){ cout << d << " "; });
}

void cosineSimilarityTest(){
	vector<double> x{0,1};
	vector<double> y{1,0};
	cout << cosineSimilarity(x, y) << endl;
}
