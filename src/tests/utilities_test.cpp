#include <vector>
#include <set>
#include <map>
#include <string>
#include <iostream>
#include <iterator>

#include "../utilities.hpp"
#include "utilities_test.hpp"

using namespace std;

void splitTest() {
	string test = "TCGA-A6-5661-01.genes.normalized.results";
	vector<string> strs = split(test, vector<char> { '.', '-' });
	for (const string &s : strs) {
		cout << s << endl;
	}
}

void copy_if_two_ranges_test() {
	vector<double> data { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	vector<bool> mustCopy { true, true, false, false, true, true, true, false,
			false };
	vector<double> dataCopy;
	copy_if_two_ranges(data.cbegin(), data.cend(), mustCopy.cbegin(),
			back_inserter(dataCopy));
	for_each(dataCopy.cbegin(), dataCopy.cend(),
			[](double d) {cout << d << " ";});
}

void cosineSimilarityTest() {
	vector<double> x { 1, 1 };
	vector<double> y { 1, 0 };
	cout << cosineSimilarity(x, y) << endl;
}

void buildIndexMap_test() {
	set<int> s { 5, 8, 10, -4, 3 };
	map<int, unsigned int> map = buildIndexMap<int>(s.begin(), s.end());
	print_map(map);
}
