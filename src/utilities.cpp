#include <iostream>
#include <vector>
#include <string>

#include "utilities.hpp"

using namespace std;

void printAdvancement(unsigned int currentCount, unsigned int totalCount){
	cout << (100*currentCount)/(totalCount) << "% \r" << flush;
}

vector<string> split(const string &s, const vector<char> &delimiters){
	vector<string> strs;
	string currentString;

	for(const char &c : s){
		if(find(delimiters.cbegin(), delimiters.cend(), c) == delimiters.cend()){
			currentString.push_back(c);
		}
		else{
			strs.push_back(currentString);
			currentString.erase();
		}
	}

	strs.push_back(currentString);

	return strs;
}
