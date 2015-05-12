#include <iostream>
#include "utilities.hpp"

using namespace std;

void printAdvancement(unsigned int currentCount, unsigned int totalCount){
	cout << (100*currentCount)/(totalCount) << "% \r" << flush;
}
