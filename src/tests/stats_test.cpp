#include <iostream>
#include <vector>

#include "stats_test.hpp"
#include "../stats.hpp"
#include "../utilities.hpp"

using namespace std;

//R output is 4.666667
void meanTest(){
	vector<double> x { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	double mean = computeMean(x);
	cout << "Mean = " << mean << endl;
}

void meanVectorTest(){
	std::vector<double> x { 0,0,0 };
	std::vector<double> y { 2,0,3 };
	std::vector<double> z { 1,1,6 };
	std::vector<std::vector<double>> data{x, y, z};
	std::vector<double> mean = computeMeanVector(data);
	print_vector(mean);
}

//R output is 3.872983
void standardDeviationTest(){
	std::vector<double> x { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	double std_dev = computeStandardDeviation(x);
	std::cout << "Standard Deviation = " << std_dev << std::endl;
}

//R output is 0.335578
void pearsonCorrelationTest1(){
	std::vector<double> x { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	std::vector<double> y { 4, 0, 5, 0, 5, 3, -5, 1, 0 };
	std::cout << computePearsonCorrelation(x, y) << std::endl;
}

//Test with ties
//R output is 0.136201
//Output is 0.136201
void spearmanCorrelationTest1() {
	std::vector<double> x { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	std::vector<double> y { 4, 0, 5, 0, 5, 3, -5, 1, 0 };
	std::cout << computeSpearmanCorrelation(x, y) << std::endl;
}
