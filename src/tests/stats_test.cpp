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

void pearsonCorrelationTest2(){
	int dim = 2;
	std::vector<double> x { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	std::vector<double> y { 4, 0, 5, 0, 5, 3, -5, 1, 0 };
	std::vector<std::vector<double>> M { x, y };
	std::vector<double> correlationMatrix = computePearsonCorrelation(M);
	for (int l = 0; l < dim; ++l) {
		for (int c = 0; c < dim; ++c) {
			std::cout << correlationMatrix[dim * l + c] << "\t";
		}
		std::cout << std::endl;
	}
}

//Test with ties
//R output is 0.136201
//Output is 0.136201
void spearmanCorrelationTest1() {
	std::vector<double> x { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	std::vector<double> y { 4, 0, 5, 0, 5, 3, -5, 1, 0 };
	std::cout << computeSpearmanCorrelation(x, y) << std::endl;
}

//Same output as R :
/*
 M = matrix(c(1,2,4,5,6,7,-2,9,10,4,0,5,0,5,3,-5,1,0,0,0,-5,9,-4,1,1,2,1), nrow=9, ncol=3)
 cor(M, method="spearman")
 */
/*
1	0.136201	0.33199
0.136201	1	-0.652174
0.33199	-0.652174	1
 */
void spearmanCorrelationTest2() {
	std::vector<double> x1 { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	std::vector<double> x2 { 4, 0, 5, 0, 5, 3, -5, 1, 0 };
	std::vector<double> x3 { 0, 0, -5, 9, -4, 1, 1, 2, 1 };
	std::vector<std::vector<double>> M { x1, x2, x3 };
	std::vector<double> correlationMatrix = computeSpearmanCorrelation(M);

	for (int l = 0; l < 3; ++l) {
		for (int c = 0; c < 3; ++c) {
			std::cout << correlationMatrix[3 * l + c] << "\t";
		}
		std::cout << std::endl;
	}
}
