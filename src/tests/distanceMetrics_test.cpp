#include <vector>
#include <iostream>
#include "distanceMetrics_test.hpp"
#include "../distanceMetrics.hpp"

void pearsonCorrelationTest2(){
	int dim = 2;
	std::vector<double> x { 1, 2, 4, 5, 6, 7, -2, 9, 10 };
	std::vector<double> y { 4, 0, 5, 0, 5, 3, -5, 1, 0 };
	std::vector<std::vector<double>> M { x, y };
	std::vector<double> correlationMatrix = computePairwisePearsonCorrelation(M);
	for (int l = 0; l < dim; ++l) {
		for (int c = 0; c < dim; ++c) {
			std::cout << correlationMatrix[dim * l + c] << "\t";
		}
		std::cout << std::endl;
	}
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
	std::vector<double> correlationMatrix = computePairwiseSpearmanCorrelation(M);

	for (int l = 0; l < 3; ++l) {
		for (int c = 0; c < 3; ++c) {
			std::cout << correlationMatrix[3 * l + c] << "\t";
		}
		std::cout << std::endl;
	}
}
