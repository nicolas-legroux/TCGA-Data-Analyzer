#include <vector>
#include <iostream>
#include "distanceMetrics_test.hpp"
#include "../distanceMetrics.hpp"
#include "../typedefs.hpp"

void pearsonCorrelationTest2() {
	int dim = 9;
	MatrixX matrix(dim, 2);
	matrix << 1, 4, 2, 0, 4, 5, 5, 0, 6, 5, 7, 3, -2, -5, 9, 1, 10, 0;

	std::cout << "Computing pearson correlation matrix of :" << std::endl
			<< matrix << std::endl << "Result is : ";

	MatrixX correlationMatrix = computePairwisePearsonCorrelation(matrix);
	std::cout << correlationMatrix;
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
	int dim = 9;
	MatrixX matrix(dim, 9);
	matrix << 1, 4, 0, 2, 0, 0, 4, 5, -5, 5, 0, 9, 6, 5, -4, 7, 3, 1, -2, -5, 1, 9, 1, 2, 10, 0, 1;
	std::cout << "Computing spearman correlation matrix of :" << std::endl
			<< matrix << std::endl << "Result is : ";

	MatrixX correlationMatrix = computePairwiseSpearmanCorrelation(matrix);
	std::cout << correlationMatrix;
}
