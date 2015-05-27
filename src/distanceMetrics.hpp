#ifndef SRC_DISTANCEMETRICS_HPP_
#define SRC_DISTANCEMETRICS_HPP_

#include <vector>

//Compute the pearson correlation matrix (output in 1D) between a range of vectors
std::vector<double> computePairwisePearsonCorrelation(
		const std::vector<std::vector<double>> &M);

//Compute the spearman correlation matrix (output in 1D) between a range of vectors
std::vector<double> computePairwiseSpearmanCorrelation(
		const std::vector<std::vector<double>> &M);

//Compute the euclidean distance matrix between a range of vectors
std::vector<double> computePairwiseEuclideanDistance(
		const std::vector<std::vector<double>> &M);

//Compute the manhattan distance matrix between a range of vectors
std::vector<double> computePairwiseManhattanDistance(
		const std::vector<std::vector<double>> &M);

std::vector<double> computePairwiseCosineSimilarity(
		const std::vector<std::vector<double>> &M);





#endif /* SRC_DISTANCEMETRICS_HPP_ */
