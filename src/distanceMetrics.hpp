#ifndef SRC_DISTANCEMETRICS_HPP_
#define SRC_DISTANCEMETRICS_HPP_

#include "typedefs.hpp"

//Compute the pearson correlation matrix (output in 1D) between a range of vectors
MatrixX computePairwisePearsonCorrelation(const MatrixX &data);

//Compute the spearman correlation matrix (output in 1D) between a range of vectors
MatrixX computePairwiseSpearmanCorrelation(const MatrixX &data);

//Compute the euclidean distance matrix between a range of vectors
MatrixX computePairwiseEuclideanDistance(const MatrixX &data);

//Compute the manhattan distance matrix between a range of vectors
MatrixX computePairwiseManhattanDistance(const MatrixX &data);

MatrixX computePairwiseCosineSimilarity(const MatrixX &data);

#endif /* SRC_DISTANCEMETRICS_HPP_ */
