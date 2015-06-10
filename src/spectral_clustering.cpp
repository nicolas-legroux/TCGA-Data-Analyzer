#include "spectral_clustering.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <Eigen/Eigenvalues>
#include "utilities.hpp"
#include "k_means.hpp"
#include "distanceMatrix.hpp"
#include "typedefs.hpp"
#include "clustering.hpp"

using std::vector;
using std::cout;
using std::endl;
using Eigen::MatrixXd;

Spectral_Clustering::Spectral_Clustering(const MatrixX &_distanceMatrix,
		const DistanceMetric &_distanceMetric,
		SimilarityGraphTransformation _similarityGraphTransformation,
		SpectralClusteringParameters _similarityGraphTransformationParameters) :
		distanceMatrix(_distanceMatrix), matrixType(
				getMatrixType(_distanceMetric)), similarityGraphTransformation(
				_similarityGraphTransformation), spectralClusteringParameters(
				_similarityGraphTransformationParameters) {
	unsigned int n = distanceMatrix.cols();
	distanceMatrixCopy.resize(n, n);
}

void Spectral_Clustering::initializeSimilarityMatrix() {
	unsigned int n = distanceMatrix.cols();

	//The original data is stored in row-major order
	//This computes W
	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j < n; ++j) {
			double dist_ij = distanceMatrix(i, j);
			if (i == j) {
				if (matrixType == MatrixType::SIMILARITY) {
					distanceMatrixCopy(i, j) = 0;
				} else if (matrixType == MatrixType::DISTANCE) {
					distanceMatrixCopy(i, j) =
							std::numeric_limits<double>::infinity();
				} else {
					assert(false);
				}
			} else {
				if (matrixType == MatrixType::SIMILARITY) {
					distanceMatrixCopy(i, j) = std::fabs(dist_ij); //Correlations can be negative
				}

				else if (matrixType == MatrixType::DISTANCE
						&& similarityGraphTransformation
								== SimilarityGraphTransformation::NO_TRANSFORMATION) {
					//Use Gaussian similarity function
					distanceMatrixCopy(i, j) =
							std::exp(
									-1.0 * dist_ij
											/ (2.0
													* std::pow(
															spectralClusteringParameters.gaussianStandardDeviation,
															2.0)));
				}

				else {
					assert(
							matrixType == MatrixType::DISTANCE
									&& similarityGraphTransformation
											== SimilarityGraphTransformation::K_NEAREST_NEIGHBORS);
					distanceMatrixCopy(i, j) = std::fabs(dist_ij);
				}
			}
		}
	}
}

void Spectral_Clustering::computeLaplacianMatrix() {
	unsigned int n = distanceMatrix.cols();
//Assume W has been computed in Matrix
//Compute Laplacian, L = D-W
	for (unsigned int i = 0; i < n; ++i) {
		distanceMatrixCopy(i, i) = 0.0;
		double sum = 0;
		for (unsigned int j = 0; j < n; ++j) {
			double d = distanceMatrixCopy(i, j);
			sum += d;
			distanceMatrixCopy(i, j) = -1.0 * d;
		}
		distanceMatrixCopy(i, i) = std::fabs(sum);
	}
}

void Spectral_Clustering::transformSimilarityMatrix() {

	unsigned int n = distanceMatrix.cols();

	MatrixXd transformedMatrix = MatrixXd::Zero(n, n);

	if (similarityGraphTransformation
			== SimilarityGraphTransformation::K_NEAREST_NEIGHBORS) {
		for (unsigned int i = 0; i < n; ++i) {
			vector<double> vec;
			for (unsigned int j = 0; j < n; ++j) {
				double d = distanceMatrixCopy(j, i);
				vec.push_back(d);
			}
			vector<size_t> sortedIndexes;

			if (matrixType == MatrixType::SIMILARITY) {
				sortedIndexes = sort_indexes_decreasing(vec);
			}

			else if (matrixType == MatrixType::DISTANCE) {
				sortedIndexes = sort_indexes_increasing(vec);
			}

			for (unsigned int k = 0;
					k < spectralClusteringParameters.numberOfNeighbors; ++k) {
				unsigned int indexNeighbor = sortedIndexes[k];
				transformedMatrix(i, indexNeighbor) = 1.0;
				transformedMatrix(indexNeighbor, i) = 1.0;
			}
		}

		distanceMatrixCopy = transformedMatrix;
	}

	else {
		assert(
				similarityGraphTransformation
						== SimilarityGraphTransformation::NO_TRANSFORMATION);
	}
}

vector<int> Spectral_Clustering::compute(unsigned int k) {
	unsigned int n = distanceMatrix.cols();
	vector<int> clusters(n, 0);

	initializeSimilarityMatrix();
	transformSimilarityMatrix();
	computeLaplacianMatrix();

	Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver;
	eigenSolver.compute(distanceMatrixCopy);

	MatrixXd eigenvectors = eigenSolver.eigenvectors();
	MatrixX dataForKMeans(k, n);

	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j < k; ++j) {
			dataForKMeans(j, i) = eigenvectors(i, j + 1);
		}
	}

	ClusteringParameters kMeansParameters;
	kMeansParameters.setKMeansParameters(k, 1000);

	K_Means K_MeansClusterer(dataForKMeans, kMeansParameters, clusters);
	K_MeansClusterer.compute();

	return clusters;
}
