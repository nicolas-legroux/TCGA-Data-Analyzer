#include "spectral_clustering.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <Eigen/Eigenvalues>
#include "utilities.hpp"
#include "k_means.hpp"
#include "normedVectorSpace.hpp"
#include "distanceMatrix.hpp"

using std::vector;
using std::cout;
using std::endl;
using Eigen::MatrixXd;

Spectral_Clustering::Spectral_Clustering(const vector<double> &_originalData,
		const DistanceMetric &_distanceMetric,
		SimilarityGraphTransformation _similarityGraphTransformation,
		SpectralClusteringParameters _similarityGraphTransformationParameters) :
		originalData(_originalData), matrixType(getMatrixType(_distanceMetric)), similarityGraphTransformation(
				_similarityGraphTransformation), spectralClusteringParameters(
				_similarityGraphTransformationParameters) {
	unsigned int n = std::sqrt(_originalData.size());
	matrix.resize(n, n);
}

void Spectral_Clustering::initializeSimilarityMatrix() {
	unsigned int n = std::sqrt(originalData.size());

	//The original data is stored in row-major order
	//This computes W
	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j < n; ++j) {
			if (i == j) {
				if (matrixType == MatrixType::SIMILARITY) {
					matrix(i, j) = 0;
				} else if (matrixType == MatrixType::DISTANCE) {
					matrix(i, j) = std::numeric_limits<double>::infinity();
				} else {
					assert(false);
				}
			} else {
				if (matrixType == MatrixType::SIMILARITY) {
					matrix(i, j) = std::fabs(originalData[i * n + j]); //Correlations can be negative
				}

				else if (matrixType == MatrixType::DISTANCE
						&& similarityGraphTransformation
								== SimilarityGraphTransformation::NO_TRANSFORMATION) {
					//Use Gaussian similarity function
					matrix(i, j) =
							std::exp(
									-1.0 * originalData[i * n + j]
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
					matrix(i, j) = std::fabs(originalData[i * n + j]);
				}
			}
		}
	}
}

void Spectral_Clustering::computeLaplacianMatrix() {
	unsigned int n = std::sqrt(originalData.size());
//Assume W has been computed in Matrix
//Compute Laplacian, L = D-W
	for (unsigned int i = 0; i < n; ++i) {
		matrix(i, i) = 0.0;
		double sum = 0;
		for (unsigned int j = 0; j < n; ++j) {
			double d = matrix(i, j);
			sum += d;
			matrix(i, j) = -1.0 * d;
		}
		matrix(i, i) = std::fabs(sum);
	}
}

void Spectral_Clustering::transformSimilarityMatrix() {

	unsigned int n = std::sqrt(originalData.size());

	MatrixXd transformedMatrix = MatrixXd::Zero(n, n);

	if (similarityGraphTransformation
			== SimilarityGraphTransformation::K_NEAREST_NEIGHBORS) {
		for (unsigned int i = 0; i < n; ++i) {
			vector<double> vec;
			for (unsigned int j = 0; j < n; ++j) {
				double d = matrix(j, i);
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

		matrix = transformedMatrix;
	}

	else {
		assert(
				similarityGraphTransformation
						== SimilarityGraphTransformation::NO_TRANSFORMATION);
	}
}

vector<int> Spectral_Clustering::compute(unsigned int k) {
	unsigned int n = std::sqrt(originalData.size());
	vector<int> clusters(n, 0);

	initializeSimilarityMatrix();
	transformSimilarityMatrix();
	computeLaplacianMatrix();

	Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver;
	eigenSolver.compute(matrix);

	MatrixXd eigenvectors = eigenSolver.eigenvectors();

	vector<vector<double>> dataForKMeans;
	for (unsigned int i = 0; i < n; ++i) {
		vector<double> v;
		for (unsigned int j = 0; j < k; ++j) {
			double value = eigenvectors(i, j + 1);
			v.push_back(value);
		}
		dataForKMeans.push_back(v);
	}

	K_Means<vector<double>> K_MeansClusterer(dataForKMeans, clusters, k, 1000,
			EuclideanSpace<double>(k));
	K_MeansClusterer.compute();

	return clusters;
}
