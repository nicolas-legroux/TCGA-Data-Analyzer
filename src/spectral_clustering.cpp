#include "spectral_clustering.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <Eigen/Eigenvalues>
#include "utilities.hpp"
#include "k_means.hpp"
#include "normedVectorSpace.hpp"

using std::vector;
using std::cout;
using std::endl;
using Eigen::MatrixXd;

Spectral_Clustering::Spectral_Clustering(const vector<double> &_originalData,
		SimilarityGraphTransformation _similarityGraphTransformation,
		SimilarityGraphTransformationParameters _similarityGraphTransformationParameters) :
		originalData(_originalData), similarityGraphTransformation(
				_similarityGraphTransformation), similarityGraphTransformationParameters(
				_similarityGraphTransformationParameters) {
	unsigned int n = std::sqrt(_originalData.size());
	matrix.resize(n, n);
}

void Spectral_Clustering::initializeSimilarityMatrix() {
	unsigned int n = std::sqrt(originalData.size());

	//The original data is stored in row-major order
	//This computes -W
	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j < n; ++j) {
			if (i == j) {
				matrix(i, j) = 0;
			} else {
				matrix(i, j) = std::fabs(originalData[i * n + j]);
			}
		}
	}
}

void Spectral_Clustering::computeLaplacianMatrix() {
	unsigned int n = std::sqrt(originalData.size());
	//Assume W has been computed in Matrix
	//Compute Laplacian, L = D-W
	for (unsigned int i = 0; i < n; ++i) {
		double sum = 0;
		for (unsigned int j = 0; j < n; ++j) {
			double d = matrix(i, j);
			sum += d;
			matrix(i, j) = -1.0 * d;
		}
		matrix(i, i) = std::fabs(sum);
	}
}

void Spectral_Clustering::transformSimilarityGraph() {
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
			vector<size_t> sortedIndexes = sort_indexes_decreasing(vec);
			for (unsigned int k = 0;
					k
							< similarityGraphTransformationParameters.numberOfNeighbors;
					++k) {
				unsigned int indexNeighbor = sortedIndexes[k];
				transformedMatrix(i, indexNeighbor) = 1.0;
				transformedMatrix(indexNeighbor, i) = 1.0;
			}
		}

		matrix = transformedMatrix;
		return;
	}

	assert(
			similarityGraphTransformation
					== SimilarityGraphTransformation::NO_TRANSFORMATION);
}

vector<int> Spectral_Clustering::compute(unsigned int k) {
	unsigned int n = std::sqrt(originalData.size());
	vector<int> clusters(n, 0);

	initializeSimilarityMatrix();
	transformSimilarityGraph();
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
