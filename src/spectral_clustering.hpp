#ifndef SRC_SPECTRAL_CLUSTERING_HPP_
#define SRC_SPECTRAL_CLUSTERING_HPP_

#include <Eigen/Dense>
#include "distanceMatrix.hpp"
#include "clustering.hpp"

enum SimilarityGraphTransformation {
	NO_TRANSFORMATION, K_NEAREST_NEIGHBORS
};

struct SpectralClusteringParameters {
	//Use this paramater for K Nearest neighbor computation
	unsigned int numberOfNeighbors;

	//Use this parameter for the Gaussian Similarity function
	double gaussianStandardDeviation;


	SpectralClusteringParameters() {
		numberOfNeighbors = 0;
		gaussianStandardDeviation = 1.0;
	}
	void setKNearestNeighborsParameters(unsigned int k) {
		numberOfNeighbors = k;
	}

	void setGussianSimilarityParameters(double standard_dev){
		gaussianStandardDeviation = standard_dev;
	}
};

class Spectral_Clustering {
private:
	const MatrixX &distanceMatrix;
	MatrixX distanceMatrixCopy;
	MatrixType matrixType;
	SimilarityGraphTransformation similarityGraphTransformation;
	SpectralClusteringParameters spectralClusteringParameters;

	// Utility functions
	void initializeSimilarityMatrix();
	void computeLaplacianMatrix();
	void computeEigenvectorsOfLaplacian(unsigned int k);
	void transformSimilarityMatrix();

public:
	Spectral_Clustering(const MatrixX &_distanceMatrix,
			const DistanceMetric &_distanceMetric,
			SimilarityGraphTransformation _similarityGraphTransformation =
					SimilarityGraphTransformation::NO_TRANSFORMATION,
			SpectralClusteringParameters _similarityGraphTransformationParameters =
					SpectralClusteringParameters());
	std::vector<int> compute(unsigned int k);
};

#endif /* SRC_SPECTRAL_CLUSTERING_HPP_ */
