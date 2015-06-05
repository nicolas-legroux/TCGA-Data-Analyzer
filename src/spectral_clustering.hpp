#ifndef SRC_SPECTRAL_CLUSTERING_HPP_
#define SRC_SPECTRAL_CLUSTERING_HPP_

#include <Eigen/Dense>

enum SimilarityGraphTransformation {
	NO_TRANSFORMATION, K_NEAREST_NEIGHBORS
};

struct SimilarityGraphTransformationParameters {
	unsigned int numberOfNeighbors;
	SimilarityGraphTransformationParameters() {
		numberOfNeighbors = 0;
	}
	void setKNearestNeighborsParameters(unsigned int k) {
		numberOfNeighbors = k;
	}
};

class Spectral_Clustering {
private:
	const std::vector<double> &originalData;
	Eigen::MatrixXd matrix;
	SimilarityGraphTransformation similarityGraphTransformation;
	SimilarityGraphTransformationParameters similarityGraphTransformationParameters;

	// Utility functions
	void initializeSimilarityMatrix();
	void computeLaplacianMatrix();
	void computeEigenvectorsOfLaplacian(unsigned int k);
	void transformSimilarityGraph();

public:
	Spectral_Clustering(const std::vector<double> &originalData,
			SimilarityGraphTransformation _similarityGraphTransformation =
					SimilarityGraphTransformation::NO_TRANSFORMATION,
			SimilarityGraphTransformationParameters _similarityGraphTransformationParameters =
					SimilarityGraphTransformationParameters());
	std::vector<int> compute(unsigned int k);
};

#endif /* SRC_SPECTRAL_CLUSTERING_HPP_ */
