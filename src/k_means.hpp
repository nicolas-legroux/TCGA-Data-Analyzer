#ifndef SRC_K_MEANS_HPP_
#define SRC_K_MEANS_HPP_

#include <vector>
#include <random>
#include <algorithm>
#include <cassert>

#include "utilities.hpp"
#include "typedefs.hpp"
#include "clustering.hpp"
#include <Eigen/Dense>

class K_Means {
private:
	const MatrixX &data;
	ClusteringParameters clusteringParameters;
	std::vector<int> clusters;

	MatrixX medoids;
	std::vector<bool> dataToCluster;

	void initializeClustersRandomly() {
		std::default_random_engine engine(std::random_device { }());
		std::uniform_int_distribution<int> distribution(0, data.cols() - 1);
		std::vector<VectorX> initialMedoids(clusteringParameters.K);
		for (unsigned int i = 0; i < clusteringParameters.K; ++i) {
			while (true) {
				int randomInt = distribution(engine);
				if (dataToCluster[randomInt]) {
					VectorX randomPoint = data.col(randomInt);
					if (find(initialMedoids.begin(), initialMedoids.begin() + i,
							randomPoint) == (initialMedoids.begin() + i)) {
						initialMedoids[i] = randomPoint;
						break;
					}
				}
			}
		}

		for (unsigned int i = 0; i < clusteringParameters.K; ++i) {
			medoids.col(i) = initialMedoids[i];
		}
	}

	int findClosestClusterFromDataPoint(const VectorX &dataPoint) {
		int closestCluster = 0;
		double closestDistance = (dataPoint - medoids.col(0)).squaredNorm();
		for (unsigned int i = 1; i < clusteringParameters.K; ++i) {
			double distance = (dataPoint - medoids.col(i)).squaredNorm();
			if (distance < closestDistance) {
				closestCluster = i;
				closestDistance = distance;
			}
		}
		return closestCluster;
	}

	void recalculateMeans() {
		std::vector<int> clusterSize(clusteringParameters.K, 0);
		unsigned int n = data.rows();
		for (unsigned int i = 0; i < clusteringParameters.K; ++i) {
			medoids.col(i) = VectorX::Zero(n);
		}
		for (unsigned int i = 0; i != clusters.size(); ++i) {
			int cluster = clusters[i];
			if (cluster != -1) {
				++clusterSize[cluster];
				medoids.col(cluster) += data.col(i);
			}
		}

		for (unsigned int i = 0; i < clusteringParameters.K; ++i) {
			medoids.col(i) /= static_cast<double>(clusterSize[i]);
		}
	}

	int kMeansIteration() {
		int numberClusterChange = 0;
		// Assign clusters
		for (unsigned int i = 0; i != data.cols(); ++i) {
			int oldCluster = clusters[i];
			if (oldCluster != -1) {
				int newCluster = findClosestClusterFromDataPoint(data.col(i));
				if (newCluster != clusters[i]) {
					numberClusterChange++;
					clusters[i] = newCluster;
				}
			}
		}
		// Recalculate means
		recalculateMeans();
		return numberClusterChange;
	}

public:
	K_Means(const MatrixX &_data, const ClusteringParameters &parameters,
			const std::vector<int> &_clusters) :
			data(_data), clusteringParameters(parameters), clusters(_clusters) {
		unsigned int dim = data.rows();
		medoids.resize(dim, clusteringParameters.K);
		dataToCluster = std::vector<bool>(data.cols());
		transform(clusters.cbegin(), clusters.cend(), dataToCluster.begin(),
				[](int cluster) {
					return (cluster >= 0);
				});
	}

	K_Means(const MatrixX &_data, ClusteringParameters parameters) :
			K_Means(_data, parameters, std::vector<int>(_data.cols(), 0)) {
	}

	std::vector<int> compute() {

		initializeClustersRandomly();

		// Iterate
		if (clusteringParameters.verbose) {
			std::cout
					<< "Percentage of points that switched between clusters : "
					<< std::endl;
		}

		int numberOfPoints = kMeansIteration();
		unsigned int iterations = 1;
		while (iterations < clusteringParameters.maxIterations
				&& numberOfPoints > 0) {
			++iterations;
			numberOfPoints = kMeansIteration();
			if (clusteringParameters.verbose) {
				printAdvancement(numberOfPoints, data.cols());
			}
		}

		if (iterations == clusteringParameters.maxIterations) {
			std::cout << "K-Means did not converge." << std::endl;
		}

		return clusters;
	}

	MatrixX getMedoids() {
		return medoids;
	}
};

#endif /* SRC_K_MEANS_HstPP_ */
