#ifndef SRC_K_MEANS_HPP_
#define SRC_K_MEANS_HPP_

#include <vector>
#include <random>
#include <algorithm>

#include "normedVectorSpace.hpp"
#include "utilities.hpp"

template<typename T, typename K = double>
class K_Means {
private:
	const std::vector<T> &data;
	std::vector<int> &clusters;
	unsigned int K_param;
	int Nmax;

	NormedVectorSpace<T, K> normedVectorSpace;

	bool verbose;

	std::vector<T> means;
	std::vector<bool> dataToCluster;

	void initializeClustersRandomly() {
		std::default_random_engine engine(std::random_device{}());
		std::uniform_int_distribution<int> distribution(0, data.size() - 1);
		for (unsigned int i = 0; i < K_param; i++) {
			while (true) {
				int randomInt = distribution(engine);
				if (dataToCluster[randomInt]) {
					T randomData(data[randomInt]);
					if (find(means.begin(), means.begin() + i, randomData)
							== (means.begin() + i)) {
						means[i] = randomData;
						break;
					}
				}
			}
		}
	}

	int findClosestClusterFromDataPoint(const T &dataPoint) {
		int closestCluster = 0;
		K closestDistance = normedVectorSpace.distance(dataPoint, means[0]);
		for (unsigned int i = 1; i < K_param; ++i) {
			K d = normedVectorSpace.distance(dataPoint, means[i]);
			if (d < closestDistance) {
				closestCluster = i;
				closestDistance = d;
			}
		}
		return closestCluster;
	}

	void recalculateMeans() {
		std::vector<int> clusterSize(K_param, 0);
		fill(means.begin(), means.end(), normedVectorSpace.zero());
		for (unsigned int i = 0; i != clusters.size(); ++i) {
			int cluster = clusters[i];
			if (cluster != -1) {
				++clusterSize[cluster];
				normedVectorSpace.addTo(means[cluster], data[i]);
			}
		}
		for (unsigned int i = 0; i < K_param; ++i) {
			normedVectorSpace.multiplyByConstant(means[i],
					1.0 / static_cast<double>(clusterSize[i]));
		}
	}

	int kMeansIteration() {
		bool numberClusterChange = 0;
		// Assign clusters
		for (unsigned int i = 0; i != data.size(); ++i) {
			int oldCluster = clusters[i];
			if (oldCluster != -1) {
				int newCluster = findClosestClusterFromDataPoint(data[i]);
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

	void assignSortedClusters(const std::vector<size_t> &clusterRanks) {
		for (unsigned int i = 0; i != clusters.size(); ++i) {
			int currentCluster = clusters[i];
			if (currentCluster != -1) {
				clusters[i] = clusterRanks.at(currentCluster);
			}
		}
	}

public:
	K_Means(const std::vector<T> &_data, std::vector<int> &_clusters,
			unsigned int _K, int _Nmax,
			const NormedVectorSpace<T, K> &_normedVectorSpace, bool _verbose = false) :
			data(_data), clusters(_clusters), K_param(_K), Nmax(_Nmax), normedVectorSpace(
					_normedVectorSpace), verbose(_verbose) {
	}

	std::vector<T> compute() {

		means = std::vector<T>(K_param, normedVectorSpace.zero());
		dataToCluster = std::vector<bool>(data.size());
		transform(clusters.cbegin(), clusters.cend(), dataToCluster.begin(),
				[](int cluster) {
					return (cluster != -1);
				});

		initializeClustersRandomly();

		// Iterate
		if(verbose){
			std::cout << "Percentage of points that switched between clusters : ";
		}

		int numberOfPoints = kMeansIteration();
		int iterations = 1;
		while (iterations < Nmax && numberOfPoints > 0) {
			++iterations;
			numberOfPoints = kMeansIteration();
			if(verbose){
				printAdvancement(numberOfPoints, data.size());
			}
		}

		// Sort the clusters
		std::vector<size_t> clusterRanks(get_rank_increasing(means));
		assignSortedClusters(clusterRanks);
		sort(means.begin(), means.end());

		if (iterations == Nmax) {
			std::cout << "K-Means did not converge." << std::endl;
		}

		return means;
	}

	void computeIteratedBinaryKMeans(int Niteration) {
		assert(K_param == 2);

		for (int i = 0; i < Niteration; ++i) {
			compute();   // 1000 should be enough
			transform(clusters.cbegin(), clusters.cend(), clusters.begin(),
					[](int cluster) {
						return (cluster == 0)? 0 : -1;
					});
		}
		transform(clusters.cbegin(), clusters.cend(), clusters.begin(),
				[](int cluster) {
					return (cluster == 0)? 0 : 1;
				});
	}
};

#endif /* SRC_K_MEANS_HstPP_ */
