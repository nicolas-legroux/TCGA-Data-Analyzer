#ifndef SRC_K_MEANS_HPP_
#define SRC_K_MEANS_HPP_

#include <vector>

typedef double (*DistanceFunction)(const std::vector<double>&, const std::vector<double>&);

/*
 *
 * ONE - DIMENSIONAL
 *
 */

std::vector<double> computeKMeans(const std::vector<double> &data,
		std::vector<int> &clusters, int K, int Nmax);
void iteratedBinaryKMeans(const std::vector<double> &data,
		std::vector<int> &clusters, int N_iter);

/*
 *
 * MULTI-DIMENSIONAL
 *
 */

std::vector<std::vector<double>> computeKMeans(const std::vector<std::vector<double>> &data,
		std::vector<int> &clusters, int K, int Nmax, DistanceFunction distance);

/*
 *
 * DISTANCE MEASURES
 *
 */

double euclidianNorm(const std::vector<double> &a, const std::vector<double> &b);
double norm1(const std::vector<double> &a, const std::vector<double> &b);

#endif /* SRC_K_MEANS_HPP_ */
