#include <vector>
#include <cmath>

#include "distanceMetrics.hpp"
#include "utilities.hpp"
#include "stats.hpp"

using namespace std;

vector<double> computePairwisePearsonCorrelation(
		const vector<vector<double>> &M) {
	unsigned int N = M.size();
	cout << endl << "Pearson correlation to be computed for " << N
			<< " vectors..." << endl;
	vector<double> correlationMatrix(N * N);
	vector<double> means(N);
	vector<double> standard_deviations(N);

	cout << "Computing all means and standard deviations... " << flush;
	transform(M.cbegin(), M.cend(), means.begin(), computeMean);
	transform(M.cbegin(), M.cend(), standard_deviations.begin(),
			[](const vector<double> &vec) {
				return computeStandardDeviation(vec, false);
			});
	cout << "Done." << endl;

	unsigned int count = 0;

	cout << "Computing correlations for all pairs of vectors... " << endl;
	for (unsigned int i = 0; i < N; ++i) {
		printAdvancement(count, (N * (N - 1) / 2));
		correlationMatrix[i + N * i] = 1.0;
		for (unsigned int j = i + 1; j < N; ++j) {
			double cor = computePearsonCorrelation(M[i], M[j], means[i],
					standard_deviations[i], means[j], standard_deviations[j]);
			correlationMatrix[i + N * j] = cor;
			correlationMatrix[j + N * i] = cor;
		}
		count += N - i + 1;
	}
	cout << "Done. " << endl << flush;

	return correlationMatrix;
}

vector<double> computePairwiseSpearmanCorrelation(
		const vector<vector<double>> &M) {
	vector<vector<double>> M_copy(M);
	for_each(M_copy.begin(), M_copy.end(), computeRank);
	return computePairwisePearsonCorrelation(M_copy);
}

vector<double> computePairwiseEuclideanDistance(
		const vector<vector<double>> &M) {
	unsigned int N = M.size();
	cout << endl << "Pairwise Euclidean Distance to be computed for " << N
			<< " vectors..." << endl;
	vector<double> distanceMatrix(N * N);

	unsigned int count = 0;

	cout << "Computing distances for all pairs of vectors... " << endl;
	for (unsigned int i = 0; i < N; ++i) {
		printAdvancement(count, (N * (N - 1) / 2));
		distanceMatrix[i + N * i] = 0.0;
		for (unsigned int j = i + 1; j < N; ++j) {
			double distance = euclideanDistance(M[i], M[j]);
			distanceMatrix[i + N * j] = distance;
			distanceMatrix[j + N * i] = distance;
		}
		count += N - i + 1;
	}
	cout << "Done. " << endl << flush;

	return distanceMatrix;
}

vector<double> computePairwiseManhattanDistance(
		const vector<vector<double>> &M) {
	unsigned int N = M.size();
	cout << endl << "Pairwise Manhattan Distance to be computed for " << N
			<< " vectors..." << endl;
	vector<double> distanceMatrix(N * N);

	unsigned int count = 0;

	cout << "Computing distances for all pairs of vectors... " << endl;
	for (unsigned int i = 0; i < N; ++i) {
		printAdvancement(count, (N * (N - 1) / 2));
		distanceMatrix[i + N * i] = 0.0;
		for (unsigned int j = i + 1; j < N; ++j) {
			double distance = manhattanDistance(M[i], M[j]);
			distanceMatrix[i + N * j] = distance;
			distanceMatrix[j + N * i] = distance;
		}
		count += N - i + 1;
	}
	cout << "Done. " << endl << flush;

	return distanceMatrix;
}

vector<double> computePairwiseCosineSimilarity(
		const vector<vector<double>> &M) {
	unsigned int N = M.size();
	cout << endl << "Pairwise Cosine Similarity to be computed for " << N
			<< " vectors..." << endl;
	vector<double> distanceMatrix(N * N);
	vector<double> euclideanNorms(N);
	cout << "Computing all euclidean norms... " << flush;
	transform(M.cbegin(), M.cend(), euclideanNorms.begin(), euclideanNorm);
	cout << "Done." << endl;
	unsigned int count = 0;

	cout << "Computing cosine similarity for all pairs of vectors... " << endl;
	for (unsigned int i = 0; i < N; ++i) {
		printAdvancement(count, (N * (N - 1) / 2));
		distanceMatrix[i + N * i] = 1.0;
		for (unsigned int j = i + 1; j < N; ++j) {
			double dist = cosineSimilarity(M[i], M[j], euclideanNorms[i],
					euclideanNorms[j]);
			distanceMatrix[i + N * j] = dist;
			distanceMatrix[j + N * i] = dist;
		}
		count += N - i + 1;
	}
	cout << "Done. " << endl << flush;
	return distanceMatrix;
}
