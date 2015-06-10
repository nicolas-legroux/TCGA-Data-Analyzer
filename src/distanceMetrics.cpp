#include <vector>
#include <cmath>

#include "distanceMetrics.hpp"
#include "utilities.hpp"
#include "stats.hpp"

using namespace std;

MatrixX computePairwisePearsonCorrelation(const MatrixX &data) {
	unsigned int N = data.cols();
	cout << endl << "Pearson correlation to be computed for " << N
			<< " vectors..." << endl;
	MatrixX correlationMatrix(N, N);
	vector<double> means(N);
	vector<double> standard_deviations(N);

	cout << "Computing all means and standard deviations... " << flush;
	for (unsigned int i = 0; i < N; ++i) {
		means[i] = computeMean(data.col(i));
		standard_deviations[i] = computeStandardDeviation(data.col(i));
	}
	cout << "Done." << endl;

	unsigned int count = 0;

	cout << "Computing correlations for all pairs of vectors... " << endl;
	for (unsigned int i = 0; i < N; ++i) {
		printAdvancement(count, (N * (N + 1) / 2));
		correlationMatrix(i, i) = 1.0;
		for (unsigned int j = i + 1; j < N; ++j) {
			double cor = computePearsonCorrelation(data.col(i), data.col(j),
					means[i], standard_deviations[i], means[j],
					standard_deviations[j]);
			correlationMatrix(i, j) = cor;
			correlationMatrix(j, i) = cor;
		}
		count += N - i;
	}
	cout << "Done. " << endl << flush;

	return correlationMatrix;
}

MatrixX computePairwiseSpearmanCorrelation(const MatrixX &data) {
	vector<vector<double>> data_copy_stl(data.cols());
	for (unsigned int i = 0; i < data.cols(); ++i) {
		for (unsigned int j = 0; j < data.rows(); ++j) {
			double d = data(j, i);
			data_copy_stl[i].push_back(d);
		}
	}
	for_each(data_copy_stl.begin(), data_copy_stl.end(), computeRank);
	MatrixX data_copy(data.rows(), data.cols());
	for (unsigned int i = 0; i < data.cols(); ++i) {
		for (unsigned int j = 0; j < data.rows(); ++j) {
			data_copy(j, i) = data_copy_stl[i][j];
		}
	}
	return computePairwisePearsonCorrelation(data_copy);
}

MatrixX computePairwiseEuclideanDistance(const MatrixX &data) {
	unsigned int N = data.cols();
	cout << endl << "Pairwise Euclidean Distance to be computed for " << N
			<< " vectors..." << endl;
	MatrixX distanceMatrix(N, N);

	unsigned int count = 0;

	cout << "Computing distances for all pairs of vectors... " << endl;
	for (unsigned int i = 0; i < N; ++i) {
		printAdvancement(count, (N * (N + 1) / 2));
		distanceMatrix(i, i) = 0.0;
		for (unsigned int j = i + 1; j < N; ++j) {
			double distance = euclideanDistance(data.col(i), data.col(j));
			distanceMatrix(i, j) = distance;
			distanceMatrix(j, i) = distance;
		}
		count += N - i;
	}
	cout << "Done. " << endl << flush;

	return distanceMatrix;
}

MatrixX computePairwiseManhattanDistance(const MatrixX &data) {
	unsigned int N = data.cols();
	cout << endl << "Pairwise Manhattan Distance to be computed for " << N
			<< " vectors..." << endl;
	MatrixX distanceMatrix(N, N);

	unsigned int count = 0;

	cout << "Computing distances for all pairs of vectors... " << endl;
	for (unsigned int i = 0; i < N; ++i) {
		printAdvancement(count, (N * (N + 1) / 2));
		distanceMatrix(i, i) = 0.0;
		for (unsigned int j = i + 1; j < N; ++j) {
			double distance = manhattanDistance(data.col(i), data.col(j));
			distanceMatrix(i, j) = distance;
			distanceMatrix(j, i) = distance;
		}
		count += N - i;
	}
	cout << "Done. " << endl << flush;

	return distanceMatrix;
}

MatrixX computePairwiseCosineSimilarity(const MatrixX &data) {
	unsigned int N = data.cols();
	cout << endl << "Pairwise Cosine Similarity to be computed for " << N
			<< " vectors..." << endl;
	MatrixX distanceMatrix(N, N);
	vector<double> euclideanNorms(N);
	cout << "Computing all euclidean norms... " << flush;
	for (unsigned int i = 0; i < N; ++i) {
		euclideanNorms[i] = data.col(i).norm();
	}
	cout << "Done." << endl;
	unsigned int count = 0;

	cout << "Computing cosine similarity for all pairs of vectors... " << endl;
	for (unsigned int i = 0; i < N; ++i) {
		printAdvancement(count, (N * (N + 1) / 2));
		distanceMatrix(i, i) = 1.0;
		for (unsigned int j = i + 1; j < N; ++j) {
			double dist = cosineSimilarity(data.col(i), data.col(j),
					euclideanNorms[i], euclideanNorms[j]);
			distanceMatrix(i, j) = dist;
			distanceMatrix(j, i) = dist;
		}
		count += N - i;
	}
	cout << "Done. " << endl << flush;
	return distanceMatrix;
}
