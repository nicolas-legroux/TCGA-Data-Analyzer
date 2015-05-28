#include "hierarchicalClustering_test.hpp"
#include "../hierarchical_clustering.hpp"
#include "../utilities.hpp"
#include <vector>

using namespace std;

void hierarchicalClusteringTest() {
	vector<double> data { 1, 0.8, 0.9, 0.2, 0.3, 0.8, 1, 0.95, 0.1, 0.2, 0.9,
			0.95, 1, 0.3, 0.4, 0.2, 0.1, 0.4, 1, 0.7, 0.3, 0.2, 0.2, 07, 1 };

	Hierarchical_Clustering hierarchicalClustering(data, LinkageMethod::COMPLETE, MatrixType::SIMILARITY);
	vector<int> clusters = hierarchicalClustering.compute(2);
	print_vector(clusters);
}
