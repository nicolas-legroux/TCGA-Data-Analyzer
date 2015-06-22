#ifndef SRC_TCGADATACLUSTERER_HPP_
#define SRC_TCGADATACLUSTERER_HPP_

#include <memory>
#include <map>
#include <vector>
#include <ClusterXX/clustering/algorithms.hpp>
#include "TCGAData.hpp"
#include "config.hpp"

enum ClusteringMethod {
	KMEANS_CLUSTERING, SPECTRAL_CLUSTERING, HIERARCHICAL_CLUSTERING
};

class TCGADataClusterer {
public:
	TCGADataClusterer(TCGAData *_ptrToData, unsigned int _K, bool _verbose);
	virtual ~TCGADataClusterer() = 0;
	void computeClustering();
	std::vector<int> getClusters();
	virtual void printClusteringInfo();
protected:
	TCGAData *ptrToData;
	unsigned int K;
	bool verbose;
	std::map<int, std::string> realLabelsMap;
	std::vector<int> realClusters;
	std::shared_ptr<ClusterXX::ClustererParameters> clustererParameters; //To be initialized in children class
	std::shared_ptr<ClusterXX::Clusterer> clusterer; //To be initialized in children class
private:
	void buildRealClusters();
	void buildRealLabelsMap();
};

class TCGADataKMeansClusterer: public TCGADataClusterer {
public:
	TCGADataKMeansClusterer(TCGAData *_ptrToData, unsigned int _K = 0,
			unsigned int _maxIterations = K_MEANS_MAX_ITERATIONS, bool verbose =
					VERBOSE);
};

class TCGADataHierarchicalClusterer: public TCGADataClusterer {
public:
	TCGADataHierarchicalClusterer(TCGAData *_ptrToData,
			const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K =
					0,
			ClusterXX::HierarchicalParameters::LinkageMethod linkageMethod =
					DEFAULT_LINKAGE_METHOD, bool verbose = VERBOSE);
	TCGADataHierarchicalClusterer(TCGAData *_ptrToData,
			const Eigen::MatrixXd &distanceMatrix,
			const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K =
					0,
			ClusterXX::HierarchicalParameters::LinkageMethod linkageMethod =
					DEFAULT_LINKAGE_METHOD, bool verbose = VERBOSE);
};

class TCGADataSpectralClusterer: public TCGADataClusterer {
public:
	TCGADataSpectralClusterer(TCGAData *_ptrToData,
			const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K =
					0,
			std::pair<
					ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName,
					double> transformationParameters =
					DEFAULT_GRAPH_TRANSFORMATION, bool _verbose = VERBOSE);
	TCGADataSpectralClusterer(TCGAData *_ptrToData,
			const Eigen::MatrixXd &_distanceMatrix,
			const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K =
					0,
			std::pair<
					ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName,
					double> transformationParameters =
					DEFAULT_GRAPH_TRANSFORMATION, bool _verbose = VERBOSE);
};

#endif /* SRC_TCGADATACLUSTERER_HPP_ */
