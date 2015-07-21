#ifndef SRC_TCGADATACLUSTERER_HPP_
#define SRC_TCGADATACLUSTERER_HPP_

#include <memory>
#include <map>
#include <vector>
#include <ClusterXX/clustering/algorithms.hpp>

#include "../tcga-analyzer/TCGAData.hpp"

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
	double getAdjustedRandIndex() {
		return clusterer->computeAdjustedRandIndex(realClusters);
	}
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
	TCGADataKMeansClusterer(TCGAData *_ptrToData, unsigned int _K,
			unsigned int _maxIteration, bool verbose);
};

class TCGADataHierarchicalClusterer: public TCGADataClusterer {
public:
	TCGADataHierarchicalClusterer(TCGAData *_ptrToData,
			const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K,
			ClusterXX::HierarchicalParameters::LinkageMethod linkageMethod,
			bool verbose);
	TCGADataHierarchicalClusterer(TCGAData *_ptrToData,
			const Eigen::MatrixXd &distanceMatrix,
			const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K,
			ClusterXX::HierarchicalParameters::LinkageMethod linkageMethod,
			bool verbose);
};

class TCGADataSpectralClusterer: public TCGADataClusterer {
public:
	TCGADataSpectralClusterer(TCGAData *_ptrToData,
			const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K,
			std::pair<
					ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName,
					double> transformationParameters, bool _verbose);
	TCGADataSpectralClusterer(TCGAData *_ptrToData,
			const Eigen::MatrixXd &_distanceMatrix,
			const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K,
			std::pair<
					ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName,
					double> transformationParameters, bool _verbose);
};

#endif /* SRC_TCGADATACLUSTERER_HPP_ */
