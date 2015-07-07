#include "TCGADataClusterer.hpp"

TCGADataClusterer::TCGADataClusterer(TCGAData *_ptrToData, unsigned int _K,
		bool _verbose) :
		ptrToData(_ptrToData), K(_K), verbose(_verbose) {
	buildRealLabelsMap();
	buildRealClusters();
	ptrToData->transposeData(verbose);
	if (K == 0) {
		K = realLabelsMap.size();
	}
}

TCGADataClusterer::~TCGADataClusterer() {

}

void TCGADataClusterer::buildRealClusters() {

	Sample prev = ptrToData->getSamplesHandler()[0];
	int currentCluster = 0;
	realClusters.push_back(currentCluster);

	for (auto it = ptrToData->getSamplesHandler().begin() + 1;
			it != ptrToData->getSamplesHandler().end(); ++it) {
		Sample next = *it;
		if (prev.cancerName != next.cancerName
				|| prev.isTumor != next.isTumor) {
			++currentCluster;
		}
		realClusters.push_back(currentCluster);
		prev = next;
	}
}

void TCGADataClusterer::buildRealLabelsMap() {
	Sample prev = ptrToData->getSamplesHandler()[0];
	std::string label = prev.cancerName + "-"
			+ ((prev.isTumor) ? "Tumor" : "Control");
	int currentCluster = 0;
	realLabelsMap.insert(std::make_pair(currentCluster, label));

	for (auto it = ptrToData->getSamplesHandler().begin() + 1;
			it != ptrToData->getSamplesHandler().end(); ++it) {
		Sample next = *it;
		if (prev.cancerName != next.cancerName
				|| prev.isTumor != next.isTumor) {
			++currentCluster;
			std::string label = next.cancerName + "-"
					+ ((next.isTumor) ? "Tumor" : "Control");
			realLabelsMap.insert(std::make_pair(currentCluster, label));
		}
		prev = next;
	}
}

void TCGADataClusterer::computeClustering() {
	clusterer->compute();
}

std::vector<int> TCGADataClusterer::getClusters() {
	return clusterer->getClusters();
}

void TCGADataClusterer::printClusteringInfo() {
	clusterer->printClustering(realLabelsMap, realClusters);
	std::cout << std::endl << "Adjusted Rand Index : " << clusterer->computeAdjustedRandIndex(realClusters) << std::endl;
}

TCGADataKMeansClusterer::TCGADataKMeansClusterer(TCGAData *_ptrToData,
		unsigned int _K, unsigned int _maxIterations, bool _verbose) :
		TCGADataClusterer(_ptrToData, _K, _verbose) {
	clustererParameters = std::make_shared<ClusterXX::KMeansParameters>(K,
			_maxIterations, _verbose);
	clusterer = std::make_shared<ClusterXX::KMeans_Clusterer>(
			ptrToData->getDataMatrixHandler(), clustererParameters);
}

TCGADataHierarchicalClusterer::TCGADataHierarchicalClusterer(
		TCGAData *_ptrToData, const std::shared_ptr<ClusterXX::Metric> &_metric,
		unsigned int _K,
		ClusterXX::HierarchicalParameters::LinkageMethod linkageMethod,
		bool _verbose) :
		TCGADataClusterer(_ptrToData, _K, verbose) {
	clustererParameters = std::make_shared<ClusterXX::HierarchicalParameters>(K,
			_metric, linkageMethod, _verbose);
	clusterer = std::make_shared<ClusterXX::Hierarchical_Clusterer>(
			ptrToData->getDataMatrixHandler(), clustererParameters);
}

TCGADataHierarchicalClusterer::TCGADataHierarchicalClusterer(
		TCGAData *_ptrToData, const Eigen::MatrixXd &_distanceMatrix,
		const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K,
		ClusterXX::HierarchicalParameters::LinkageMethod linkageMethod,
		bool _verbose) :
		TCGADataClusterer(_ptrToData, _K, verbose) {
	clustererParameters = std::make_shared<ClusterXX::HierarchicalParameters>(K,
			_metric, linkageMethod, _verbose);
	clusterer = std::make_shared<ClusterXX::Hierarchical_Clusterer>(
			_distanceMatrix, clustererParameters, true);

}

TCGADataSpectralClusterer::TCGADataSpectralClusterer(TCGAData *_ptrToData,
		const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K,
		std::pair<
				ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName,
				double> transformationParameters, bool _verbose) :
		TCGADataClusterer(_ptrToData, _K, _verbose) {
	clustererParameters = std::make_shared<ClusterXX::SpectralParameters>(K,
			_metric,
			ClusterXX::SpectralParameters::GraphTransformationMethod(
					transformationParameters.first,
					transformationParameters.second), _verbose);
	clusterer = std::make_shared<ClusterXX::Spectral_Clusterer>(
			ptrToData->getDataMatrixHandler(), clustererParameters);
}

TCGADataSpectralClusterer::TCGADataSpectralClusterer(TCGAData *_ptrToData,
		const Eigen::MatrixXd &_distanceMatrix,
		const std::shared_ptr<ClusterXX::Metric> &_metric, unsigned int _K,
		std::pair<
				ClusterXX::SpectralParameters::GraphTransformationMethod::GraphTransformationMethodName,
				double> transformationParameters, bool _verbose) :
		TCGADataClusterer(_ptrToData, _K, _verbose) {
	clustererParameters = std::make_shared<ClusterXX::SpectralParameters>(K,
			_metric,
			ClusterXX::SpectralParameters::GraphTransformationMethod(
					transformationParameters.first,
					transformationParameters.second), _verbose);
	clusterer = std::make_shared<ClusterXX::Spectral_Clusterer>(_distanceMatrix,
			clustererParameters, true);
}
