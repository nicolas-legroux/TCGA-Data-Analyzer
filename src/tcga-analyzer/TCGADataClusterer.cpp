#include "TCGADataClusterer.hpp"

TCGADataClusterer::TCGADataClusterer(TCGAData *_ptrToData, unsigned int _K,
		bool _verbose) :
		ptrToData(_ptrToData), K(_K), verbose(_verbose) {
	realClusters.resize(ptrToData->getNumberOfSamples());
	buildRealClasses();
	ptrToData->buildDataMatrix({}, verbose);
	if (K == 0) {
		K = ptrToData->getClassMapHandler().size();
	}
}

TCGADataClusterer::~TCGADataClusterer() {

}

void TCGADataClusterer::buildRealClasses(){
	int currentCluster = 0;
	for(const auto &kv : ptrToData->getClassMapHandler()){
		realLabels.push_back( kv.first);
		for(auto i : kv.second){
			realClusters[i] = currentCluster;
		}
		++currentCluster;
	}
}

void TCGADataClusterer::computeClustering() {
	clusterer->compute();
}

std::vector<int> TCGADataClusterer::getClusters() {
	return clusterer->getClusters();
}

void TCGADataClusterer::printClusteringInfo() {
	clusterer->printClusteringMatrix(realLabels, realClusters);
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

TCGADataUnnormalizedSpectralClusterer::TCGADataUnnormalizedSpectralClusterer(TCGAData *_ptrToData,
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
	clusterer = std::make_shared<ClusterXX::UnnormalizedSpectralClustering>(
			ptrToData->getDataMatrixHandler(), clustererParameters);
}

TCGADataUnnormalizedSpectralClusterer::TCGADataUnnormalizedSpectralClusterer(TCGAData *_ptrToData,
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
	clusterer = std::make_shared<ClusterXX::UnnormalizedSpectralClustering>(_distanceMatrix,
			clustererParameters, true);
}

TCGADataNormalizedSpectralClusterer::TCGADataNormalizedSpectralClusterer(TCGAData *_ptrToData,
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
	clusterer = std::make_shared<ClusterXX::NormalizedSpectralClustering>(
			ptrToData->getDataMatrixHandler(), clustererParameters);
}

TCGADataNormalizedSpectralClusterer::TCGADataNormalizedSpectralClusterer(TCGAData *_ptrToData,
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
	clusterer = std::make_shared<ClusterXX::NormalizedSpectralClustering>(_distanceMatrix,
			clustererParameters, true);
}

TCGADataNormalizedSpectralClusterer_RandomWalk::TCGADataNormalizedSpectralClusterer_RandomWalk(TCGAData *_ptrToData,
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
	clusterer = std::make_shared<ClusterXX::NormalizedSpectralClustering_RandomWalk>(
			ptrToData->getDataMatrixHandler(), clustererParameters);
}

TCGADataNormalizedSpectralClusterer_RandomWalk::TCGADataNormalizedSpectralClusterer_RandomWalk(TCGAData *_ptrToData,
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
	clusterer = std::make_shared<ClusterXX::NormalizedSpectralClustering_RandomWalk>(_distanceMatrix,
			clustererParameters, true);
}
