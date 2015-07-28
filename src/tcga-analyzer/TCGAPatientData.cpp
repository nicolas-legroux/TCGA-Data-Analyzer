/*
 * TCGAPatientData.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: nicolas
 */

#include "TCGAPatientData.hpp"

TCGAPatientData::TCGAPatientData(const std::string &_id, const std::string &_cancerName, bool _isTumor) :
		identifier(_id), cancerName(_cancerName), tumor(_isTumor) {
}

std::string TCGAPatientData::getId() {
	return identifier;
}

std::string TCGAPatientData::getCancerName(){
	return cancerName;
}

bool TCGAPatientData::isTumor(){
	return tumor;
}

void TCGAPatientData::setClinicalData(const std::string &key,
		const std::string &value) {
	clinicalData[key] = value;
}

bool TCGAPatientData::existsClinicalData(const std::string &key) const {
	return (clinicalData.find(key) != clinicalData.end());
}

std::string TCGAPatientData::getClinicalData(const std::string &key) const {
	auto result = clinicalData.find(key);
	if (result == clinicalData.end()) {
		throw tcga_data_exception(
				"Could not find clinical data '" + key + "' for patient "
						+ identifier + ".");
	}
	else{
		return result->second;
	}
}

std::string TCGAPatientData::toString() const {
	return cancerName + "_" + ((tumor) ? "Tumor" : "Control") + "_"
			+ identifier;
}

std::string TCGAPatientData::toClassString(std::set<std::string> keys) const {
	std::string className = cancerName + "_" + ((tumor) ? "Tumor" : "Control");
	for(const auto &key : keys){
		className += "_" + getClinicalData(key);
	}
	return className;
}

