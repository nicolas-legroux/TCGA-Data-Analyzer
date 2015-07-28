/*
 * TCGAPatientData.hpp
 *
 *  Created on: Jul 27, 2015
 *      Author: nicolas
 */

#ifndef SRC_TCGA_ANALYZER_TCGAPATIENTDATA_HPP_
#define SRC_TCGA_ANALYZER_TCGAPATIENTDATA_HPP_

#include <map>
#include <set>

class TCGAPatientData {
public:
	TCGAPatientData(const std::string &_identifier, const std::string &_cancerName, bool _isTumor);
	std::string getPatientName();
	std::string getCancerName();
	bool isTumor();
	void setClinicalData(const std::string &key, const std::string &value);
	bool existsClinicalData(const std::string &key) const;
	std::string getClinicalData(const std::string &key) const;
	std::string toString() const;
	std::string toClassString(std::set<std::string> keys = std::set<std::string>()) const;
private:
	std::string patientName;
	std::string cancerName;
	bool tumor;
	std::map<std::string, std::string> clinicalData;
};

#endif /* SRC_TCGA_ANALYZER_TCGAPATIENTDATA_HPP_ */
