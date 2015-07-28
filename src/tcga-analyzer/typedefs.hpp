/*
 * typedefs.hpp
 *
 *  Created on: May 11, 2015
 *      Author: nicolas
 */

#ifndef SRC_TYPEDEFS_HPP_
#define SRC_TYPEDEFS_HPP_

#include <unordered_map>
#include <vector>

class tcga_data_exception: public std::exception {
public:
	 tcga_data_exception(std::string _msg = "TCGA Data Exception") :
			msg(_msg) {
	}
	~ tcga_data_exception() throw () {
	}
	const char* what() const throw () {
		return msg.c_str();
	}
private:
	std::string msg;
};

// cancerName -> geneID -> patientID -> RNASeq value
typedef std::vector<std::vector<double>> RNASeqData;
// geneID -> (HNSC Symbol, Entrez ID)
typedef std::vector<std::pair<std::string, int>> GeneList;
// class -> list of Patient IDs
typedef std::map<std::string, std::vector<int>> ClassMap;

#endif /* SRC_TYPEDEFS_HPP_ */
