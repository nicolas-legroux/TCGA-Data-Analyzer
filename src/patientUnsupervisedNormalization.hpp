/*
 * patientUnsupervisedNormalization.hpp
 *
 *  Created on: Apr 30, 2015
 *      Author: nicolas
 */

#ifndef SRC_PATIENTUNSUPERVISEDNORMALIZATION_HPP_
#define SRC_PATIENTUNSUPERVISEDNORMALIZATION_HPP_

typedef std::unordered_map<std::string, std::vector<std::vector<double>>>RNASeqData;

void normalizeKMeans(RNASeqData &controlData, RNASeqData &tumorData, int K, int Nmax);

#endif /* SRC_PATIENTUNSUPERVISEDNORMALIZATION_HPP_ */
