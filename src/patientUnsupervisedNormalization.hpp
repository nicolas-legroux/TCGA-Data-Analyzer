/*
 * patientUnsupervisedNormalization.hpp
 *
 *  Created on: Apr 30, 2015
 *      Author: nicolas
 */

#ifndef SRC_PATIENTUNSUPERVISEDNORMALIZATION_HPP_
#define SRC_PATIENTUNSUPERVISEDNORMALIZATION_HPP_

#include "typedefs.hpp"

void normalizeKMeans(RNASeqData &controlData, RNASeqData &tumorData, int K, int Nmax);
void printMaxExpressedGenes(const RNASeqData &controlNormalized, const RNASeqData &tumorNormalized, const GeneList &geneList);

#endif /* SRC_PATIENTUNSUPERVISEDNORMALIZATION_HPP_ */
