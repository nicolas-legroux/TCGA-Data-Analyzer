/*
 * patientUnsupervisedNormalization.hpp
 *
 *  Created on: Apr 30, 2015
 *      Author: nicolas
 */

#ifndef SRC_UNSUPERVISEDNORMALIZATION_HPP_
#define SRC_UNSUPERVISEDNORMALIZATION_HPP_

#include "typedefs.hpp"

void normalizeKMeans(RNASeqData &controlData, RNASeqData &tumorData, int K,
		int Nmax);
void normalizeQuantile(RNASeqData &controlData, RNASeqData &tumorData,
		double cutPercentage);
void printMaxExpressedGenes(const RNASeqData &controlNormalized,
		const RNASeqData &tumorNormalized, const GeneList &geneList);

#endif /* SRC_UNSUPERVISEDNORMALIZATION_HPP_ */
