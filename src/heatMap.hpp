#ifndef SRC_HEATMAP_HPP_
#define SRC_HEATMAP_HPP_

#include <vector>
#include "typedefs.hpp"

std::vector<unsigned int> buildClassDivision(
		std::vector<SampleIdentifier> &sampleIdentifiers);

void makeHeatMap(const MatrixX &matrix, const char* filename,
		std::vector<unsigned int> classDivision = std::vector<unsigned int>(),
		unsigned int divisionLineWidth = 0);

#endif /* SRC_HEATMAP_HPP_ */
