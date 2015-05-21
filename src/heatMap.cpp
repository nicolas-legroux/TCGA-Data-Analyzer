#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <climits>

#include "lodePNG/lodepng.h"
#include "typedefs.hpp"

using namespace std;

vector<unsigned int> buildClassDivision(vector<SampleIdentifier> &sampleIdentifiers){
	vector<unsigned int> classDivision;
	SampleIdentifier prev = sampleIdentifiers[0];
	unsigned int count = 1;

	for(auto it = sampleIdentifiers.begin()+1; it != sampleIdentifiers.end(); ++it){
		SampleIdentifier next = *it;
		if(prev.cancerName != next.cancerName || prev.isTumor != next.isTumor){
			classDivision.push_back(count);
		}
		++count;
		prev = next;
	}

	return classDivision;
}

void makeHeatMap(const vector<double> &matrix, const char* filenameC,
		vector<unsigned int> classDivision, unsigned int divisionLineWidth) {
	vector<unsigned char> image;

	double min = *min_element(matrix.cbegin(), matrix.cend());
	double max = *max_element(matrix.cbegin(), matrix.cend());
	double range = max - min;

	unsigned int n = sqrt(matrix.size());

	unsigned int lineSeparator = UINT_MAX;
	auto iteratorLineSeparator = classDivision.begin();
	if (iteratorLineSeparator != classDivision.end()) {
		lineSeparator = *iteratorLineSeparator;
	}

	for (unsigned int line = 0; line < n; ++line) {

		if (line == lineSeparator) {
			for (unsigned int i = 0; i < divisionLineWidth; ++i) {
				for (unsigned int j = 0; j < n + classDivision.size() * divisionLineWidth; ++j) {
					image.push_back(255);
					image.push_back(105);
					image.push_back(0);
					image.push_back(255);
				}
			}

			++iteratorLineSeparator;
			if (iteratorLineSeparator == classDivision.end()) {
				lineSeparator = UINT_MAX;
			} else {
				lineSeparator = *iteratorLineSeparator;
			}
		}

		unsigned int columnSeparator = UINT_MAX;
		auto iteratorColumnSeparator = classDivision.begin();
		if (iteratorColumnSeparator != classDivision.end()) {
			columnSeparator = *iteratorColumnSeparator;
		}

		for (unsigned int column = 0; column < n; ++column) {

			if (column == columnSeparator) {
				for (unsigned int i = 0; i < divisionLineWidth; ++i) {
					image.push_back(255);
					image.push_back(105);
					image.push_back(0);
					image.push_back(255);
				}

				++iteratorColumnSeparator;
				if (iteratorColumnSeparator == classDivision.end()) {
					columnSeparator = UINT_MAX;
				} else {
					columnSeparator = *iteratorColumnSeparator;
				}
			}

			double d = matrix[n * line + column];
			image.push_back(255 * ((d - min) / range));
			image.push_back(255 * ((d - min) / range));
			image.push_back(255 * ((d - min) / range));
			image.push_back(255);
		}
	}

	string fileName("export/");
	fileName += filenameC;

	unsigned error = lodepng::encode(fileName.c_str(), image,
			n + classDivision.size() * divisionLineWidth,
			n + classDivision.size() * divisionLineWidth);

	//if there's an error, display it
	if (error)
		std::cout << "encoder error " << error << ": "
				<< lodepng_error_text(error) << std::endl;
}
