#include <vector>
#include "lodePNG/lodepng.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

void makeHeatMap(const vector<double> &matrix, const char* filenameC){
	vector<unsigned char>image;

	double min = *min_element(matrix.cbegin(), matrix.cend());
	double max = 1.0;
	double range = max-min;

	unsigned int n = sqrt(matrix.size());

	for(double d : matrix){
		image.push_back(255*((d-min)/range));
		image.push_back(255*((d-min)/range));
		image.push_back(255*((d-min)/range));
		image.push_back(255);
	}

	string fileName("export/");
	fileName += filenameC;

	 unsigned error = lodepng::encode(fileName.c_str(), image, n, n);

	 //if there's an error, display it
	 if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}
