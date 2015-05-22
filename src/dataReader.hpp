#ifndef DATAREADER_H_INCLUDED
#define DATAREADER_H_INCLUDED

#include <vector>
#include <string>
#include <unordered_map>
#include "typedefs.hpp"

void readData(const std::vector<std::string> &cancers, Data &data, int maxControl = 100, int maxTumor = 1000);
void readData(const std::string &filename, Data &data, int maxControl = 100, int maxTumor = 1000);

void exportDataToFile(const Data &data, const std::string &filename);
void importDataFromFile(Data &data, const std::string &filename);

void exportToMatrix(const Data &data, const std::string &matrixFilename,
		const std::string &patientListFilename);

#endif // DATAREADER_H_INCLUDED
