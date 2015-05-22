#include "dataReader_test.hpp"

#include "../dataReader.hpp"
#include "../typedefs.hpp"

void exportData_test(){
	//Read the data
	std::string filenameCancers = "cancer.list";
	Data data;
	readData(filenameCancers, data);
	exportToMatrix(data, "data_matrix.out", "patients.out");
}
