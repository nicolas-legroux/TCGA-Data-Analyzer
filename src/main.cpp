// Copyright Nicolas Legroux 2015

#include "command_line_processor.hpp"
#include "config.hpp"
#include "utilities.hpp"
#include <fstream>
#include <vector>
#include <iostream>
#include <map>
#include "heinz-analyzer/heinzModuleAnalyzer.hpp"
#include "heinz-analyzer/heinzOutputAnalyzer.hpp"
#include "heinz-analyzer/typedefs.hpp"

int main(int argc, char *argv[]) {
	/*
	 CommandLineProcessor clp(argc, argv);
	 clp.runProgram();
	 */

	HeinzOutputAnalyzer outputAnalyzer("negative-weights.txt", "samples.txt");
	outputAnalyzer.analyze();
}
