// Copyright Nicolas Legroux 2015

#include "command_line_processor.hpp"

int main(int argc, char *argv[]) {

//	HeinzAnalyzer heinzAnalyzer("0.5-sample0", "biogrid-edges.txt");
//	heinzAnalyzer.analyze();

	CommandLineProcessor clp(argc, argv);
	clp.runProgram();
}
