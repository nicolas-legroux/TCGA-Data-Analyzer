// Copyright Nicolas Legroux 2015

#include "command_line_processor.hpp"

int main(int argc, char *argv[]) {
	CommandLineProcessor clp(argc, argv);
	clp.runProgram();
}
