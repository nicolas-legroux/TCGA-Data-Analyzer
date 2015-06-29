/*
 * CommandLineProcessor.hpp
 *
 *  Created on: Jun 29, 2015
 *      Author: nicolas
 */

#ifndef SRC_COMMAND_LINE_PROCESSOR_HPP_
#define SRC_COMMAND_LINE_PROCESSOR_HPP_

#include <set>
#include <stdexcept>
#include <map>

class wrong_usage_exception: public std::exception {
public:
	wrong_usage_exception(std::string _msg = "Wrong Usage Exception") :
			msg(_msg) {
	}
	~wrong_usage_exception() throw () {
	}
	const char* what() const throw () {
		return msg.c_str();
	}
private:
	std::string msg;
};

const std::set<std::string> allowedOptions = { "mode", "cancers", "maxcontrol",
		"maxtumor", "negativeweights", "binaryquantileparam" };

class CommandLineProcessor {
public:
	CommandLineProcessor(int argc, char *argv[]);
	void runProgram();
private:
	void process(const std::string &optionName, const std::string &optionValue);
	std::map<std::string, std::string> environment;
};

#endif /* SRC_COMMAND_LINE_PROCESSOR_HPP_ */
