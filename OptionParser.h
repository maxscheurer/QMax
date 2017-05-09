// Author: Adrian Dempwolff

#ifndef OPTIONPARSER_H
#define OPTIONPARSER_H

#include <iostream>
#include <sstream>
#include <cstdlib> // std::exit
#include <stdexcept>
#include <string>
using std::string;
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>
#include <utility> // pair

#include "Option.h"

using Options_t = std::map<string, Option>;

class OptionParser {
public:
	OptionParser() {
		init_help_option();
	}
	OptionParser(string u): usage(u) {
		init_help_option();
	}

	// member functions to access the stored options
	Option &get_by_name(const string name);
	Option &get_by_specifier(const string option_specifier);
	Options_t &get_options() { return options; }

	// member function to parse the command line arguments
	Options_t &parse_args(int argc, char **argv, std::deque<string> &arguments);
	
	// member functions to print out info, especially the help
	std::ostream &print_info();
	std::ostream &print_help();

	// member functions to add options
	void add_bool_option(const string name, const string short_option,
		const string long_option, const string &help, bool default_value);
	void add_int_option(const string name, const string short_option,
		const string long_option, const string &help, int default_value,
		const string metavar);
	void add_float_option(const string name, const string short_option,
		const string long_option, const string &help, double default_value,
		const string metavar);
	void add_count_option(const string name, const string short_option,
		const string long_option, const string &help, unsigned start_value);
	void add_string_option(const string name, const string short_option,
		const string long_option, const string &help, string default_value,
		const string metavar);
	void add_vector_option(const string name, const string short_option,
		const string long_option, const string &help, const string metavar);
private:
	std::map<string, std::vector<string>> opt_map;
	Options_t options;
	string program_name = "";
	const string program_name_subst = "%prog";
	string usage = "";
	std::vector<string> opt_vector;

	void map_option(const string short_option, const string long_option,
		const string name);
	void check_option(const string option_specifier, const string option_name);
	void init_help_option();
	std::deque<string> wrap_lines(const string &input, const unsigned max_width,
		const unsigned start_spaces);
	std::deque<string> wrap_lines(std::deque<string> &input, const unsigned max_width,
		const unsigned start_spaces);
	void init_option(Option &opt, const string option_arg);
};


void OptionParser::init_option(Option &opt, const string option_arg) {
	// is called by parse_args() and sets the option corresponding to the command
	// line input
	if (opt.name == "help") {
		// This is a special case:  call print_help() and exit
		print_help();
		std::exit(0);
	}
	if (opt.type == "bool") {
		opt.set_true();
		return;
	}
	if (opt.type == "count") {
		opt.inc_count();
		return;
	}
	if (opt.type == "double" || opt.type == "int" || opt.type == "string" ||
			opt.type == "vector") {
		if (option_arg.size() == 0) {
			// This is an error, since we need an argument for these option types
			throw std::invalid_argument(string("Need an argument for the \"")
				+ opt.name + "\" option.");
		}
	}
	if (opt.type == "double") {
		opt.set_float(std::stod(option_arg));
		return;
	}
	if (opt.type == "int") {
		opt.set_int(std::stoi(option_arg));
		return;
	}
	if (opt.type == "string") {
		opt.set_string(option_arg);
		return;
	}
	if (opt.type == "vector") {
		opt.add_vec_item(option_arg);
		return;
	}
}

Options_t &OptionParser::parse_args(int argc, char **argv, std::deque<string> &arguments) {
	// parses the command line
	std::deque<string> args;
	// save the program name
	program_name = string(*argv);
	// store the remaining command line arguments in args
	for (auto it = argv+1; it != argv + argc; ++it) {
		args.push_back(string(*it));
	}
	// catch calls to the help option in advance; init_option will call
	// std::exit(0)
	for (string arg : args) {
		if (find(opt_map["help"].cbegin(), opt_map["help"].cend(), arg) != opt_map["help"].cend()) {
			// arg is in opt_map["help"], i. e. arg is "-h" or "--help"
			init_option(get_by_name("help"), "");
		}
	}
	// process args
	while (args.size()) {
		string curArg = args[0];
		string optArg = "";
		args.pop_front();
		if (curArg.substr(0, 2) == string("--")) {
			// this is a long-style option
			string curOptString;
			auto ind = curArg.find('=');
			if (ind == string::npos) {
				// no '=' found -> this seems to be a bool  or count option
				curOptString = curArg;
			} else {
				curOptString = curArg.substr(0, ind);
				optArg = curArg.substr(ind+1);
			}
			Option &curOpt = get_by_specifier(curOptString);
			// init option curOpt (eventually with optArg) here
			init_option(curOpt, optArg);
		} else {
			if (curArg[0] == '-') {
				// this is a short-style option
				Option &curOpt = get_by_specifier(curArg);
				if (curOpt.type != "bool" && curOpt.type != "count") {
					if (args.size() == 0) {
						// we need a value but there is no, this is a faulty command line
						throw std::invalid_argument(string("Need an argument for the \"")
							+ curArg + "\" option.");
					}
					optArg = args[0];
					args.pop_front();
				}
				// init option curOpt (eventually with optArg) here
				init_option(curOpt, optArg);
			} else {
				// this is an argument, so store it and continue with the while body
				arguments.push_back(curArg);
				continue;
			}
		}
	}
	return options;
}

Option &OptionParser::get_by_name(const string option_name) {
	// returns a reference to the Option object denoted by name option_name
	try {
		return options.at(option_name);
	} catch (std::out_of_range err) {
		throw std::out_of_range(string("No such option name defined: ") + option_name);
	}
}

Option &OptionParser::get_by_specifier(const string option_specifier) {
	// returns a reference to the Option object denoted by option_specifier.
	// look for defined option specifiers whose substrings of the length of
	// option_specifier match option_specifier and vice versa.
	std::vector<std::pair<string, string>> matches;
	for (auto it = opt_map.cbegin(); it != opt_map.cend(); ++it) {
		for (auto iit = it->second.cbegin(); iit != it->second.cend(); ++iit) {
			if (iit->substr(0, option_specifier.size()) == option_specifier
					|| *iit == option_specifier.substr(0, iit->size())) {
				// construct a pair<string, string> and push it at the back of the
				// matches vector;  *iit is the matched option specifier, it->first the
				// corresponding option name
				//matches.push_back(std::pair<string, string>(*iit, it->first));
				matches.emplace_back(*iit, it->first);
			}
		}
	}
	if (matches.size() == 0) {
		// no matching option specifier found in opt_map; this exception is caught
		// when called by check_option()
		throw std::out_of_range(string("No such option specifier defined: ") + option_specifier);
	}
	if (matches.size() > 1) {
		// ambiguous result
		throw std::invalid_argument(string("Ambigous option specification: ")
			+ option_specifier);
	}
	return options[matches[0].second];
}

std::ostream &OptionParser::print_info() {
	// prints some raw info on available options.  Mainly for debugging purposes.
	std::cout << string(80, '-') << std::endl;
	for (auto it = opt_map.cbegin(); it != opt_map.cend(); ++it) {
		std::cout << "Option: ";
		for (std::vector<string>::const_iterator iit = it->second.cbegin();
				iit != it->second.cend(); ++iit) {
			std::cout << *iit << ", ";
		}
		std::cout << "name: " << it->first << std::endl;
		std::cout << options[it->first].info() << std::endl;
		std::cout << string(80, '-') << std::endl;
	}
	return std::cout;
}

std::deque<string> OptionParser::wrap_lines(std::deque<string> &input,
		const unsigned max_width, const unsigned start_spaces) {
	// takes a deque<string>, joins the strings and rewraps them to contain at
	// most max_with characters per line; adds start_spaces space characters to
	// the beginning of each line; the wrapped lines are returned as deque<string>
	// using the overloaded function which takes a string as input
	string input_string;
	for (auto it = input.cbegin(); it != input.cend(); ++it) {
		input_string += *it;
		if (it != input.cend() - 1) {
			input_string += " ";
		}
	}
	return wrap_lines(input_string, max_width, start_spaces);
}

std::deque<string> OptionParser::wrap_lines(const string &input,
		const unsigned max_width, const unsigned start_spaces) {
	// takes a string and rewraps it to contain at most max_with characters per
	// line; adds start_spaces space characters to the beginning of each line; the
	// wrapped lines are returned as deque<string>
	std::deque<string> output;
	string word;
	string line(start_spaces, ' ');
	std::istringstream iss(input);
	while (iss >> word) {
		if (line.size() && line.size() + 1 + word.size() > max_width) {
			// lines containing only one word longer than max_width are not wrapped
			output.push_back(line);
			line = string(start_spaces, ' ');
			line += word;
		} else {
			if (line.size() > start_spaces) {
				line += string(1, ' ');
			}
			line += word;
		}
	}
	if (line.size()) {
		output.push_back(line);
	}
	return output;
}

std::ostream &OptionParser::print_help() {
	// prints out usage and command line options help
	unsigned max_line_width = 85;
	unsigned start_spaces = 24;

	// print usage information
	string usage_string = "Usage: ";
	if (program_name.size()) {
		auto program_name_subst_start = usage.find(program_name_subst);
		if (program_name_subst_start == string::npos) {
			usage_string += usage;
		} else {
			usage_string += usage.replace(program_name_subst_start,
				program_name_subst.size(), program_name);
		}
	} else {
		usage_string += usage;
	}
	std::deque<string> usage_output = wrap_lines(usage_string, max_line_width, 0);
	for (string el : usage_output) {
		std::cout << el << std::endl;
	}

	// print options information
	std::cout << "\n" << "Options:" << std::endl;
	// prepare the option help output
	for (auto opt_vec_it = opt_vector.cbegin(); opt_vec_it != opt_vector.cend();
			++opt_vec_it) {
		auto it = opt_map.find(*opt_vec_it);
		const string &help_string = options[it->first].help;
		// cycle through the options in opt_map; prepare an option_string from the
		// option specifiers
		string option_string(2, ' ');
		for (auto iit = it->second.cbegin(); iit != it->second.cend(); ++iit) {
			option_string += *iit;
			// care abount metavar stuff
			if (options[it->first].type != "bool" && options[it->first].type != "count") {
				// we don't have a metavar for bool or count type options
				if (iit->substr(0, 2) == string("--")) {
					// this is a long option specifier
					option_string.append("=");
				} else {
					// this is a short option specifier
					option_string.append(" ");
				}
				option_string.append(options[it->first].metavar);
			}
			if (iit != it->second.cend()-1) {
				option_string += string(", ");
			}
		}
		string first_help_line;
		std::deque<string> wrapped_help_for_first_line;
		if (option_string.size() < start_spaces) {
			// fill option_string up with spaces
			option_string += string(start_spaces - option_string.size(), ' ');
			wrapped_help_for_first_line = wrap_lines(
				help_string, max_line_width - option_string.size(), 0);
			first_help_line = option_string + wrapped_help_for_first_line[0];
		} else {
			// non-whitespace characters till the end, so start with the help text in
			// the following line
			std::cout << option_string << std::endl;
			wrapped_help_for_first_line = wrap_lines(
				help_string, max_line_width, start_spaces);
			first_help_line = wrapped_help_for_first_line[0];
		}
		// print out the first help line, pop this first line
		std::cout << first_help_line << std::endl;
		wrapped_help_for_first_line.pop_front();
		if (wrapped_help_for_first_line.size() == 0) {
			continue;
		}
		std::deque<string> wrapped_help_rest = wrap_lines(
			wrapped_help_for_first_line, max_line_width, start_spaces);
		for (string el : wrapped_help_rest) {
			std::cout << el << std::endl;
		}
	}
	std::cout << std::endl;
	return std::cout;
}

void OptionParser::check_option(const string option_specifier, const string option_name) {
	// ensures that no option with the same ore incompatible option specifiers
	// was previously defined
	try {
		Option opt = get_by_specifier(option_specifier);
		if (opt.name != option_name) {
			// This is a problem:  Option spcifier already used for another option
			string specifier_string = "\"";
			for (auto nit = opt_map[opt.name].cbegin(); nit !=
					opt_map[opt.name].cend(); ++nit) {
				specifier_string += *nit;
				if (nit != opt_map[opt.name].cend() - 1) {
					specifier_string += string(", ");
				}
			}
			specifier_string += string("\"");
			throw std::invalid_argument(string("Option clash: Option \"")
				+ option_name + string("\" with specifier \"")
				+ option_specifier + string("\" and option \"")
				+ opt.name + string("\" with specifiers ") + specifier_string
				+ string(" have incompatible specifiers."));
		}
	} catch (std::out_of_range err) {
		// if we are here, everything is ok since this option specifier is not
		// known till now
	}
}

void OptionParser::map_option(const string short_option,
		const string long_option, const string name) {
	// adds a mapping between option name and specifiers to opt_map.
	// it is ensured that a) no option with the same name and b)
	// by calling check_option() that no option with an incompatible
	// specifier is already defined.
	if (short_option.empty() && long_option.empty()) {
		throw std::invalid_argument(
			string("No option specifier supplied for option: \"")
			+ name + string("\"."));
	}
	// at least one of short_option and long_option is not empty
	std::vector<string> option_specifiers;
	if (!short_option.empty()) {
		if (short_option[0] != '-') {
			throw std::invalid_argument(
				string("Short option specifier should start with a hyphen '-'."));
		}
		check_option(short_option, name);
		option_specifiers.push_back(short_option);
	}
	if (!long_option.empty()) {
		if (long_option.substr(0, 2) != string("--")) {
			throw std::invalid_argument(
				string("Long option specifier should start with two hyphens \"--\"."));
		}
		check_option(long_option, name);
		option_specifiers.push_back(long_option);
	}
	// try to insert the {name, option_specifier} pair vector into opt_map
	auto insert_result = opt_map.insert({name, option_specifiers});
	if (!insert_result.second) {
		// insert didn't succeed, because an option with the same name is
		// already defined in opt_map
		throw std::invalid_argument(string("Option \"")
			+ name + "\" is already defined (is defined twice).");
	} else {
		// inserting succedded, so append this option to opt_vector
		opt_vector.push_back(name);
	}
}

void OptionParser::add_bool_option(const string name, const string short_option = "",
		const string long_option = "", const string &help = "", bool default_value = false) {
	// creates a bool Option object and adds it to options
	if (opt_map.count(name)) {
		throw std::invalid_argument(string("Option clash: identifier \"")
			+ name + string("\" already defined."));
	}
	map_option(short_option, long_option, name);
	Option o(name, "bool", help, "", default_value, 0, 0.0, 0, "");
	options[name] = o;
}

void OptionParser::add_count_option(const string name, const string short_option = "",
		const string long_option = "", const string &help = "", unsigned start_value = 0) {
	// creates a count Option object and adds it to options
	map_option(short_option, long_option, name);
	Option o(name, "count", help, "", false, start_value, 0.0, 0, "");
	options[name] = o;
}

void OptionParser::add_float_option(const string name, const string short_option = "",
		const string long_option = "", const string &help = "", double default_value = 0.0,
		const string metavar = "") {
	// creates a double Option object and adds it to options
	map_option(short_option, long_option, name);
	Option o(name, "double", help, metavar, false, 0, default_value, 0, "");
	options[name] = o;
}

void OptionParser::add_int_option(const string name, const string short_option = "",
		const string long_option = "", const string &help = "", int default_value = 0,
		const string metavar = "") {
	// creates a int Option object and adds it to options
	map_option(short_option, long_option, name);
	Option o(name, "int", help, metavar, false, 0, 0.0, default_value, "");
	options[name] = o;
}

void OptionParser::add_string_option(const string name, const string short_option = "",
		const string long_option = "", const string &help = "", string default_value = "",
		const string metavar = "") {
	// creates a string Option object and adds it to options
	map_option(short_option, long_option, name);
	Option o(name, "string", help, metavar, false, 0, 0.0, 0, default_value);
	options[name] = o;
}

void OptionParser::add_vector_option(const string name, const string short_option,
		const string long_option, const string &help, const string metavar) {
	// creates a vector Option object and adds it to options
	map_option(short_option, long_option, name);
	Option o(name, "vector", help, metavar, false, 0, 0.0, 0, "");
	options[name] = o;
}

void OptionParser::init_help_option() {
	// initializes the default -h / --help option
	map_option("-h", "--help", "help");
	Option o("help", "bool", "Print this help message and exit.", "", false, 0, 0.0, 0, "");
	options["help"] = o;
}

#endif

