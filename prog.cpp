#include <iostream>
using std::cout;
using std::endl;
#include <deque>
using std::deque;
#include <string>
using std::string;

#include "OptionParser.h"

Options_t parseCommandline(int argc, char **argv, deque<string> &arguments) {
	OptionParser parser("%prog [options]");
	parser.add_string_option("inputFilename", "-i", "--input-filename",
		"The input file to evaluate.", "/home/panarchos/default.input",
		"INPUT_FILENAME");
	parser.add_string_option("outputFilename", "-o", "--output-filename",
		"The output filename.  Existing files will be overwritten.  The output"
		" format is documented in the README in the installation directory.",
		"/home/panarchos/kamelhaardecke", "OUTPUT_FILENAME");
	parser.add_int_option("multiplier", "-m", "", "The multiplier to use"
		" for the calculation.", 3, "MULTIPLIER_VALUE");
	parser.add_float_option("initialValue", "", "--initial-value",
		"The initial value to use for the calculation.", 1.1);
	parser.add_count_option("verbosity", "-v", "--verbose", "The verbosity level."
		" Repeating this option increases verbosity.", 0);
	parser.add_vector_option("stringList", "-s", "--string-item",
		"A string to process.  This option may be specified more than once.",
		"STRING");
	parser.add_bool_option("fancy", "-f", "--fancy", "Switch on fancy formatting.");
	//cout << parser.get_by_option("-i").info() << endl;
	//cout << parser.get_by_name("help").info() << endl;
	//cout << parser.get_by_name("counter").info() << endl;
	Options_t options = parser.parse_args(argc, argv, arguments);
	//cout << parser.get_by_name("counter").info() << endl;
	parser.print_info() << endl;
	//parser.print_help();
	//cout << "arguments: ";
	//for (auto el : arguments) {
		//cout << el << ", ";
	//}
	//cout << endl;
	return options;
}

int main(int argc, char **argv) {
	deque<string> arguments;
	Options_t options = parseCommandline(argc, argv, arguments);
	//for (auto opt : options) {
		//cout << opt.first << endl;
	//}
	for (const string s : options["stringList"].get_vector()) {
		cout << '"' << s << "\", ";
	}
	cout << endl;
	return 0;
}

