#ifndef OPTION_H
#define OPTION_H

#include <stdexcept>
#include <string>
using std::string;
#include <sstream>

class OptionParser;

class Option {
	friend class OptionParser;
public:
	Option() = default;
	Option(const string n, const string t, const string h, const string m,
			bool defB, unsigned defC, double defD, int defI, string defS):
			bool_value(defB), count_value(defC), double_value(defD), int_value(defI),
			string_value(defS), name(n), type(t), help(h) {
		metavar = m;
		if (type != "bool" && type != "count" && metavar.size() == 0) {
			// if no metavar was given, generate one from the supplied name
			for (const char c : name) {
				metavar += std::toupper(c);
			}
		}
	}
	string name = "";
	string help = "";
	string type = "";
	string metavar = "";

	string info();

	// access function for bool_value
	bool is_true();
	bool is_false();
	bool get_bool();
	bool bval();
	// access functions for count_value
	unsigned get_count();
	unsigned count();
	unsigned cval();
	// access functions for double_value
	double get_float();
	double get_double();
	double fval();
	double dval();
	// access functions for int_value
	int get_int();
	int ival();
	// access functions for string_value
	string get_string();
	string str();
	string sval();
	// access functions for vector_value
	std::vector<string> get_vector();
	std::vector<string> vec();
	std::vector<string> vval();
	void set_double_value(double d) { set_double(d); }
private:
	// data members
	bool bool_value = false;
	unsigned count_value = 0;
	double double_value = 0.0;
	int int_value = 0;
	string string_value = "";
	std::vector<string> vector_value;
	// members to write to the data members
	void set_true() { bool_value = true; }
	void set_false() { bool_value = false; }
	void inc_count() { ++count_value; }
	void set_float(double d) { double_value = d; }
	void set_double(double d) { double_value = d; }
	void set_int(int i) { int_value = i; }
	void set_string(string s) { string_value = s; }
	void add_vec_item(string s) { vector_value.push_back(s); }
};


string Option::info() {
	std::ostringstream oss;
	oss << "Option: " << name << ", type: " << type << "\n";
	oss << "\t" << help << "\n";
	oss << "\tstring_value: \"" << string_value << "\", ";
	oss << "int_value: " << int_value << ", ";
	oss << "double_value: " << double_value << ", ";
	oss << "count_value: " << count_value << ", ";
	oss << "bool_value: " << (bool_value ? "true" : "false") << ", ";
	oss << "vector_value: {";
	for (auto it = vector_value.cbegin(); it != vector_value.cend(); ++it) {
		oss << '"' << *it << '"';
		if (it != vector_value.cend() - 1) {
			oss << ", ";
		}
	}
	oss << "}\n";
	return oss.str();
}

// Option's data member's access functions
//
// access function for bool_value
inline bool Option::get_bool() {
	if (type == "bool") {
		return bool_value;
	} else {
		throw std::runtime_error(
			string("Trying to access a bool value, but option type is ") + type
			+ string("."));
	}
}

inline bool Option::is_true() { return get_bool() == true; }
inline bool Option::is_false() { return get_bool() == false; }
inline bool Option::bval() { return get_bool(); }

// access functions for count_value
inline unsigned Option::get_count() {
	if (type == "count") {
		return count_value;
	} else {
		throw std::runtime_error(
			string("Trying to access a count value, but option type is ") + type
			+ string("."));
	}
}

inline unsigned Option::count() { return get_count(); }
inline unsigned Option::cval() { return get_count(); }

// access functions for double_value
inline double Option::get_double() {
	if (type == "double") {
		return double_value;
	} else {
		throw std::runtime_error(
			string("Trying to access a double value, but option type is ") + type
			+ string("."));
	}
}

inline double Option::get_float() { return get_double(); }
inline double Option::dval() { return get_double(); }
inline double Option::fval() { return get_double(); }

// access functions for int_value
inline int Option::get_int() {
	if (type == "int") {
		return int_value;
	} else {
		throw std::runtime_error(
			string("Trying to access a int value, but option type is ") + type
			+ string("."));
	}
}

inline int Option::ival() { return get_int(); }

// access functions for string_value
inline string Option::get_string() {
	if (type == "string") {
		return string_value;
	} else {
		throw std::runtime_error(
			string("Trying to access a string value, but option type is ") + type
			+ string("."));
	}
}

inline string Option::str() { return get_string(); }
inline string Option::sval() { return get_string(); }

// access functions for vector_value
inline std::vector<string> Option::get_vector() {
	if (type == "vector") {
		return vector_value;
	} else {
		throw std::runtime_error(
			string("Trying to access a vector value, but option type is ") + type
			+ string("."));
	}
}

inline std::vector<string> Option::vec() { return get_vector(); }
inline std::vector<string> Option::vval() { return get_vector(); }

#endif

