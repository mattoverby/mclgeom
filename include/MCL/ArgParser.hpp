// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_ARGPARSER_HPP
#define MCL_ARGPARSER_HPP 1

#include <sstream>
#include <string>
#include <unordered_map>
#include <fstream>

//
// mcl::ArgParser does simple space-delimited argument parsing
//
// To run a sample program:
//
//	./test_parser --testDouble 0.3 --someflag --testString helloworld --testInt 3 --testFloat 0.33
//
// To load your arguments in the source:
//
//	mcl::ArgParser aParser(argc, argv);
//	double testDouble = aParser.get<double>("--testDouble");
//	float testFloat = aParser.get<float>("--testFloat");
//	int testInt = aParser.get<int>("--testInt");
//	std::string testString = aParser.get<std::string>("--testString");
//	bool someflag_set = aParser.exists("--someflag");
//
// Or if you have defaults that you want to overwrite only if the argument is given:
//
//	double overwriteMe = 3.0;
//	aParser.get<double>("--overwriteMe", &overwriteMe);
//
// And each value will be:
//
//	testDouble: 0.3
//	testFloat: 0.33
//	testInt: 3
//	testString: helloworld
//	someflag_set: true
//	overwriteMe: 3.0
//

namespace mcl
{

class ArgParser
{
public:
	ArgParser(const int &argc, char** argv)
	{
		if (argc < 1) { return; }
		std::stringstream ss_full;
		ss_full << argv[0];
		for (int i=1; i<argc-1; ++i)
		{
			args[argv[i]] = argv[i+1];
			ss_full << ' ' << argv[i];
		}
		if (argc>1)
		{
			args[argv[argc-1]] = "0";
			ss_full << ' ' << argv[argc-1];
		}
		full = ss_full.str();
	}

	// Return whether or not an argument exists
	inline bool exists(const std::string label) const { return (args.count(label) > 0); }

	// Return the value of the argument (or zero, if not exists)
	template<typename T> const T get(const std::string label) const
	{
		if (exists(label))
		{
			std::stringstream ss; ss << args.at(label);
			T value; ss >> value;
			return value;
		}
		return T();
	} // end getter

	// Return true on exists and overwrites value, false otherwise.
	template<typename T> bool get(const std::string label, T *result) const
	{
		if(exists(label))
		{
			std::stringstream ss; ss << args.at(label);
			ss >> *result; return true;
		}
		return false;
	} // end getter to reference

	void save_to_file(const std::string &fn) const
	{
		std::ofstream ofout(fn.c_str());
		if (!ofout.good()) { return; }
		ofout << full;
		ofout.close();
	}

protected:
	std::unordered_map<std::string, std::string> args;
	std::string full;
};


}; // end namespace mcl

#endif
