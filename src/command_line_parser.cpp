#include "command_line_parser.h"

void parseArguments(int argc, char **argv,
	int& dims,
	double& epsPar, double& rPar, int& Nmax, double& epsErr,
	bool& useMPI,
	bool& useThreads, int& threadsNum) {

	std::vector<std::string> gotArgs;
	for (int i = 0; i < argc; ++i)
		gotArgs.push_back(argv[i]);
	
	std::string delimiter = "=";
	for (auto x : gotArgs) {
		const char *output = NULL;
		
		output = strstr(x.data(), "dims");
		if (output) {
			size_t pos = 0;
			std::string token;
			while ((pos = x.find(delimiter)) != std::string::npos) {
				token = x.substr(0, pos);
				x.erase(0, pos + delimiter.length());
			}
			dims = std::stoi(x);
			continue;
		}

		output = strstr(x.data(), "epsPar");
		if (output) {
			size_t pos = 0;
			std::string token;
			while ((pos = x.find(delimiter)) != std::string::npos) {
				token = x.substr(0, pos);
				x.erase(0, pos + delimiter.length());
			}
			epsPar = std::stod(x);
			continue;
		}

		output = strstr(x.data(), "rPar");
		if (output) {
			size_t pos = 0;
			std::string token;
			while ((pos = x.find(delimiter)) != std::string::npos) {
				token = x.substr(0, pos);
				x.erase(0, pos + delimiter.length());
			}
			rPar = std::stod(x);
			continue;
		}

		output = strstr(x.data(), "Nmax");
		if (output) {
			size_t pos = 0;
			std::string token;
			while ((pos = x.find(delimiter)) != std::string::npos) {
				token = x.substr(0, pos);
				x.erase(0, pos + delimiter.length());
			}
			Nmax = std::stoi(x);
			continue;
		}

		output = strstr(x.data(), "epsErr");
		if (output) {
			size_t pos = 0;
			std::string token;
			while ((pos = x.find(delimiter)) != std::string::npos) {
				token = x.substr(0, pos);
				x.erase(0, pos + delimiter.length());
			}
			epsErr = std::stod(x);
			continue;
		}

		output = strstr(x.data(), "useMPI");
		if (output) {
			size_t pos = 0;
			std::string token;
			while ((pos = x.find(delimiter)) != std::string::npos) {
				token = x.substr(0, pos);
				x.erase(0, pos + delimiter.length());
			}
			std::istringstream(x) >> useMPI;
			continue;
		}

		output = strstr(x.data(), "useThreads");
		if (output) {
			size_t pos = 0;
			std::string token;
			while ((pos = x.find(delimiter)) != std::string::npos) {
				token = x.substr(0, pos);
				x.erase(0, pos + delimiter.length());
			}
			std::istringstream(x) >> useThreads;
			continue;
		}

		output = strstr(x.data(), "threadsNum");
		if (output) {
			size_t pos = 0;
			std::string token;
			while ((pos = x.find(delimiter)) != std::string::npos) {
				token = x.substr(0, pos);
				x.erase(0, pos + delimiter.length());
			}
			threadsNum = std::stoi(x);
			continue;
		}
	}
		
}