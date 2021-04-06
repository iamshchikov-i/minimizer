#include <vector>
#include <string.h>
#include <iostream>
#include <sstream>

void parseArguments(int argc, char **argv,
	int& dims,
	double& epsPar, double& rPar, int& Nmax, double& epsErr,
	bool& useMPI,
	bool& useThreads, int& threadsNum);