#include "functions.h"

point::point(double _x, double _y): x(_x), y(_y) {}

double load() {
	const char* env_p = std::getenv("NCOUNT");
	int n = std::atoi(env_p);
	double h = 0;
	for (int i = 1; i < n; i++) {

		h += std::pow(sin(i), 2) + std::pow(cos(i), 2);
		h -= std::log2((i / i) * 2);
	}
	return h;
}

double f1(double x, double y) {
	double res;
	//double h;
	//h = load();
	res = std::pow(x, 2) + std::pow(y - 1, 2);
	//res += h;
    return res;
}

double f2(double x, double y) {
	double res;
	//double h;
    //h = load();
	res = std::pow(x, 2) - x * y + std::pow(y, 2) - 2 * x + y;
	//res += h;
	return res;
}

double f3(double x, double y) {
	double res;
	//double h;
	//h = load();
	res = std::pow(x, 2) + x * y + std::pow(y, 2) - 4;
	//res += h;
	return res;
}

