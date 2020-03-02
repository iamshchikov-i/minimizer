#include "functions.h"

point::point(double _x, double _y): x(_x), y(_y) {}

double load() {
	int n = 500;
	double h = 0;
	for (int i = 1; i < n; i++) {

		h += pow(sin(i), 2) + pow(cos(i), 2);
		h -= log2((i / i) * 2);
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

