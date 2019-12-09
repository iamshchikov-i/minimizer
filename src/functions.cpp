#include "functions.h"

double load() {
	int size = 1000000;
	double h = 0;

	for (int i = 1; i < size / 15; i++) {

		h += pow(sin(i), 2) + pow(cos(i), 2);
		h -= log2((i / i) * 2);
	}
	return h;
}

double f1(double x) {
	//double h = load();
	double res = pow(x - 5.46, 2) + 9;
	//res += h;
	return res;
}

double f2(double x) {
	//double h = load();
	double res = 5 + x + pow(x, 2);
	//res += h;
	return res;
}

double f3(double x) {
	//double h = load();
	double res = pow(x, 3) - 6 * pow(x, 2) + 9 * x - 4;
	//res += h;
	return res;
}

double f4(double x) {
	//double h = load();
	double res = (16.4*x) / (1 + pow(x, 2));
	//res += h;
	return res;
}

double f5(double x) {
	//double h = load();
	double res = sin(5.87*x);
	//res += h;
	return res;
}

double f6(double x) {
	//double h = load();
	double res = sqrt(x)*log(0.1*x);
	//res += h;
	return res;
}