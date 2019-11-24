#include "functions.h"

double f1(double x) {
	int size = 1000000;
	double h = 0;
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));
	std::vector<int> v(size);
	v.assign(size, 2);
	
	for (int i = 0; i < size/20; i++) {
		h += pow(sin(v.at(gen() % size)), 2) + pow(cos(v.at(gen() % size)), 2);
		h -= log2(v.at(gen() % size));
	}

	return pow(x - 5.46, 2) + 9 + h;
}

double f2(double x) {
	return 5 + x + pow(x, 2);
}

double f3(double x) {
	return pow(x, 3) - 6 * pow(x, 2) + 9 * x - 4;
}

double f4(double x) {
	return (16.4*x) / (1 + pow(x, 2));
}

double f5(double x) {
	return sin(5.87*x);
}

double f6(double x) {
	return sqrt(x)*log(0.1*x);
}