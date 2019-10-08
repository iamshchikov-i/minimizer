#include "functions.h"

double f1(double x) {
	return pow(x - 5.46, 2) + 9;
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