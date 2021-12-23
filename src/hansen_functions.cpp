#include "hansen_functions.h"

double load() {
	int size = 1000000;
	double h = 0;

	for (int i = 1; i < size / 1000; i++) {

		h += std::pow(std::sin(i), 2) + std::pow(std::cos(i), 2);
		h -= log2((i / i) * 2);
	}
	return h;
}

double hfunc1(double x) {
	double h = load();
	return std::pow(x, 6) / 6.0 - 52.0 / 25.0 * std::pow(x, 5) + 39.0 / 80.0 * std::pow(x, 4) +
		71.0 / 10.0 * std::pow(x, 3) - 79.0 / 20.0 * std::pow(x, 2) - x + 0.1 + h;
}

double hfunc2(double x) {
	double h = load();
	return std::sin(x) + std::sin(10 * x / 3) + h;
}

double hfunc3(double x) {
	double h = load();
	double res = 0;
	for (int i = 1; i < 6; i++)
		res += i * std::sin((i + 1) * x + i);
	return -res + h;
}

double hfunc4(double x) {
	double h = load();
	return (-16 * x * x + 24 * x - 5) * exp(-x) + h;
}

double hfunc5(double x) {
	double h = load();
	return -(-3 * x + 1.4) * std::sin(18 * x) + h;
}

double hfunc6(double x) {
	double h = load();
	return -(x + std::sin(x)) * exp(-x * x) + h;
}

double hfunc7(double x) {
	double h = load();
	return std::sin(x) + std::sin(10 * x / 3) + log(x) - 0.84 * x + 3 + h;
}

double hfunc8(double x) {
	double h = load();
	double res = 0;
	for (int i = 1; i < 6; i++)
		res += i * std::cos((i + 1) * x + i);
	return -res + h;
}

double hfunc9(double x) {
	double h = load();
	return std::sin(x) + std::sin(2.0 / 3.0 * x) + h;
}

double hfunc10(double x) {
	double h = load();
	return -x * std::sin(x) + h;
}

double hfunc11(double x) {
	double h = load();
	return 2 * std::cos(x) + std::cos(2 * x) + h;
}

double hfunc12(double x) {
	double h = load();
	return std::pow(std::sin(x), 3) + std::pow(std::cos(x), 3) + h;
}

double hfunc13(double x) {
	double h = load();
	double sgn = 0.0;
	if (x * x - 1 < 0)
		sgn = -1.0;
	else
		sgn = 1.0;
	return -std::pow(x * x, 1.0 / 3.0) + sgn * std::pow(sgn * (x * x - 1.0), 1.0 / 3.0) + h;
}

double hfunc14(double x) {
	double h = load();
	return -exp(-x) * std::sin(2 * std::acos(-1.0) * x) + h;
}

double hfunc15(double x) {
	double h = load();
	return (x * x - 5 * x + 6) / (x * x + 1) + h;
}

double hfunc16(double x) {
	double h = load();
	return 2 * (x - 3) * (x - 3) + exp(x * x / 2) + h;
}

double hfunc17(double x) {
	double h = load();
	return std::pow(x, 6) - 15 * std::pow(x, 4) + 27 * x * x + 250 + h;
}

double hfunc18(double x) {
	double h = load();
	if (x <= 3)
		return (x - 2) * (x - 2) + h;
	else
		return 2 * log(x - 2) + 1 + h;
}

double hfunc19(double x) {
	double h = load();
	return -x + std::sin(3 * x) - 1 + h;
}

double hfunc20(double x) {
	double h = load();
	return -(x - std::sin(x)) * exp(-x * x) + h;
}

double hpfunc1(double x) {
	return std::pow(x, 5) - 10.4 * std::pow(x, 4) + 1.95 * std::pow(x, 3) + 21.3 * x * x -
		7.9 * x - 1.0;
}

double hpfunc2(double x) {
	return std::cos(x) + 10.0 * std::cos(10.0 * x / 3.0) / 3.0;
}

double hpfunc3(double x) {
	double res = 0.0;
	for (int i = 1; i < 6; i++)
		res += i * (i + 1) * std::cos((i + 1) * x + i);
	return -res;
}

double hpfunc4(double x) {
	return (16.0 * x * x - 56.0 * x + 29.0) * exp(-x);
}

double hpfunc5(double x) {
	return 3.0 * std::sin(18.0 * x) - 18.0 * (-3.0 * x + 1.4) * std::cos(18.0 * x);
}

double hpfunc6(double x) {
	return (2.0 * x * (x + std::sin(x)) - std::cos(x) - 1) * exp(-x * x);
}

double hpfunc7(double x) {
	return std::cos(x) + 10.0 * std::cos(10.0 * x / 3.0) / 3.0 + 1 / x - 0.84;
}

double hpfunc8(double x) {
	double res = 0.0;
	for (int i = 1; i < 6; i++)
		res += i * (i + 1) * std::sin((i + 1) * x + i);
	return res;
}

double hpfunc9(double x) {
	return std::cos(x) + 2.0 * std::cos(2.0 * x / 3.0) / 3.0;
}

double hpfunc10(double x) {
	return -std::sin(x) - x * std::cos(x);
}

double hpfunc11(double x) {
	return -2.0 * (std::sin(x) + std::sin(2.0 * x));
}

double hpfunc12(double x) {
	return 3.0 * std::cos(x) * std::sin(x) * (std::sin(x) - std::cos(x));
}

double hpfunc13(double x) {
	double st = (1.0 / 3.0);
	if (x == 0.0)
		return 0.0;
	return (2.0 * x / std::pow((x * x - 1) * (x * x - 1), st) - 2.0 * std::pow(x, -st)) / 3.0;
}

double hpfunc14(double x) {
	double pi = std::acos(-1.0);
	return exp(-x) * (std::sin(2.0 * pi * x) - 2.0 * pi * std::cos(2 * pi * x));
}

double hpfunc15(double x) {
	return (2.0 * x * (-x * x + 5.0 * x - 6.0) - (x * x + 1.0) * (5.0 - 2.0 * x))
		/ ((x * x + 1.0) * (x * x + 1.0));
}

double hpfunc16(double x) {
	return 4.0 * (x - 3.0) + x * exp(x * x / 2.0);
}

double hpfunc17(double x) {
	return 6.0 * std::pow(x, 5) - 60.0 * x * x * x + 54.0 * x;
}

double hpfunc18(double x) {
	if (x <= 3)
		return 2.0 * x - 4.0;
	else
		return 2.0 / (x - 2.0);
}

double hpfunc19(double x) {
	return 3.0 * std::cos(3.0 * x) - 1.0;
}

double hpfunc20(double x) {
	return exp(-x * x) * (2.0 * x * (x - std::sin(x)) - 1.0 + std::cos(x));
}