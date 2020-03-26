#ifndef __HANSEN_FUNCTIONS_H__
#define __HANSEN_FUNCTIONS_H__

#include <cmath>

struct point {
	double x;
	double y;
	point(double _x, double _y);
};

double (*pfn[])(double x) = {hfunc1, hfunc2, hfunc3, hfunc4, hfunc5,
hfunc6, hfunc7, hfunc8, hfunc9, hfunc10,
hfunc11, hfunc12, hfunc13, hfunc14, hfunc15,
hfunc16, hfunc17, hfunc18, hfunc19, hfunc20};

double intervals[][2] = { {-1.5, 11}, {2.7, 7.5}, {-10.0, 10.0}, {1.9, 3.9}, {0.0, 1.2},
{-10.0, 10.0}, {2.7, 7.5}, {-10.0, 10.0}, {3.1, 20.4}, {0.0, 10.0},
{-1.57, 6.28}, {0.0, 6.28}, {0.001, 0.99}, {0.0, 4.0}, {-5.0, 5.0},
{-3.0, 3.0}, {-4.0, 4.0}, {0.0, 6.0}, {0.0, 6.5}, {-10.0, 10.0} };

double hfunc1(double x);
double hfunc2(double x);
double hfunc3(double x);
double hfunc4(double x);
double hfunc5(double x);
double hfunc6(double x);
double hfunc7(double x);
double hfunc8(double x);
double hfunc9(double x);
double hfunc10(double x);
double hfunc11(double x);
double hfunc12(double x);
double hfunc13(double x);
double hfunc14(double x);
double hfunc15(double x);
double hfunc16(double x);
double hfunc17(double x);
double hfunc18(double x);
double hfunc19(double x);
double hfunc20(double x);
double hpfunc1(double x);
double hpfunc2(double x);
double hpfunc3(double x);
double hpfunc4(double x);
double hpfunc5(double x);
double hpfunc6(double x);
double hpfunc7(double x);
double hpfunc8(double x);
double hpfunc9(double x);
double hpfunc10(double x);
double hpfunc11(double x);
double hpfunc12(double x);
double hpfunc13(double x);
double hpfunc14(double x);
double hpfunc15(double x);
double hpfunc16(double x);
double hpfunc17(double x);
double hpfunc18(double x);
double hpfunc19(double x);
double hpfunc20(double x);

funcs = {hfunc1, hfunc2, hfunc3, hfunc4, hfunc5,
hfunc6, hfunc7, hfunc8, hfunc9, hfunc10,
hfunc11, hfunc12, hfunc13, hfunc14, hfunc15,
hfunc16, hfunc17, hfunc18, hfunc19, hfunc20}
ders = {hpfunc1, hpfunc2, hpfunc3, hpfunc4, hpfunc5,
hpfunc6, hpfunc7, hpfunc8, hpfunc9, hpfunc10,
hpfunc11, hpfunc12, hpfunc13, hpfunc14, hpfunc15,
hpfunc16, hpfunc17, hpfunc18, hpfunc19, hpfunc20}

#endif  
