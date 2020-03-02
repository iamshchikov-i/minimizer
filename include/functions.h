#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <cmath>

struct point {
	double x;
	double y;
	point(double _x, double _y);
};


double load();
double f1(double x, double y);
double f2(double x, double y);
double f3(double x, double y);

#endif  
