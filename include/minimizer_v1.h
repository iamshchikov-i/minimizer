#ifndef __MINIMIZER_v_1__H__
#define __MINIMIZER_v_1__H__

#include <cmath>
#include <map>
#include <algorithm>
#include <iostream>

class Minimizer_v1 {
private:
	int k; 
	double min;  
	const double a;  
	const double b;
	double m; 
	const double r_p; 
	const double eps; 
	const int N_max;  
	std::map<double, double>* values; 
	std::map<double, double>::iterator l;
	std::map<double, double>::iterator r;
	double(*function) (double _x);
protected:
	bool stop1();
	bool stop2();
	void go_Next_Interval();
	void reset();
	double get_M();
	double get_M_Max();
	double get_R(double m);
public:
	Minimizer_v1(double _a, double _b, double (*f)(double x), double _eps = 0.001, int _N_max = 500, double _r_par = 2.0);
	bool isEnd();
	double get_m();
	std::pair<double, double> find_R_Max(double m);
	double get_new_point(std::pair<double, double> p, double m);
	double find_point();
};

#endif