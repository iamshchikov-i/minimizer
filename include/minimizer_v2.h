#ifndef __MINIMIZER_v_2_H__
#define __MINIMIZER_v_2_H__

#include <cmath>
#include <map>
#include <algorithm>
#include <iostream>

struct result
{
	int k;
	double x;
	double y;
};

class Minimizer_v2 {
private:
	result res;
	double a;  
	double b;
	double m; 
	double r_p; 
	double eps; 
	int N_max;
	double min_interval_length;
	std::map<double, double>* values; 
	std::map<double, double>::iterator left_point;
	std::map<double, double>::iterator right_point;
	double(*function) (double _x);
protected:
	bool stop1();
	bool stop2();
	void go_Next_Interval();
	void reset();
	double get_M();
	double get_M_Max();
	double get_R(double m);
	double get_m();
	std::pair<double, double> find_R_Max(double m);
	double get_new_point(std::pair<double, double> p, double m);
public:
	Minimizer_v2(double _a, double _b, double (*f)(double x), double _eps = 0.001, int _N_max = 500, double _r_par = 2.0);
	void set_experiment(const double _a, const double _b, double(*f)(double x), const double _eps = 0.01, const int _N_max = 1000, const double _r_par = 1.5);
	result solve();
};

#endif