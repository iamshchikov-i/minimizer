#ifndef __ONE_DIMENSIONAL_MINIMIZER_H__
#define __ONE_DIMENSIONAL_MINIMIZER_H__

#include <cmath>
#include <map>
#include <queue>
#include <algorithm>
#include <iostream>
#include <utility>

struct result {
	double x;
	double z;
	int k;
};

struct characteristics {
	double z, R, u;
	int g, delta;
	characteristics();
	characteristics(double _z, double _R);
};

struct interval {
	std::pair<double, characteristics> first_point;
	std::pair<double, characteristics> second_point;
	interval();
	interval(std::pair<double, characteristics> _f_point,
		std::pair<double, characteristics> _s_point);
};

class One_Dimensional_Minimizer {
protected:
	double a, b, m, M_Max, eps, r_p, min_interval_length;
	double(*function) (double _x);
	std::map<double, characteristics>* points;
	std::map<double, characteristics>::iterator left_point;
	std::map<double, characteristics>::iterator right_point;
	result res;
	void go_Next_Interval(), go_new_left_interval(double new_point), reset(),
		compare_interval_len(double new_point);
	
	virtual bool isEnd() = 0;
	virtual double get_M() = 0, get_m() = 0, get_R() = 0, get_new_point(interval i) = 0;
	virtual void compute_R() = 0,
		insert_to_map(double _x, double _z, double _R) = 0,
		compare_M(double new_point) = 0, perform_first_iteration() = 0,
		delete_containers() = 0;
public:
	One_Dimensional_Minimizer(double _a, double _b, double(*f)(double x),
		double _eps = 0.001, double _r_par = 2.0);
	virtual ~One_Dimensional_Minimizer();
	result get_result();
	virtual result solve() = 0;
	double get_r();
};

#endif