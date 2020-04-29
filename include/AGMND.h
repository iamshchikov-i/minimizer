#ifndef __AGMND_H__
#define __AGMND_H__

// AGMND sequential version

#include <cmath>
#include <map>
#include <queue>
#include <algorithm>
#include <iostream>

struct _result {
	double x;
	double y;
	int k;
};

struct _characteristics {
	double y, R, num_estimation, supposed_x1, supposed_x2, supposed_x3,
		   auxiliary_function_x2, auxiliary_function_x3;
	_characteristics();
	_characteristics(double _y, double _R, double _num_estimation);
	_characteristics(double _R);
};

struct _interval {
	std::pair<double, _characteristics> first_point;
	std::pair<double, _characteristics> second_point;
	_interval();
	_interval(std::pair<double, _characteristics> _f_point,
		     std::pair<double, _characteristics> _s_point);
};

struct CompareR_max {
	bool operator()(_interval const& i1, _interval const& i2) {
		return i1.second_point.second.R > i2.second_point.second.R;
	}
};

class Minimizer {
protected:
	bool recalc;
	double a, b, m, new_m, M_Max, eps, r_p, min_interval_length;
	double(*function) (double _x);
	std::map<double, _characteristics>* points;
	std::map<double, _characteristics>::iterator left_point;
	std::map<double, _characteristics>::iterator right_point;
	std::priority_queue<_interval, std::vector<_interval>, CompareR_max >* pq;
	_result res;
	bool isEnd();
	void go_Next_interval(), go_new_left_interval(double new_point), reset();
	double get_M(), get_A(),
		   get_B(double supposed_x2), get_d(), get_R(), get_m(),
		   auxiliary_function(double x),
		   get_new_point(_interval i);
	void compute_R(double new_point), compute_num_estimation(),
		insert_to_map(double _x, double _y, double _R, double _num_estimation),
		compare_interval_len(double new_point), compare_M(double new_point),
		delete_containers(), perform_first_iteration(), compute_supposed_x(),
		recalc_characteristics(), check_supposed_x(), update_m();
public:
	Minimizer(double _a, double _b, double (*f)(double x), double _eps = 0.001,
		      double _r_par = 2.0);
	void set_experiment(const double _a, const double _b, double(*f)(double x),
		                const double _eps = 0.001, const double _r_par = 2.0);
	_result solve();
	_result get_result();
};

#endif