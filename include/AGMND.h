#ifndef __AGMND_H__
#define __AGMND_H__

// AGMND sequential version

#include <cmath>
#include <map>
#include <queue>
#include <algorithm>
#include <iostream>

struct result {
	double x;
	double y;
};

struct characteristics {
	double y, R, num_estimation, supposed_x1, supposed_x2, supposed_x3,
		   auxiliary_function_x2, auxiliary_function_x3;
	characteristics();
	characteristics(double _y, double _R, double _num_estimation);
	characteristics(double _R);
};

struct interval {
	std::pair<double, characteristics> first_point;
	std::pair<double, characteristics> second_point;
	interval();
	interval(std::pair<double, characteristics> _f_point,
		     std::pair<double, characteristics> _s_point);
};

struct CompareR {
	bool operator()(interval const& i1, interval const& i2) {
		return i1.second_point.second.R > i2.second_point.second.R;
	}
};

class Minimizer {
protected:
	double a, b, m, M_Max, eps, r_p, min_interval_length;
	double(*function) (double _x);
	std::map<double, characteristics>* points;
	std::map<double, characteristics>::iterator left_point;
	std::map<double, characteristics>::iterator right_point;
	std::priority_queue<interval, std::vector<interval>, CompareR>* pq;
	result res;

	bool isEnd();
	void go_Next_Interval(), go_new_left_interval(double new_point), reset();
	double get_num_estimation(), get_M(), get_A(),
		   get_B(double supposed_x2), get_d(), get_R(), get_m(),
		   auxiliary_function(double x),
		   get_new_point(interval i);
	void compute_R(double new_point, double new_m),
		 insert_to_map(double _x, double _y, double _R, double _num_estimation),
		 compare_interval_len(double new_point), compare_M(double new_point),
		 delete_containers(), perform_first_iteration(), compute_supposed_x();
public:
	Minimizer(double _a, double _b, double (*f)(double x), double _eps = 0.001,
		      double _r_par = 2.0);
	void set_experiment(const double _a, const double _b, double(*f)(double x),
		                const double _eps = 0.001, const double _r_par = 2.0);
	result solve();
	result get_result();
};

#endif