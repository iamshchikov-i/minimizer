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
	double y;
	double z;
	int k;
	int k_on_x;
	int k_max_on_y;
};

struct characteristics {
	double z, R, num_estimation, supposed_x1, supposed_x2, supposed_x3,
		auxiliary_function_x2, auxiliary_function_x3;
	characteristics();
	characteristics(double _z, double _R, double _num_estimation);
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
	double a, b, m, M_Max, eps, r_p, min_interval_length, curr_x;
	double(*function) (double _x, double _y);
	std::map<double, characteristics>* points;
	std::map<double, characteristics>::iterator left_point;
	std::map<double, characteristics>::iterator right_point;
	result res;
	void go_Next_Interval(), go_new_left_interval(double new_point), reset(),
		compare_interval_len(double new_point);
	
	virtual bool isEnd() = 0;
	virtual double get_M() = 0, get_m() = 0, get_R() = 0, get_new_point(interval i) = 0;
	virtual void compute_R(double new_point, double new_m) = 0,
		insert_to_map(double _y, double _z, double _R, double _num_estimation) = 0,
		compare_M(double new_point) = 0,
		perform_first_iteration() = 0, delete_containers() = 0;
public:
	One_Dimensional_Minimizer(double _a, double _b, double _curr_x,
		double(*f)(double x, double y),
		double _eps = 0.001, double _r_par = 2.0);
	virtual ~One_Dimensional_Minimizer();
	virtual void set_experiment(double _a, double _b,
		double _curr_x, double(*f)(double x, double y), double _r_p) = 0;
	result get_result();
	virtual result solve() = 0;
	
};

#endif
