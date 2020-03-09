#ifndef __AGMND_H__
#define __AGMND_H__

// AGMND sequential version

#include <cmath>
#include <map>
#include <queue>
#include <algorithm>
#include <iostream>

const double r = 4.0;

struct _result {
	double x;
	double y;
	double z;
	int k;
};

struct _characteristics {
	double z, R, num_estimation, supposed_x1, supposed_x2, supposed_x3,
		auxiliary_function_x2, auxiliary_function_x3;
	_characteristics();
	_characteristics(double _z, double _R, double _num_estimation);
};

struct _interval {
	std::pair<double, _characteristics> first_point;
	std::pair<double, _characteristics> second_point;
	_interval();
	_interval(std::pair<double, _characteristics> _f_point,
		std::pair<double, _characteristics> _s_point);
};

struct _CompareR {
	bool operator()(_interval const& i1, _interval const& i2) {
		return i1.second_point.second.R > i2.second_point.second.R;
	}
};

class One_Dimensional_AGMND {
protected:
	double a, b, m, M_Max, eps, r_p, min_interval_length, curr_x;
	double(*function) (double _x, double _y);
	std::map<double, _characteristics>* points;
	std::map<double, _characteristics>::iterator left_point;
	std::map<double, _characteristics>::iterator right_point;
	std::priority_queue<_interval, std::vector<_interval>, _CompareR>* pq;
	_result res;
	bool isEnd();
	void go_Next_Interval(), go_new_left_interval(double new_point), reset();
	double get_num_estimation(), get_M(), get_A(),
		get_B(double supposed_x2), get_d(), get_R(), get_m(),
		auxiliary_function(double x),
		get_new_point(_interval i);
	void compute_R(double new_point, double new_m),
		insert_to_map(double _y, double _z, double _R, double _num_estimation),
		compare_interval_len(double new_point), compare_M(double new_point),
		delete_containers(), perform_first_iteration(), compute_supposed_x();
public:
	One_Dimensional_AGMND(double _a, double _b, double _curr_x,
		double(*f)(double x, double y),
		double _eps = 0.001, double _r_par = r);
	void set_experiment(double _a, double _b,
		double _curr_x, double(*f)(double x, double y), double _eps = 0.001, double _r_par = r);
	_result solve();
	_result get_result();
};

class AGMND : public One_Dimensional_AGMND {
private:
	double lower_y;
	double upper_y;
	double(*function) (double _x, double _y);
	void do_first_iteration(One_Dimensional_AGMND* odm, _result* tmp_res);
public:
	void insert_to_map(double _x, double _y, double _z, double _R, double _num_estimation);
	_result get_result();
	AGMND(double _a, double _b, double _lower_y, double _upper_y,
		double(*f)(double x, double y), double _eps = 0.001, double _r_par = r);
	void set_experiment(const double _a, const double _b, double _lower_y,
		double _upper_y, double(*f)(double x, double y),
		const double _eps = 0.001, const double _r_par = r);
	void solve();
};

#endif