#ifndef __MINIMIZER_v_2_H__
#define __MINIMIZER_v_2_H__

// improved sequential version

#include <cmath>
#include <map>
#include <queue>
#include <algorithm>
#include <iostream>

struct result {
	int k;
	double x;
	double y;
};

struct characteristics {
	double y;
	double R;
	characteristics();
	characteristics(double _y, double _R);
	characteristics(double _R);
};

struct interval {
	std::pair<double, double> first_point;
	std::pair<double, double> second_point;
	characteristics _ch;
	interval();
	interval(std::pair<double, double> _f_point, std::pair<double, double> _s_point, double _R);
};

struct CompareR {
	bool operator()(interval const& i1, interval const& i2) {
		return i1._ch.R < i2._ch.R;
	}
};

class Minimizer_v2 {
protected:
	std::map<double, characteristics>::iterator left_point;
	std::map<double, characteristics>::iterator right_point;
	double r_p;
	int N_max;
	double eps;
	result res;
	double a;
	double b;
	double(*function) (double _x);
	double min_interval_length;
	double m;
	double M_Max;
	std::map<double, characteristics>* values;
	std::priority_queue<interval, std::vector<interval>, CompareR>* pq;
	bool stop1();
	bool stop2();
	bool isEnd();
	void go_Next_Interval();
	void go_new_left_interval(double new_point);
	void reset();
	double get_M();
	double get_R();
	void calculate_R(double new_point, double new_m);
	double get_m();
	void insert_to_map(double _x, double _y, double _R);
	void compare_interval_len(double new_point);
	void compare_M(double new_point);
	double get_new_point(interval i);
	void delete_containers();
public:
	Minimizer_v2(double _a, double _b, double (*f)(double x), double _eps = 0.001, int _N_max = 500, double _r_par = 2.0);
	void set_experiment(const double _a, const double _b, double(*f)(double x), const double _eps = 0.001, const int _N_max = 500, const double _r_par = 2.0);
	result solve();
	result get_result();
};

#endif