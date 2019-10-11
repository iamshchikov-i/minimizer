#ifndef __MINIMIZER_v_2_H__
#define __MINIMIZER_v_2_H__

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
	characteristics(double _y, double _R);
	characteristics(double _R);
};

struct interval {
	double first_point;
	double second_point;
	characteristics _ch;
	interval(double _f_point, double _s_point, double _R);
};

struct CompareR {
	bool operator()(interval const& i1, interval const& i2) {
		return i1._ch.R < i2._ch.R;
	}
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
	std::map<double, characteristics>* values;
	std::map<double, characteristics>::iterator left_point;
	std::map<double, characteristics>::iterator right_point;
	std::priority_queue<interval, std::vector<interval>, CompareR>* pq;
	double(*function) (double _x);
protected:
	bool stop1();
	bool stop2();
	bool isEnd();
	void go_Next_Interval();
	void go_new_left_interval(double new_point);
	void reset();
	double get_M();
	double get_M_Max();
	double get_R();
	void calculate_R(double new_point, double new_m);
	double get_m();
	void insert_to_map(double _x, double _y, double _R);
	void compare_interval_len(double new_point);
	double get_new_point(interval i);
	void delete_containers();
public:
	Minimizer_v2(double _a, double _b, double (*f)(double x), double _eps = 0.01, int _N_max = 500, double _r_par = 2.0);
	void set_experiment(const double _a, const double _b, double(*f)(double x), const double _eps = 0.01, const int _N_max = 500, const double _r_par = 2.0);
	result solve();
};

#endif