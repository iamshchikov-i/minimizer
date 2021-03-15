#ifndef __ONE_DIMENSIONAL_MINIMIZER_H__
#define __ONE_DIMENSIONAL_MINIMIZER_H__

#include <cmath>
#include <map>
#include <queue>
#include <algorithm>
#include <iostream>
#include <utility>
#include <set>

#include "mpi.h"

template <class T>
void print(std::vector<T> vec) {
	for (double item : vec)
		std::cout << item << " ";
	std::cout << std::endl;
}

struct result {
	std::vector<double> coords;
	double z;
	std::vector<int> k;
	int k_on_x;
	int k_max_on_y;
};

struct characteristics {
	double z, R, num_estimation, supposed_x1, supposed_x2, supposed_x3,
		auxiliary_function_x2, auxiliary_function_x3, u;
	int g, delta;
	characteristics();
	characteristics(double _z, double _R, double _num_estimation);
};

struct interval {
	std::pair<std::vector<double>, characteristics> first_point;
	std::pair<std::vector<double>, characteristics> second_point;
	interval();
	interval(std::pair<std::vector<double>, characteristics> _f_point,
		std::pair<std::vector<double>, characteristics> _s_point);
};

struct CompareR_min {
	bool operator()(interval const& i1, interval const& i2) {
		return i1.second_point.second.R < i2.second_point.second.R;
	}
};

class One_Dimensional_Minimizer {
protected:
	int procrank, procnum;
	int curr_dim, range, Nmax;
	double m, M_Max, eps, r_p, min_interval_length;
	std::vector<double> curr_x;
	std::vector<One_Dimensional_Minimizer*> odm;
	std::vector<std::pair<double, double>> bounds;
	double(*function) (std::vector<double> x);
	std::map<std::vector<double>, characteristics>* points;
	std::map<std::vector<double>, characteristics>::iterator left_point;
	std::map<std::vector<double>, characteristics>::iterator right_point;
	result res;
	void go_Next_Interval(), go_new_left_interval(std::vector<double> new_point), reset(),
		compare_interval_len(std::vector<double> new_point);
	
	virtual bool isEnd() = 0;
	virtual double get_M() = 0, get_m() = 0, get_R() = 0, get_new_point(interval i) = 0;
	virtual void compute_R(std::vector<double> new_point, double new_m) = 0,
		insert_to_map(std::vector<double> res_coords, double _z, double _R,
			double _num_estimation = 0.0) = 0,
		compare_M(std::vector<double> new_point) = 0,
		perform_first_iteration() = 0, delete_containers() = 0;
public:
	One_Dimensional_Minimizer(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
		std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
		double(*f)(std::vector<double> x),
		double _eps = 0.001, int _Nmax = 1000, double _r_par = 2.0);
	virtual void set_experiment(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
		std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
		double(*f)(std::vector<double> x),
		double _eps = 0.001, int _Nmax = 1000, double _r_par = 2.0) = 0;
	virtual ~One_Dimensional_Minimizer();

	result get_result();
	virtual result solve() = 0;
	double get_r();
};

#endif
