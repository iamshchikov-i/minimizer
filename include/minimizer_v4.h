#ifndef __MINIMIZER_v_4_H__
#define __MINIMIZER_v_4_H__

// parallel version using threads

#include <thread>
#include <cmath>
#include <map>
#include <queue>
#include <algorithm>
#include <iostream>
#include "minimizer_v2.h"

class Minimizer_v4 : public Minimizer_v2 {
private:
	int sem;
	int th_num;
protected:
	void prepare_data_for_parallel_job(std::vector<interval>& buf, std::pair<double, double>& new_point);
	void do_shared_job(std::vector<interval>& buf, std::vector<double>& data, int th_rank);
	void do_parallel_job(std::thread *pth, std::vector<interval>& buf, std::vector<double>& data,
	                     std::pair<double, double>& new_point, int rank);
	void insert_data_to_map(std::vector<double>& data);
	void do_first_iteration();
	void do_first_parallel_iteration(std::thread *pth, std::vector<interval>& buf, std::vector<double>& data,
		                             std::pair<double, double>& new_point, int rank);
public:
	int get_procnum();
	result get_result();
	Minimizer_v4(double _a, double _b, double(*f)(double x), int _th_num, double _eps = 0.001, int _N_max = 500, double _r_par = 2.0);
	void set_experiment(const double _a, const double _b, double(*f)(double x), int _th_num, const double _eps = 0.001, const int _N_max = 500, const double _r_par = 2.0);
	void solve();
};

#endif