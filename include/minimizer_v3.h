#ifndef __MINIMIZER_v_3_H__
#define __MINIMIZER_v_3_H__

// parallel version using MPI

#include "mpi.h"
#include <cmath>
#include <map>
#include <queue>
#include <algorithm>
#include <iostream>
#include "minimizer_v2.h"

class Minimizer_v3 : public Minimizer_v2 {
private:
	int procrank;
	int procnum;
protected:
	void do_first_iteration();
	void do_first_parallel_iteration(std::vector<double>& recvbuf, std::vector<double>& data, std::pair<double, double>& new_point);
public:
	int get_procnum();
	int get_rank();
	result get_result();
	Minimizer_v3(double _a, double _b, double(*f)(double x), double _eps = 0.001, int _N_max = 500, double _r_par = 2.0);
	void set_experiment(const double _a, const double _b, double(*f)(double x), const double _eps = 0.001, const int _N_max = 500, const double _r_par = 2.0);
	void solve();
};

#endif