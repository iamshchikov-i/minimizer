#ifndef __MINIMIZER_v_3_H__
#define __MINIMIZER_v_3_H__

#include <mpi.h>
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
	int curr_involved_procnum;
	MPI_Comm involved_comm, waiting_comm;
protected:
	void calc_involved_procnum(int queue_size);
	void set_comm();
public:
	Minimizer_v3(double _a, double _b, double(*f)(double x), double _eps = 0.001, int _N_max = 500, double _r_par = 2.0);
	result solve();
};

#endif