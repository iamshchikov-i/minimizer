#ifndef __MULTI_DIMENSIONAL_H__
#define __MULTI_DIMENSIONAL_H__

#include "one_dimensional_agmnd.h"
#include "one_dimensional_agp.h"

#include "mpi.h"

enum class Upper_method { AGP, AGPD, AGMND };

class Multi_Dimensional_Minimizer {
private:
	int range;
	std::vector<One_Dimensional_Minimizer*> odm;
	std::vector<double> lower_bound;
	std::vector<double> upper_bound;
	double(*function)(std::vector<double> coords);
	double eps;
	double r_par;
	int Nmax;
 public:
	Multi_Dimensional_Minimizer(int _range, std::vector<double>& _lower_bound,
		 std::vector<double>& _upper_bound,
		 double(*f)(std::vector<double> coords), Upper_method up, double _eps = 0.01, int _Nmax = 1000,
			  double _r_par = 2.0);
	void set_experiment(int _range, std::vector<double>& _lower_bound,
		std::vector<double>& _upper_bound,
		double(*f)(std::vector<double> coords), double _eps = 0.01,
		double _r_par = 2.0);
	result solve();
};


#endif