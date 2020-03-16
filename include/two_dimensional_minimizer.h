#ifndef __TWO_DIMENSIONAL_AGMND_H__
#define __TWO_DIMENSIONAL_AGMND_H__

#include "one_dimensional_agmnd.h"
#include "one_dimensional_agp.h"

class Two_Dimensional_Minimizer : public One_Dimensional_AGP {
private:
	One_Dimensional_Minimizer* odm;
	double lower_y;
	double upper_y;
	double(*function) (double _x, double _y);
	void perform_first_iteration(result* tmp_res), insert_to_map(double _x, double _y, double _z, double _R);
 public:
	Two_Dimensional_Minimizer(One_Dimensional_Minimizer* _odm, double _a, double _b, double _lower_y, double _upper_y,
			  double(*f)(double x, double y), double _eps = 0.001,
			  double _r_par = 2.0);
	void set_experiment(One_Dimensional_Minimizer* _odm, const double _a, const double _b, double _lower_y,
						double _upper_y, double(*f)(double x, double y),
						const double _eps = 0.001,
						const double _r_par = 2.0);
	result solve();
};


#endif