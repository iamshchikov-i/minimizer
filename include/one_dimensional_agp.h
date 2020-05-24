#ifndef __ONE_DIMENSIONAL_AGP_H__
#define __ONE_DIMENSIONAL_AGP_H__

#include "one_dimensional_minimizer.h"

class One_Dimensional_AGP : public One_Dimensional_Minimizer {
protected:
	std::priority_queue<interval, std::vector<interval>, CompareR_min>* pq;
	virtual bool isEnd();
	virtual double get_M(), get_m(), get_R(), get_new_point(interval i);
	virtual void compute_R(double new_point, double new_m),
		insert_to_map(double _y, double _z, double _R, double _num_estimation = 0.0),
		compare_M(double new_point), perform_first_iteration(), delete_containers();
public:
	 void set_experiment(double _a, double _b, double _curr_x,
		 double(*f)(double x, double y),
		 double _eps = 0.001, double _r_par = 2.0);
	One_Dimensional_AGP(double _a, double _b, double _curr_x,
		double(*f)(double x, double y),
		double _eps = 0.001, double _r_par = 2.0);
	virtual result solve();
};

#endif  
