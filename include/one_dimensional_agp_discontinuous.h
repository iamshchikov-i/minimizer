#ifndef __ONE_DIMENSIONAL_AGP_DISCONTINUOUS_H__
#define __ONE_DIMENSIONAL_AGP_DISCONTINUOUS_H__

#include "one_dimensional_minimizer.h"

class One_Dimensional_AGP_D : public One_Dimensional_Minimizer {
private:
	int k_par, p;
	double q, Q, u_max;
	std::map<double, std::map<double, characteristics>::iterator, 
		std::greater<double>> u_sorted;
	std::map<double, std::map<double, characteristics>::iterator,
		std::greater<double>>::iterator u_left, u_right;
protected:
	std::priority_queue<interval, std::vector<interval>, CompareR_min>* pq;
	virtual bool isEnd();
	virtual double get_M(), get_m(), get_R(), get_new_point(interval i), get_u();
	virtual void compute_R(double new_point, double new_m),
		insert_to_map(double _y, double _z, double _R, double _num_estimation = 0.0),
		compare_M(double new_point), perform_first_iteration(),
		delete_containers(), compute_u(), set_u_and_delta_null(double new_point),
		compute_p(), set_delta(double new_point), set_experiment(double _a, double _b, double _curr_x,
			double(*f)(double x, double y),
			double _eps = 0.001, double _r_par = 2.0);
public:
	void set_experiment_d(double _a, double _b, double _curr_x, 
		double(*f)(double x, double y), int _k_par = 10, double _q = 0.1,
		double _Q = 100, double _eps = 0.001, double _r_p = 2.0);
	One_Dimensional_AGP_D(double _a, double _b, double _curr_x,
		double(*f)(double x, double y), int _k_par = 10, double _q = 0.1,
		double _Q = 100, double _eps = 0.001, double _r_p = 2.0);
	virtual result solve();
};

#endif