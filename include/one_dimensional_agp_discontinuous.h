#ifndef __ONE_DIMENSIONAL_AGP_DISCONTINUOUS_H__
#define __ONE_DIMENSIONAL_AGP_DISCONTINUOUS_H__

#include "one_dimensional_minimizer.h"

struct CompareR_min {
	bool operator()(interval const& i1, interval const& i2) {
		return i1.second_point.second.R < i2.second_point.second.R;
	}
};

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
	virtual void compute_R(),
		insert_to_map(double _x, double _z, double _R),
		compare_M(double new_point), perform_first_iteration(),
		delete_containers(), compute_u(), set_u_and_delta_null(double new_point),
		compute_p(), set_delta(double new_point);
public:
		void set_experiment(double _a, double _b, double(*f)(double x),
		int _k_par = 5, double _q = 0.1, double _Q = 100,
		double _eps = 0.001, double _r_par = 1.1);
	One_Dimensional_AGP_D(double _a, double _b, double(*f)(double x),
		int _k_par = 5, double _q = 0.1, double _Q = 100,
		double _eps = 0.001, double _r_par = 1.1);
	virtual result solve();
};

#endif