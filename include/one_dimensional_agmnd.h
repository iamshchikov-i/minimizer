#ifndef __ONE_DIMENSIONAL_AGMND_H__
#define __ONE_DIMENSIONAL_AGMND_H__

// AGMND sequential version
#include "one_dimensional_minimizer.h"

struct CompareR_max {
	bool operator()(interval const& i1, interval const& i2) {
		return i1.second_point.second.R > i2.second_point.second.R;
	}
};

class One_Dimensional_AGMND: public One_Dimensional_Minimizer {
protected:
	bool recalc;
	std::priority_queue<interval, std::vector<interval>, CompareR_max>* pq;
	virtual bool isEnd();
	virtual double get_M() , get_m(), get_R(), get_new_point(interval i);
	virtual void compute_R(std::vector<double> new_point, double new_m),
		insert_to_map(std::vector<double> res_coords, double _z, double _R,
			double _num_estimation),
		compare_M(std::vector<double> new_point),
		perform_first_iteration(), delete_containers();
	double new_m, get_A(), get_B(double supposed_x2), get_d(),
		auxiliary_function(double x);
	void update_m(), compute_num_estimation(), compute_supposed_x(),
		check_supposed_x(), recalc_characteristics(),
		check_new_intervals(std::vector<double> new_point);
public:
	One_Dimensional_AGMND(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
		std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
		double(*f)(std::vector<double> x),
		double _eps = 0.001, int _Nmax = 1000, double _r_par = 2.0);
	virtual result solve();
	void set_experiment(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
		std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
		double(*f)(std::vector<double> x),
		double _eps = 0.001, int _Nmax = 1000, double _r_par = 2.0);
};

#endif
