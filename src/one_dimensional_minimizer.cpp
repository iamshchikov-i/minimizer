#include "one_dimensional_minimizer.h"

characteristics::characteristics() {}

characteristics::characteristics(double _z, double _R, double _num_estimation) :
	z(_z), R(_R), num_estimation(_num_estimation) {}

interval::interval() {}

interval::interval(std::pair<double, characteristics> _f_point,
	std::pair<double, characteristics> _s_point) :
	first_point(_f_point), second_point(_s_point) {}

void One_Dimensional_Minimizer::go_Next_Interval() {
	left_point++; right_point++;
}

void One_Dimensional_Minimizer::go_new_left_interval(double new_point) {
	right_point = points->find(new_point);
	left_point = right_point; left_point--;
}

void One_Dimensional_Minimizer::reset() {
	left_point = points->begin();
	right_point = left_point; right_point++;
}

void One_Dimensional_Minimizer::compare_interval_len(double new_point) {
	go_new_left_interval(new_point);
	double interval_length;
	for (int i = 0;i < 2;++i, go_Next_Interval()) {
		interval_length = (*right_point).first - (*left_point).first;
		if (interval_length < min_interval_length)
			min_interval_length = interval_length;
	}
}

One_Dimensional_Minimizer::One_Dimensional_Minimizer(double _a, double _b,
	double _curr_x, double(*f)(double x, double y), double _eps,
	double _r_par) : a(_a), b(_b), curr_x(_curr_x), function(f),
	eps(_eps), r_p(_r_par) {
	if (a > b)
		throw "b is a right border, must be more than a";
	points = new std::map<double, characteristics>;
}

One_Dimensional_Minimizer::~One_Dimensional_Minimizer() {
	if (points != nullptr)
		delete points;
	points = nullptr;
}

result One_Dimensional_Minimizer::get_result() {
	return res;
}