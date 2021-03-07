#include "one_dimensional_minimizer.h"

characteristics::characteristics() {}

characteristics::characteristics(double _z, double _R, double _num_estimation) :
	z(_z), R(_R), num_estimation(_num_estimation) {}

interval::interval() {}

interval::interval(std::pair<std::vector<double>, characteristics> _f_point,
	std::pair<std::vector<double>, characteristics> _s_point) :
	first_point(_f_point), second_point(_s_point) {}

void One_Dimensional_Minimizer::go_Next_Interval() {
	left_point++; right_point++;
}

void One_Dimensional_Minimizer::go_new_left_interval(std::vector<double> new_point) {
	right_point = points->find(new_point);
	left_point = right_point; left_point--;
}

void One_Dimensional_Minimizer::reset() {
	left_point = points->begin();
	right_point = left_point; right_point++;
}

void One_Dimensional_Minimizer::compare_interval_len(std::vector<double> new_point) {
	go_new_left_interval(new_point);
	double interval_length;
	for (int i = 0;i < 2;++i, go_Next_Interval()) {
		interval_length = (*right_point).first[curr_dim] - (*left_point).first[curr_dim];
		/*if (curr_dim == 0)
			std::cout<< interval_length << std::endl;*/
		if (interval_length < min_interval_length)
			min_interval_length = interval_length;
	}
}

One_Dimensional_Minimizer::One_Dimensional_Minimizer(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
	std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
	double(*f)(std::vector<double> x),
	double _eps, double _r_par) : range(_range), curr_dim(_curr_dim), odm(_odm),
	bounds(_bounds), curr_x(_curr_x), function(f),
	eps(_eps), r_p(_r_par) {
	
	points = new std::map<std::vector<double>, characteristics>;
}

One_Dimensional_Minimizer::~One_Dimensional_Minimizer() {
	if (points != nullptr)
		delete points;
	points = nullptr;
}

result One_Dimensional_Minimizer::get_result() {
	return res;
}

double One_Dimensional_Minimizer::get_r() {
	return r_p;
}