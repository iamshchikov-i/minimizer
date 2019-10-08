#include "minimizer_v2.h"

Minimizer_v2::Minimizer_v2(double _a, double _b, double(*f)(double x), double _eps, int _N_max, double _r_par): a(_a), b(_b), function(f),
																										eps(_eps), N_max(_N_max), r_p(_r_par) {
	res.k = 0;
	values = new std::map<double, double>;
}

bool Minimizer_v2::stop1() {
	return res.k >= N_max;
}

bool Minimizer_v2::stop2() {
	return min_interval_length <= eps;
}

void Minimizer_v2::reset() {
	left_point = values->begin();
	right_point = left_point; right_point++;
}

void Minimizer_v2::go_Next_Interval() {
	left_point++; right_point++;
}

double Minimizer_v2::get_M() {
	return abs(((*right_point).second - (*left_point).second) / ((*right_point).first - (*left_point).first));
}

double Minimizer_v2::get_M_Max() {
	double M_max, tmp;
	reset();
	double interval_length = (*right_point).first - (*left_point).first;
	min_interval_length = interval_length;
	M_max = get_M();
	go_Next_Interval();

	for (; right_point != values->end(); go_Next_Interval()) {
		tmp = get_M();
		if (tmp >= M_max)
			M_max = tmp;

		interval_length = (*right_point).first - (*left_point).first;
		if (interval_length < min_interval_length)
			min_interval_length = interval_length;
	}
	return M_max;
}

double Minimizer_v2::get_R(double m) {
	double tmp = m * ((*right_point).first - (*left_point).first);
	return tmp + (pow((*right_point).second - (*left_point).second, 2) / tmp) - 2 * ((*right_point).second + (*left_point).second);
}

double Minimizer_v2::get_m() {
	double tmp = get_M_Max();
	if (tmp > 0)
		return r_p * tmp;
	else if (tmp == 0)
		return 1;
	else
		throw - 1;
}

std::pair<double, double> Minimizer_v2::find_R_Max(double m) {
	std::pair<double, double> res;
	double R_max, tmp;
	reset();

	res = std::pair<double, double>((*left_point).first, (*right_point).first);
	R_max = get_R(m);
	go_Next_Interval();

	for (; right_point != values->end(); go_Next_Interval()) {
		tmp = get_R(m);
		if (tmp >= R_max) {
			R_max = tmp;
			res.first = (*left_point).first;
			res.second = (*right_point).first;
		}
	}
	return res;
}

double Minimizer_v2::get_new_point(std::pair<double, double> p, double m) {
	return 0.5*(p.first + p.second) - ((*function)(p.second) - (*function)(p.first)) / (2 * m);
}

void Minimizer_v2::set_experiment(double _a, double _b, double(*f)(double x),
	double _eps, int _N_max, double _r_par) {
	a = _a;
	b = _b;
	function = f;
	eps = _eps;
	N_max = _N_max;
	r_p = _r_par;
	res.k = 0;
	if(values == nullptr)
		values = new std::map<double, double>;
	else throw -1;
}

result Minimizer_v2::solve() {
	std::pair<double, double> new_point;
	double _min;

	//1я итерация
	min_interval_length = b - a;
	values->insert(std::pair<double, double>(a, (*function)(a)));
	values->insert(std::pair<double, double>(b, (*function)(b)));
	if ((*function)(a) <= (*function)(b)) {
		res.x = a;
		_min = (*function)(a);
	}
	else {
		res.x = b;
		_min = (*function)(b);
	}		
	res.k = 2;
	
	while (!stop1()) {
		m = get_m();
		if (stop2()) break;
		new_point.first = get_new_point(find_R_Max(m), m); 
		new_point.second = (*function)(new_point.first); 
		values->insert(new_point); 
		res.k++; 
		if (new_point.second <= _min) {
			res.x = new_point.first;
			res.y = new_point.second;
		}
	}
	delete values;
	values = nullptr;
	return res;
}