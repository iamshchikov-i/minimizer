#include "minimizer_v1.h"

Minimizer_v1::Minimizer_v1(double _a, double _b, double(*f)(double x), double _eps, int _N_max, double _r_par): a(_a), b(_b), function(f),
																										eps(_eps), N_max(_N_max), r_p(_r_par) {
	k = 0;
	values = new std::map<double, double>;
}

bool Minimizer_v1::stop1() {
	return k >= N_max;
}

bool Minimizer_v1::stop2() {
	for (reset(); r != values->end(); go_Next_Interval())
		if (abs((*r).first - (*l).first) <= eps)
			return true;
	return false;
}

void Minimizer_v1::reset() {
	l = values->begin();
	r = l; r++;
}

void Minimizer_v1::go_Next_Interval() {
	l++; r++;
}

double Minimizer_v1::get_M() {
	return abs(((*r).second - (*l).second) / ((*r).first - (*l).first));
}

double Minimizer_v1::get_M_Max() {
	double M_max, tmp;
	reset();
	M_max = get_M();
	go_Next_Interval();

	for (; r != values->end(); go_Next_Interval()) {
		tmp = get_M();
		if (tmp >= M_max)
			M_max = tmp;
	}
	return M_max;
}

double Minimizer_v1::get_R(double m) {
	double tmp = m * ((*r).first - (*l).first);
	return tmp + (pow((*r).second - (*l).second, 2) / tmp) - 2 * ((*r).second + (*l).second);
}

bool Minimizer_v1::isEnd() {
	return stop1() || stop2();
}

double Minimizer_v1::get_m() {
	double tmp = get_M_Max();
	if (tmp > 0)
		return r_p * tmp;
	else if (tmp == 0)
		return 1;
	else
		throw - 1;
}

std::pair<double, double> Minimizer_v1::find_R_Max(double m) {
	std::pair<double, double> res;
	double R_max, tmp;
	reset();

	res = std::pair<double, double>((*l).first, (*r).first);
	R_max = get_R(m);
	go_Next_Interval();

	for (; r != values->end(); go_Next_Interval()) {
		tmp = get_R(m);
		if (tmp >= R_max) {
			R_max = tmp;
			res.first = (*l).first;
			res.second = (*r).first;
		}
	}
	return res;
}

double Minimizer_v1::get_new_point(std::pair<double, double> p, double m) {
	return 0.5*(p.first + p.second) - ((*function)(p.second) - (*function)(p.first)) / (2 * m);
}

double Minimizer_v1::find_point() {
	std::pair<double, double> new_point;
	double _min, res;

	//1я итерация
	values->insert(std::pair<double, double>(a, (*function)(a)));
	values->insert(std::pair<double, double>(b, (*function)(b)));
	if ((*function)(a) <= (*function)(b)) {
		res = a;
		_min = (*function)(a);
	}
	else {
		res = b;
		_min = (*function)(b);
	}		
	k = 2;

	while (!isEnd()) {
		m = get_m();
		new_point.first = get_new_point(find_R_Max(m), m); 
		new_point.second = (*function)(new_point.first); 
		values->insert(new_point); 
		k++; 
		if (new_point.second <= _min)
			res = new_point.first;
	}
	return res;
}

int Minimizer_v1::get_k() {
	return k;
}