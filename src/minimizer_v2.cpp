#include "minimizer_v2.h"

characteristics::characteristics(double _y, double _R): y(_y), R(_R) {}

characteristics::characteristics(double _R): R(_R) {}

interval::interval(double _f_point, double _s_point, double _R) : first_point(_f_point), second_point(_s_point), _ch(_R) {}

Minimizer_v2::Minimizer_v2(double _a, double _b, double(*f)(double x), double _eps, int _N_max, double _r_par): a(_a), b(_b), function(f),
																										eps(_eps), N_max(_N_max), r_p(_r_par) {
	res.k = 0;
	values = new std::map<double, characteristics>;
	pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
}

bool Minimizer_v2::stop1() {
	return res.k >= N_max;
}

bool Minimizer_v2::stop2() {
	return min_interval_length <= eps;
}

bool Minimizer_v2::isEnd() {
	return stop1() || stop2();
}

void Minimizer_v2::go_Next_Interval() {
	left_point++; right_point++;
}

void Minimizer_v2::go_new_left_interval(double new_point) {
	right_point = values->find(new_point);
	left_point = right_point; left_point--;
}

void Minimizer_v2::reset() {
	left_point = values->begin();
	right_point = left_point; right_point++;
}

double Minimizer_v2::get_M() {
	return abs(((*right_point).second.y - (*left_point).second.y) / ((*right_point).first - (*left_point).first));
}

double Minimizer_v2::get_M_Max() {
	double M_max, tmp;
	reset();
	M_max = get_M();
	go_Next_Interval();

	for (; right_point != values->end(); go_Next_Interval()) {
		tmp = get_M();
		if (tmp >= M_max)
			M_max = tmp;
	}
	return M_max;
}

double Minimizer_v2::get_R() {
	double tmp = m * ((*right_point).first - (*left_point).first);
	return tmp + (pow((*right_point).second.y - (*left_point).second.y, 2) / tmp) - 2 * ((*right_point).second.y + (*left_point).second.y);
}

void Minimizer_v2::calculate_R(double new_point, double new_m) {
	if (new_m != m) {
		if (!pq->empty()) {
			delete pq;
			pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
		}
		m = new_m;
		reset();
		for (; right_point != values->end(); go_Next_Interval()) {
			(*left_point).second.R = get_R();
			pq->push(interval((*left_point).first, (*right_point).first, (*left_point).second.R));
		}
	}
	else {
		go_new_left_interval(new_point);
		for (int i = 0;i < 2;++i, go_Next_Interval()) {
			(*left_point).second.R = get_R();
			pq->push(interval((*left_point).first, (*right_point).first, (*left_point).second.R));
		}
	}
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

void Minimizer_v2::insert_to_map(double _x, double _y, double _R) {
	characteristics _ch(_y, _R);
	values->insert({ _x, _ch });
	if (_y < res.y) {
		res.x = _x;
		res.y = _y;
	}
}

void Minimizer_v2::compare_interval_len(double new_point) {
	go_new_left_interval(new_point);
	double interval_length;
	for (int i = 0;i < 2;++i, go_Next_Interval()) {
		interval_length = (*right_point).first - (*left_point).first;
		if (interval_length < min_interval_length)
			min_interval_length = interval_length;
	}
}

double Minimizer_v2::get_new_point(interval i) {
	return 0.5*(i.first_point + i.second_point) - ((*function)(i.second_point) - (*function)(i.first_point)) / (2 * m);
}

void Minimizer_v2::delete_containers() {
	delete values;
	values = nullptr;
	delete pq;
	pq = nullptr;
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
		values = new std::map<double, characteristics>;
	if(pq == nullptr)
		pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
}

result Minimizer_v2::solve() {
	std::pair<double, double> new_point;
	double new_m;

	//1я итерация
	min_interval_length = b - a;
	insert_to_map(a, (*function)(a), 0);
	res.x = a; res.y = (*function)(a);
	insert_to_map(b, (*function)(b), 0);
	m = -1;
	pq->push(interval(a, b, 0));
	res.k = 2;

	while (!isEnd()) {
		new_m = get_m();
		calculate_R(new_point.first, new_m);
		new_point.first = get_new_point(pq->top()); pq->pop();
		new_point.second = (*function)(new_point.first); 
		insert_to_map(new_point.first, new_point.second, 0);
		res.k++; 
		compare_interval_len(new_point.first);
	}
	delete_containers();
	return res;
}