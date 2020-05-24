#include "one_dimensional_agp_discontinuous.h"

One_Dimensional_AGP_D::One_Dimensional_AGP_D(double _a, double _b,
	double(*f)(double x), int _k_par, double _q, double _Q,
	double _eps, double _r_par) : One_Dimensional_Minimizer(_a, _b, f, _eps, _r_par),
								  k_par(_k_par), q(_q), Q(_Q) {
	pq = new std::priority_queue<interval, std::vector<interval>, CompareR_min>;
}
bool One_Dimensional_AGP_D::isEnd() {
	return min_interval_length <= eps;
}

double One_Dimensional_AGP_D::get_M() {
	return std::abs(((*right_point).second.z - (*left_point).second.z) /
		((*right_point).first - (*left_point).first));
}

double One_Dimensional_AGP_D::get_m() {
	if (M_Max > 0) {
		return r_p * M_Max;
	}
	else if (M_Max == 0) {
		return 1;
	}
	else {
		throw - 1;
	}
}

double One_Dimensional_AGP_D::get_R() {
	double dx = (*right_point).first - (*left_point).first;
	return (1 + (*right_point).second.g + (*left_point).second.g +
		abs((*right_point).second.delta)) * dx +
		(1 - (*right_point).second.g - (*left_point).second.g -
		abs((*right_point).second.delta)) * (pow((*right_point).second.z - (*left_point).second.z, 2) /
			pow(r_p * u_max, 2) * dx) - 
		2 * ((1 - (*right_point).second.g) * (1 + (*left_point).second.g)  *
			(1 - (*right_point).second.delta) * (*right_point).second.z + 
			(1 + (*right_point).second.g) * (1 - (*left_point).second.g) * 
			(1 + (*right_point).second.delta) * (*left_point).second.z) / 
			(r_p * u_max);
}

double One_Dimensional_AGP_D::get_new_point(interval i) {
	return 0.5*(i.first_point.first + i.second_point.first) -
		(1 - i.second_point.second.g - i.first_point.second.g - 
			abs(i.second_point.second.delta)) *
			(i.second_point.second.z - i.first_point.second.z) / (2 * r_p * u_max);
}

double One_Dimensional_AGP_D::get_u() {
	return (1 - (*right_point).second.g - (*left_point).second.g -
		abs((*right_point).second.delta))*
		(abs((*right_point).second.z - (*left_point).second.z)) /
			((*right_point).first - (*left_point).first);
}

void One_Dimensional_AGP_D::compute_u() {
	reset();
	u_sorted.clear();
	for (; right_point != points->end(); go_Next_Interval()) {
		(*right_point).second.u = get_u();
		u_sorted.insert({ (*right_point).second.u, right_point });
		if ((*right_point).second.u > u_max)
			u_max = (*right_point).second.u;
	}
	if (u_max == 0.0)
		u_max = 1.0;
}

void One_Dimensional_AGP_D::compute_p() {
	p = 0;
	u_left = u_sorted.begin();
	u_right = u_left; u_right++;
	for (int i = 1;u_right != u_sorted.end(); u_right++, i++) {
		if ((u_left->first) / (u_right->first) >= Q)
			if (i >= 1 && i < q * u_sorted.size()) {
				p = i;
				return;
			}
	}
	p = 0;
}

void One_Dimensional_AGP_D::set_u_and_delta_null(double new_point) {
	go_new_left_interval(new_point);
	(*left_point).second.g = (*right_point).second.g = 0;
	(*right_point).second.delta = 0;
}

void One_Dimensional_AGP_D::set_delta(double new_point) {
	go_new_left_interval(new_point);
	go_Next_Interval();
	(*left_point).second.delta = (*right_point).second.delta;
	(*right_point).second.delta = 0;
	(*left_point).second.g = abs((*left_point).second.delta);
}

void One_Dimensional_AGP_D::compute_R() {
	while (!pq->empty())
		pq->pop();
	reset();
	for (; right_point != points->end(); go_Next_Interval()) {
		if (p > 0 && p < points->size() &&
			(*right_point).second.u >= (u_left->first))
			(*right_point).second.delta =
			(*right_point).second.z >= (*left_point).second.z ? 1 : -1;
		else
			(*right_point).second.delta = 0;
			
		(*right_point).second.R = get_R();
		pq->push(interval({ (*left_point).first, (*left_point).second },
			{ (*right_point).first, (*right_point).second }));
	}
}

void One_Dimensional_AGP_D::insert_to_map(double _x, double _z,
	double _R) {
	characteristics _ch(_z, _R);
	points->insert({ _x, _ch });
	if (_z < res.z) {
		res.x = _x;
		res.z = _z;
	}
}

void One_Dimensional_AGP_D::compare_M(double new_point) {
	go_new_left_interval(new_point);
	double M;
	for (int i = 0; i < 2; ++i, go_Next_Interval()) {
		M = get_M();
		if (M >= M_Max)
			M_Max = M;
	}
}

void One_Dimensional_AGP_D::perform_first_iteration() {
	min_interval_length = b - a;
	res.x = a;
	res.z = (*function)(a);
	if (a != b) {
		insert_to_map(res.x, res.z, 0);
		insert_to_map(b, (*function)(b), 0);
		reset();
		(*left_point).second.g = (*right_point).second.g = 0;
		(*right_point).second.delta = 0;
		p = 0;
		u_max = get_u();
		pq->push(interval({ (*left_point).first, (*left_point).second },
			{ (*right_point).first, (*right_point).second }));
	}
}

void One_Dimensional_AGP_D::delete_containers() {
	delete points;
	points = nullptr;
	delete pq;
	pq = nullptr;
}

void One_Dimensional_AGP_D::set_experiment(double _a, double _b, double(*f)(double x),
	int _k_par, double _q, double _Q,
	double _eps, double _r_par) {
	if (a > b)
		throw "b is a right border, must be more than a";
	k_par = _k_par;
	q = _q;
	Q = _Q;
	res.k = 0;
	a = _a;
	b = _b;
	r_p = _r_par;
	function = f;
	if (points == nullptr)
		points = new std::map<double, characteristics>;
	if (pq == nullptr)
		pq = new std::priority_queue<interval, std::vector<interval>,
		CompareR_min>;
}

result One_Dimensional_AGP_D::solve() {
	std::pair<double, double> new_point;

	perform_first_iteration();

	while (min_interval_length > eps && points->size() <= k_par ) {
		compute_u();
		compute_R();
		new_point.first = get_new_point(pq->top()); pq->pop();
		new_point.second = (*function)(new_point.first);
		insert_to_map(new_point.first, new_point.second, 0);
		compare_interval_len(new_point.first);
		set_u_and_delta_null(new_point.first);
	}
	
	while (min_interval_length > eps) {
		compute_u();
		compute_p();
		compute_R();
		new_point.first = get_new_point(pq->top()); pq->pop();
		new_point.second = (*function)(new_point.first);
		insert_to_map(new_point.first, new_point.second, 0);
		compare_interval_len(new_point.first);
		set_delta(new_point.first);
	}
	res.k = points->size();
	delete_containers();

	return res;
}