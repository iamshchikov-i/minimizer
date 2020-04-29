#include "one_dimensional_agmnd.h"

One_Dimensional_AGMND::One_Dimensional_AGMND(double _a, double _b, double _curr_x,
	double(*f)(double x, double y),
	double _eps, double _r_par) : One_Dimensional_Minimizer(_a, _b, _curr_x, f, _eps, _r_par) {
	pq = new std::priority_queue<interval, std::vector<interval>, CompareR_max>;
	recalc = false;
}

bool One_Dimensional_AGMND::isEnd() {
	return min_interval_length <= eps;
}

double One_Dimensional_AGMND::get_M() {
	double t1, t2, t3, dx, dy;

	dx = (*right_point).first - (*left_point).first;
	dy = (*right_point).second.z - (*left_point).second.z;
	t1 = abs((*right_point).second.num_estimation -
		(*left_point).second.num_estimation) / abs(dx);
	t2 = -2.0 * (dy - (*left_point).second.num_estimation * dx) / pow(dx, 2);
	t3 = 2.0 * (dy - (*right_point).second.num_estimation * dx) / pow(dx, 2);

	return std::max({ t1, t2, t3 });
}

double One_Dimensional_AGMND::get_m() {
	if (M_Max > 0)
		return r_p * M_Max;
	else if (M_Max == 0)
		return 1;
	else
		throw - 1;
}

double One_Dimensional_AGMND::get_R() {
	if ((*right_point).second.supposed_x1 >= (*right_point).second.supposed_x2 &&
		(*right_point).second.supposed_x1 <= (*right_point).second.supposed_x3)
		return auxiliary_function((*right_point).second.supposed_x1);
	else {
		(*right_point).second.auxiliary_function_x2 =
			auxiliary_function((*right_point).second.supposed_x2);

		(*right_point).second.auxiliary_function_x3 =
			auxiliary_function((*right_point).second.supposed_x3);

		return std::min((*right_point).second.auxiliary_function_x2,
			(*right_point).second.auxiliary_function_x3);
	}
}

double One_Dimensional_AGMND::get_new_point(interval i) {
	if (i.second_point.second.supposed_x1 >= i.second_point.second.supposed_x2 &&
		i.second_point.second.supposed_x1 <= i.second_point.second.supposed_x3)
		return i.second_point.second.supposed_x1;
	else if (i.second_point.second.auxiliary_function_x2 <=
		i.second_point.second.auxiliary_function_x3)
		return i.second_point.second.supposed_x2;
	else if (i.second_point.second.auxiliary_function_x2 >
		i.second_point.second.auxiliary_function_x3)
		return i.second_point.second.supposed_x3;
	else
		throw - 1;
}

void One_Dimensional_AGMND::compute_R(double new_point, double new_m) {
	if (new_m != m) {
		m = new_m;
		recalc_characteristics();
	}
	else {
		go_new_left_interval(new_point);
		for (int i = 0;i < 2;++i, go_Next_Interval()) {
			compute_supposed_x();
			check_supposed_x();
			if (recalc) break;
			(*right_point).second.R = get_R();
			pq->push(interval({ (*left_point).first, (*left_point).second },
				{ (*right_point).first, (*right_point).second }));
		}

		if (recalc) recalc_characteristics();
	}
}

void One_Dimensional_AGMND::insert_to_map(double _y, double _z, double _R,
	double _num_estimation) {
	characteristics _ch(_z, _R, _num_estimation);
	points->insert({ _y, _ch });
	if (_z < res.z) {
		res.x = curr_x;
		res.y = _y;
		res.z = _z;
	}
}

void One_Dimensional_AGMND::compare_M(double new_point) {
	double M;
	go_new_left_interval(new_point);
	if (left_point == points->begin())
		(*left_point).second.num_estimation =
		(*right_point).second.num_estimation;
	for (int i = 0;i < 2;++i, go_Next_Interval()) {
		(*right_point).second.num_estimation = get_num_estimation();
		M = get_M();
		if (M >= M_Max)
			M_Max = M;
	}
}

void One_Dimensional_AGMND::perform_first_iteration() {
	min_interval_length = b - a;
	res.x = curr_x;
	res.y = a;
	res.z = (*function)(curr_x, a);

	if (a != b) {
		insert_to_map(res.y, res.z, 0, 0);
		insert_to_map(b, (*function)(curr_x, b), 0, 0);
		reset();
		(*left_point).second.num_estimation =
			(*right_point).second.num_estimation = get_num_estimation();
		M_Max = get_M();
		m = -1;
		pq->push(interval({ (*left_point).first, (*left_point).second },
			{ (*right_point).first, (*right_point).second }));
	}
}

double One_Dimensional_AGMND::get_num_estimation() {
	return ((*right_point).second.z - (*left_point).second.z) /
		((*right_point).first - (*left_point).first);
}


double One_Dimensional_AGMND::get_A() {
	return (*left_point).second.num_estimation -
		m * ((*right_point).second.supposed_x2 - (*left_point).first);
}

double One_Dimensional_AGMND::get_B(double x) {
	return (*left_point).second.z +
		(*left_point).second.num_estimation * (x - (*left_point).first) -
		0.5 * m * pow(x - (*left_point).first, 2);
}

double One_Dimensional_AGMND::get_d() {
	return ((*right_point).first - (*left_point).first) / 2 -
		((*right_point).second.num_estimation -
		(*left_point).second.num_estimation) / (2 * m);
}

double One_Dimensional_AGMND::auxiliary_function(double x) {
	if (x >= (*right_point).second.supposed_x2 &&
		x <= (*right_point).second.supposed_x3) {

		return get_A() * (x - (*right_point).second.supposed_x2) +
			0.5 * m * pow(x - (*right_point).second.supposed_x2, 2) +
			get_B((*right_point).second.supposed_x2);
	}

	else if (x > (*left_point).first &&
		x < (*right_point).second.supposed_x2) {

		return (*left_point).second.z + (*left_point).second.num_estimation *
			(x - (*left_point).first) - 0.5 * m * pow(x - (*left_point).first, 2);
	}

	else if (x > (*right_point).second.supposed_x2 &&
		x <= (*right_point).first) {

		return (*right_point).second.z - (*right_point).second.num_estimation *
			(x - (*right_point).first) -
			0.5 * m * pow(x - (*right_point).first, 2);
	}
	else throw - 1;
}

void One_Dimensional_AGMND::compute_supposed_x() {
	double tmp1, tmp2, tmp3;
	tmp1 = ((*left_point).second.z - (*left_point).second.num_estimation *
		(*left_point).first) - ((*right_point).second.z -
		(*right_point).second.num_estimation * (*right_point).first) +
		m * (pow((*right_point).first, 2) - pow((*left_point).first, 2)) / 2;

	tmp2 = m * ((*right_point).first - (*left_point).first) +
		((*right_point).second.num_estimation -
		(*left_point).second.num_estimation);

	tmp3 = m * pow(get_d(), 2);

	(*right_point).second.supposed_x2 = (tmp1 - tmp3) / tmp2;
	(*right_point).second.supposed_x3 = (tmp1 + tmp3) / tmp2;

	(*right_point).second.supposed_x1 = (-(*left_point).second.num_estimation +
		m * ((*right_point).second.supposed_x2 - (*left_point).first) +
		m * (*right_point).first) / m;
}

void One_Dimensional_AGMND::check_supposed_x() {
	double _x = (*right_point).second.supposed_x2, x_ = (*right_point).second.supposed_x3,
		x1 = (*left_point).first, x2 = (*right_point).first,
		z1 = (*left_point).second.z, z2 = (*left_point).second.z,
		z1_ = (*left_point).second.num_estimation, z2_ = (*right_point).second.num_estimation,
		dx = ((*right_point).first - (*left_point).first) / 2 -
		((*right_point).second.num_estimation -
		(*left_point).second.num_estimation) / (2 * m),
		phi2 = 0, phi3 = 0, A = 0, B = 0;

	if ((_x < x1) || (x_ > x2)) {
		if (_x < x1) {
			_x = x1;
			(*right_point).second.supposed_x2 = _x;
		}
			

		if (_x + dx > x2) {
			_x = x2 - dx;
			(*right_point).second.supposed_x2 = _x;
		}
		x_ = _x + dx;
		(*right_point).second.supposed_x3 = x_;

		phi2 = z2 - z2_ * (x2 - x_) - 0.5 * m * (x2 - x_) * (x2 - x_);
		B = z1 + z1_ * (_x - x1) - 0.5 * m * (_x - x1) * (_x - x1);
		A = z1_ - m * (_x - x1);
		phi3 = A * (x_ - _x) + 0.5 * m * (x_ - _x) * (x_ - _x) + B;
		
		if (fabs(phi3 - phi2) > fabs(phi3) * 0.01) {
			m += 2.0 * fabs(phi3 - phi2) / (dx * dx);
			recalc = true;
		}		
	}
}

void One_Dimensional_AGMND::check_new_intervals(double new_point) {
	if (points->size() == 2) {
		reset();
		compute_supposed_x();
		check_supposed_x();
	} else {
		go_new_left_interval(new_point);
		for (int i = 0;i < 2;++i, go_Next_Interval()) {
			compute_supposed_x();
			check_supposed_x();
			if (recalc) break;
		}
	}
}

void One_Dimensional_AGMND::recalc_characteristics() {
	do {
		recalc = false;
		while (!pq->empty())
			pq->pop();
		for (reset(); right_point != points->end(); go_Next_Interval()) {
			compute_supposed_x();
			check_supposed_x();
			if (recalc) break;
			(*right_point).second.R = get_R();
			pq->push(interval({ (*left_point).first, (*left_point).second },
				{ (*right_point).first, (*right_point).second }));
		}
	} while (right_point != points->end());
}

void One_Dimensional_AGMND::delete_containers() {
	delete points;
	points = nullptr;
	delete pq;
	pq = nullptr;
}

void One_Dimensional_AGMND::set_experiment(double _a, double _b,
	double _curr_x, double(*f)(double x, double y), double _r_p) {
	res.k = 0;
	a = _a;
	b = _b;
	function = f;
	r_p = _r_p;
	curr_x = _curr_x;
	recalc = false;
	if (points == nullptr)
		points = new std::map<double, characteristics>;
	if (pq == nullptr)
		pq = new std::priority_queue<interval, std::vector<interval>, CompareR_max>;
}

result One_Dimensional_AGMND::solve() {
	std::pair<double, double> new_point;
	double new_m;
	interval search_interval;

	perform_first_iteration(); // perform at the boundary points a and b 

	while (!isEnd()) {
		new_m = get_m();
		compute_R(new_point.first, new_m);
		search_interval = pq->top(); pq->pop();
		new_point.first = get_new_point(search_interval);

		if (new_point.first <= a || new_point.first >= b)
			printf("ahahhaha lol\n");

		new_point.second = (*function)(curr_x, new_point.first);
		insert_to_map(new_point.first, new_point.second, 0, 0);
		
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	res.k = points->size();
	delete_containers();
	
	return res;
}