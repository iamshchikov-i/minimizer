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
		(*left_point).second.num_estimation) / dx;
	t2 = 2 * abs((-dy + (*left_point).second.num_estimation * dx)) / pow(dx, 2);
	t3 = 2 * abs((dy - (*right_point).second.num_estimation * dx)) / pow(dx, 2);

	return std::max({ t1, t2, t3 });
}

double One_Dimensional_AGMND::get_m() {
	return (M_Max > 1.0e-8) ? r_p * M_Max : 1.0;
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
	m = new_m;
	recalc_characteristics();
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
	//go_new_left_interval(new_point);
	M_Max = 0.0;
	for (reset(); right_point != points->end(); go_Next_Interval()) {
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
			(*right_point).second.num_estimation = ((*right_point).second.z - (*left_point).second.z) /
			((*right_point).first - (*left_point).first);
		M_Max = get_M();
		m = -1;
		pq->push(interval({ (*left_point).first, (*left_point).second },
			{ (*right_point).first, (*right_point).second }));
	}
}

void One_Dimensional_AGMND::compute_num_estimation() {
	reset();
	std::map<double, characteristics>::iterator point_after_right = right_point;
	point_after_right++;

	double H, h1, h2, d;

	h1 = (*right_point).first - (*left_point).first;
	h2 = (*point_after_right).first - (*right_point).first;
	d = h2 / h1;
	H = h1 + h2;

	(*left_point).second.num_estimation = (-(2 + d) * (*left_point).second.z +
		(pow((1 + d), 2) * (*right_point).second.z / d) - (*point_after_right).second.z / d) / H;

	for (; point_after_right != --points->end(); go_Next_Interval(), point_after_right++) {
		h1 = (*right_point).first - (*left_point).first;
		h2 = (*point_after_right).first - (*right_point).first;
		d = h2 / h1;
		H = h1 + h2;

		(*right_point).second.num_estimation = (-d * (*left_point).second.z +
			((pow(d, 2) - 1) * (*right_point).second.z / d) + (*point_after_right).second.z / d) / H;
	}

	h1 = (*right_point).first - (*left_point).first;
	h2 = (*point_after_right).first - (*right_point).first;
	d = h2 / h1;
	H = h1 + h2;

	(*right_point).second.num_estimation = (-d * (*left_point).second.z +
		((pow(d, 2) - 1) * (*right_point).second.z / d) + (*point_after_right).second.z / d) / H;

	h1 = (*right_point).first - (*left_point).first;
	h2 = (*point_after_right).first - (*right_point).first;
	d = h2 / h1;
	H = h1 + h2;

	(*point_after_right).second.num_estimation = (d * (*left_point).second.z -
		(pow((1 + d), 2) * (*right_point).second.z / d) + (2 * d + 1) * (*point_after_right).second.z / d) / H;
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
	return ((*right_point).first - (*left_point).first) / 2 +
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
		printf("3 var\n");

		return (*right_point).second.z - (*right_point).second.num_estimation *
			((*right_point).first - x) -
			0.5 * m * pow(x - (*right_point).first, 2);
	}
	else throw - 1;
}

void One_Dimensional_AGMND::compute_supposed_x() {
	double tmp1, tmp2, tmp3, phi2, phi3;
	recalc = false;

	while (true) {
		// A
		tmp1 = ((*left_point).second.z - (*left_point).second.num_estimation *
			(*left_point).first) - ((*right_point).second.z -
			(*right_point).second.num_estimation * (*right_point).first) +
			m * (pow((*right_point).first, 2) - pow((*left_point).first, 2)) / 2;

		// B
		tmp2 = m * ((*right_point).first - (*left_point).first) +
			((*right_point).second.num_estimation -
			(*left_point).second.num_estimation);

		if (tmp2 < 0) {
			printf("tmp2 < 0\n");
			recalc = true;
			m = (1 - ((*right_point).second.num_estimation -
				(*left_point).second.num_estimation)) /
				((*right_point).first - (*left_point).first);
			continue;
		}

		tmp3 = m * pow(get_d(), 2);
		(*right_point).second.supposed_x2 = (tmp1 - tmp3) / tmp2;
		(*right_point).second.supposed_x3 = (tmp1 + tmp3) / tmp2;

		if (((*right_point).second.supposed_x2 < (*left_point).first) ||
			((*right_point).second.supposed_x3 > (*right_point).first)) {

			if ((*right_point).second.supposed_x2 < (*left_point).first)
				(*right_point).second.supposed_x2 = (*left_point).first;

			if ((*right_point).second.supposed_x2 + get_d() >
				(*right_point).first)
				(*right_point).second.supposed_x2 = (*right_point).first - get_d();

			(*right_point).second.supposed_x3 = (*right_point).second.supposed_x2 +
				get_d();

			phi2 = (*right_point).second.z - (*right_point).second.num_estimation *
				((*right_point).first - (*right_point).second.supposed_x3) -
				0.5 * m * pow((*right_point).first - (*right_point).second.supposed_x3, 2);

			phi3 = get_A() * ((*right_point).second.supposed_x3 - (*right_point).second.supposed_x2) +
				0.5 * m * pow((*right_point).second.supposed_x3 - (*right_point).second.supposed_x2, 2) +
				get_B((*right_point).second.supposed_x2);

			if (abs(phi3 - phi2) > abs(phi3) * 0.01) {
				m += 2.0 * fabs(phi3 - phi2) / pow(get_d(), 2);
				recalc = true;
				continue;
			}
		}

		break;
	}
}

void One_Dimensional_AGMND::check_supposed_x() {
	double old_M_Max, a, b, c, D, tmp_x1, tmp_x2, mult;

	old_M_Max = M_Max;
	for (reset(); right_point != points->end(); go_Next_Interval()) {
		a = 0.25 * pow((*right_point).first - (*left_point).first, 2);
		b = -(*right_point).second.z + (*left_point).second.z + 0.5 *
			((*right_point).second.num_estimation + (*left_point).second.num_estimation) *
			((*right_point).first - (*left_point).first);
		c = -pow((*right_point).second.num_estimation - (*left_point).second.num_estimation, 2) / 4;
		D = b * b - 4 * a*c;
		if (D > 0) {
			tmp_x1 = (-b - sqrt(D)) / (2 * a);
			tmp_x2 = (-b + sqrt(D)) / (2 * a);
			M_Max = std::max(M_Max, static_cast<double>(std::max(abs(tmp_x1), abs(tmp_x2))));
		}
	}
	mult = M_Max / old_M_Max;
	if (mult <= 1.0) {
		mult = 1.1;
		M_Max *= mult;
	}
	m = get_m();
}

void One_Dimensional_AGMND::check_new_intervals(double new_point) {
	if (points->size() == 2) {
		reset();
		compute_supposed_x();
		check_supposed_x();
	}
	else {
		go_new_left_interval(new_point);
		for (int i = 0;i < 2;++i, go_Next_Interval()) {
			compute_supposed_x();
			check_supposed_x();
			if (recalc) break;
		}
	}
}

void One_Dimensional_AGMND::update_m() {
	M_Max *= m / new_m;
	m = get_m();
}

void One_Dimensional_AGMND::recalc_characteristics() {
	do {
		while (!pq->empty())
			pq->pop();

		for (reset(); right_point != points->end(); go_Next_Interval()) {
			new_m = m;
			compute_supposed_x();

			if (recalc) {
				update_m();
				reset();
				break;
			}

			if (((*right_point).second.supposed_x2 <= (*left_point).first ||
				(*right_point).second.supposed_x3 >= (*right_point).first)) {
				check_supposed_x();
				reset();
				break;
			}

			(*right_point).second.supposed_x1 = (-(*left_point).second.num_estimation +
				m * ((*right_point).second.supposed_x2 - (*left_point).first) +
				m * (*right_point).second.supposed_x2) / m; // supposed x1
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
	//eps *= b - a;
	std::pair<double, double> new_point;
	double new_m;
	interval search_interval;

	perform_first_iteration(); // perform at the boundary points a and b 

	while (!isEnd()) {
		new_m = get_m();
		compute_R(new_point.first, new_m);
		search_interval = pq->top(); pq->pop();
		new_point.first = get_new_point(search_interval);
		new_point.second = (*function)(curr_x, new_point.first);
		insert_to_map(new_point.first, new_point.second, 0, 0);
		compute_num_estimation();
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	res.k = points->size();
	delete_containers();

	return res;
}
