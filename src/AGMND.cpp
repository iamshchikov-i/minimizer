#include "AGMND.h"

_characteristics::_characteristics() {}

_characteristics::_characteristics(double _z, double _R, double _num_estimation) :
	z(_z), R(_R), num_estimation(_num_estimation) {}

_interval::_interval() {}

_interval::_interval(std::pair<double, _characteristics> _f_point,
	std::pair<double, _characteristics> _s_point) :
	first_point(_f_point), second_point(_s_point) {}

One_Dimensional_AGMND::One_Dimensional_AGMND(double _a, double _b,
	double _curr_x, double(*f)(double x, double y), double _eps,
	double _r_par) : a(_a), b(_b), curr_x(_curr_x), function(f),
	eps(_eps), r_p(_r_par) {
	if (a > b)
		throw "b is a right border, must be more than a";
	points = new std::map<double, _characteristics>;
	pq = new std::priority_queue<_interval, std::vector<_interval>, _CompareR>;
}

double One_Dimensional_AGMND::get_num_estimation() {
	return ((*right_point).second.z - (*left_point).second.z) /
		((*right_point).first - (*left_point).first);
}

bool One_Dimensional_AGMND::isEnd() {
	return min_interval_length <= eps;
}

void One_Dimensional_AGMND::go_Next_Interval() {
	left_point++; right_point++;
}

void One_Dimensional_AGMND::go_new_left_interval(double new_point) {
	right_point = points->find(new_point);
	left_point = right_point; left_point--;
}

void One_Dimensional_AGMND::reset() {
	left_point = points->begin();
	right_point = left_point; right_point++;
}

double One_Dimensional_AGMND::get_M() {
	double t1, t2, t3, dx, dy;

	dx = (*right_point).first - (*left_point).first;
	dy = (*right_point).second.z - (*left_point).second.z;
	t1 = abs((*right_point).second.num_estimation -
		(*left_point).second.num_estimation) / abs(dx);
	t2 = -2 * (dy - (*left_point).second.num_estimation * dx) / pow(dx, 2);
	t3 = 2 * (dy - (*right_point).second.num_estimation * dx) / pow(dx, 2);

	return std::max({ t1, t2, t3 });
}

double One_Dimensional_AGMND::get_d() {
	return ((*right_point).first - (*left_point).first) / 2 -
		((*right_point).second.num_estimation -
		(*left_point).second.num_estimation) / (2 * m);
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

double One_Dimensional_AGMND::get_A() {
	return (*left_point).second.num_estimation -
		m * ((*right_point).second.supposed_x2 - (*left_point).first);
}

double One_Dimensional_AGMND::get_B(double x) {
	return (*left_point).second.z +
		(*left_point).second.num_estimation * (x - (*left_point).first) -
		0.5 * m * pow(x - (*left_point).first, 2);
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

void One_Dimensional_AGMND::compute_R(double new_point, double new_m) {
	if (new_m != m) {
		while (!pq->empty())
			pq->pop();
		m = new_m;
		for (reset(); right_point != points->end(); go_Next_Interval()) {
			compute_supposed_x();
			(*right_point).second.R = get_R();
			pq->push(_interval({ (*left_point).first, (*left_point).second },
				{ (*right_point).first, (*right_point).second }));
		}
	}
	else {
		go_new_left_interval(new_point);
		for (int i = 0;i < 2;++i, go_Next_Interval()) {
			compute_supposed_x();
			(*right_point).second.R = get_R();
			pq->push(_interval({ (*left_point).first, (*left_point).second },
				{ (*right_point).first, (*right_point).second }));
		}
	}
}

double One_Dimensional_AGMND::get_m() {
	if (M_Max > 0)
		return r_p * M_Max;
	else if (M_Max == 0)
		return 1;
	else
		throw - 1;
}

void One_Dimensional_AGMND::insert_to_map(double _y, double _z, double _R,
	double _num_estimation) {
	_characteristics _ch(_z, _R, _num_estimation);
	points->insert({ _y, _ch });
	if (_z < res.z) {
		res.x = curr_x;
		res.y = _y;
		res.z = _z;
	}
}

void One_Dimensional_AGMND::compare_interval_len(double new_point) {
	go_new_left_interval(new_point);
	double interval_length;
	for (int i = 0;i < 2;++i, go_Next_Interval()) {
		interval_length = (*right_point).first - (*left_point).first;
		if (interval_length < min_interval_length)
			min_interval_length = interval_length;
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

double One_Dimensional_AGMND::get_new_point(_interval i) {
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

void One_Dimensional_AGMND::delete_containers() {
	delete points;
	points = nullptr;
	delete pq;
	pq = nullptr;
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
		pq->push(_interval({ (*left_point).first, (*left_point).second },
			{ (*right_point).first, (*right_point).second }));
	}
}

void One_Dimensional_AGMND::set_experiment (double _a, double _b,
	double _curr_x, double(*f)(double x, double y), double _eps, double _r_par) {
	res.k = 0;
	a = _a;
	b = _b;
	function = f;
	eps = _eps;
	r_p = _r_par;
	curr_x = _curr_x;
	if (points == nullptr)
		points = new std::map<double, _characteristics>;
	if (pq == nullptr)
		pq = new std::priority_queue<_interval, std::vector<_interval>, _CompareR>;
}

_result One_Dimensional_AGMND::get_result() {
	return res;
}

_result One_Dimensional_AGMND::solve() {
	std::pair<double, double> new_point;
	double new_m;

	perform_first_iteration(); // perform at the boundary points a and b 

	while (!isEnd()) {
		new_m = get_m();
		compute_R(new_point.first, new_m);
		new_point.first = get_new_point(pq->top()); pq->pop();
		new_point.second = (*function)(curr_x, new_point.first);

		if (new_point.first <= a || new_point.first >= b) {
			std::cout << "Out of border" << std::endl;
			throw - 1;
		}

		insert_to_map(new_point.first, new_point.second, 0, 0);
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	res.k = points->size();
	delete_containers();
	
	return res;
}

AGMND::AGMND(double _a, double _b, double _lower_y, double _upper_y,
	double(*f)(double x, double y),
	double _eps, double _r_par) :
	One_Dimensional_AGMND(_a, _b, 0, f, _eps, _r_par),
	lower_y(_lower_y), upper_y(_upper_y), function(f) {
	if (a > b)
		throw "b is a right border, must be more than a";
	if (lower_y > upper_y)
		throw "upper_y is a upper border, must be more than lower_y";
}

void AGMND::do_first_iteration(One_Dimensional_AGMND* odm,
	_result* tmp_res) {
	min_interval_length = b - a;
	odm->set_experiment(lower_y, upper_y, a, function, eps, r_p);

	*tmp_res = odm->solve();
	res.k += tmp_res->k;
	res.x = a;
	res.y = tmp_res->y;
	res.z = tmp_res->z;
	insert_to_map(res.x, res.y, res.z, 0, 0);

	if (a != b) {
		odm->set_experiment(lower_y, upper_y, b, function, eps, r_p);
		*tmp_res = odm->solve();
		res.k += tmp_res->k;
		insert_to_map(b, tmp_res->y, tmp_res->z, 0, 0);

		reset();
		(*left_point).second.num_estimation =
		(*right_point).second.num_estimation = get_num_estimation();
		M_Max = get_M();
		m = -1;
		pq->push(_interval({ (*left_point).first, (*left_point).second },
			{ (*right_point).first, (*right_point).second }));
	}
}

void AGMND::insert_to_map(double _x, double _y, double _z, double _R, double _num_estimation) {
	_characteristics _ch(_z, _R, _num_estimation);
	points->insert({ _x, _ch });
	if (_z < res.z) {
		res.x = _x;
		res.y = _y;
		res.z = _z;
	}
}

_result AGMND::get_result() {
	return res;
}

void AGMND::set_experiment(const double _a, const double _b,
	double _lower_y, double _upper_y, double(*f)(double x, double y),
	const double _eps, const double _r_par) {
	a = _a;
	b = _b;
	lower_y = _lower_y;
	upper_y = _upper_y;
	function = f;
	eps = _eps;
	r_p = _r_par;
	if (points == nullptr)
		points = new std::map<double, _characteristics>;
	if (pq == nullptr)
		pq = new std::priority_queue<_interval, std::vector<_interval>,
		_CompareR>;
}

void AGMND::solve() {
	One_Dimensional_AGMND odm(0, 0, 0, nullptr), *podm = &odm;
	std::pair<double, double> new_point;
	double new_m;
	_result tmp_res, *ptmp_res = &tmp_res;
	res.k = 0;
	do_first_iteration(podm, ptmp_res);
	
	while (min_interval_length > eps) {
		new_m = get_m();
		compute_R(new_point.first, new_m);
		new_point.first = get_new_point(pq->top()); pq->pop();
		odm.set_experiment(lower_y, upper_y, new_point.first, function, eps, r_p);
		tmp_res = odm.solve();
		res.k += tmp_res.k;
		new_point.second = tmp_res.z;
		insert_to_map(new_point.first, tmp_res.y, new_point.second, 0, 0);
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	delete_containers();
}
