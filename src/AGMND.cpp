#include "AGMND.h"

_characteristics::_characteristics() {}

_characteristics::_characteristics(double _y, double _R, double _num_estimation): 
							     y(_y), R(_R), num_estimation(_num_estimation) {}

_characteristics::_characteristics(double _R): R(_R) {}

_interval::_interval() {}

_interval::_interval(std::pair<double, _characteristics> _f_point,
	               std::pair<double, _characteristics> _s_point) :
				   first_point(_f_point), second_point(_s_point) {}

Minimizer::Minimizer(double _a, double _b, double(*f)(double x),
	                 double _eps, double _r_par): a(_a), b(_b), function(f),
					 eps(_eps), r_p(_r_par) {
	points = new std::map<double, _characteristics>;
	pq = new std::priority_queue<_interval, std::vector<_interval>, CompareR_max >;
	recalc = false;
}

void Minimizer::compute_num_estimation() {
	reset();
	std::map<double, _characteristics>::iterator point_after_right = right_point;
	point_after_right++;

	double H, h1, h2, d;

	h1 = (*right_point).first - (*left_point).first;
	h2 = (*point_after_right).first - (*right_point).first;
	d = h2 / h1;
	H = h1 + h2;

	(*left_point).second.num_estimation = ( -(2 + d) * (*left_point).second.y + 
		(pow((1 + d), 2) * (*right_point).second.y / d) - (*point_after_right).second.y / d) / H;

	for (; point_after_right != --points->end(); go_Next_interval(), point_after_right++) {
		h1 = (*right_point).first - (*left_point).first;
		h2 = (*point_after_right).first - (*right_point).first;
		d = h2 / h1;
		H = h1 + h2;

		(*right_point).second.num_estimation = (-d * (*left_point).second.y +
			((pow(d, 2) - 1) * (*right_point).second.y / d) + (*point_after_right).second.y / d) / H;
	}

	h1 = (*right_point).first - (*left_point).first;
	h2 = (*point_after_right).first - (*right_point).first;
	d = h2 / h1;
	H = h1 + h2;

	(*point_after_right).second.num_estimation = (d * (*left_point).second.y -
		( pow((1 + d), 2) * (*right_point).second.y / d) + (2 * d + 1) * (*point_after_right).second.y / d) / H;
}

bool Minimizer::isEnd() {
	return min_interval_length <= eps;
}

void Minimizer::go_Next_interval() {
	left_point++; right_point++;
}

void Minimizer::go_new_left_interval(double new_point) {
	right_point = points->find(new_point);
	left_point = right_point; left_point--;
}

void Minimizer::reset() {
	left_point = points->begin();
	right_point = left_point; right_point++;
}

double Minimizer::get_M() {
	double t1, t2, t3, dx, dy;

	dx = (*right_point).first - (*left_point).first;
	dy = (*right_point).second.y - (*left_point).second.y;
	
	t1 = abs((*right_point).second.num_estimation - 
		     (*left_point).second.num_estimation) / abs(dx);
	t2 = -2 * (dy - (*left_point).second.num_estimation * dx) / pow(dx, 2);
	t3 = 2 * (dy - (*right_point).second.num_estimation * dx) / pow(dx, 2);
	
	return std::max({ t1, t2, t3 });
}

double Minimizer::get_d() {
	return ((*right_point).first - (*left_point).first) / 2 -
		   ((*right_point).second.num_estimation - 
		   (*left_point).second.num_estimation) / (2 * m);
}

void Minimizer::compute_supposed_x() {
	double tmp1, tmp2, tmp3;
	tmp1 = ((*left_point).second.y - (*left_point).second.num_estimation * 
		   (*left_point).first) - ((*right_point).second.y - 
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

double Minimizer::get_A() {
	return (*left_point).second.num_estimation - 
		   m * ((*right_point).second.supposed_x2 - (*left_point).first);
}

double Minimizer::get_B(double x) {
	return (*left_point).second.y +
		   (*left_point).second.num_estimation * (x - (*left_point).first) -
		   0.5 * m * pow(x - (*left_point).first, 2);
}

double Minimizer::auxiliary_function(double x) {
	if (x >= (*right_point).second.supposed_x2 &&
		x <= (*right_point).second.supposed_x3) { 
		
		return get_A() * (x - (*right_point).second.supposed_x2) +
			   0.5 * m * pow(x - (*right_point).second.supposed_x2, 2) +
			   get_B((*right_point).second.supposed_x2);
	}
		
	else if (x > (*left_point).first &&
		     x < (*right_point).second.supposed_x2) {
		
		return (*left_point).second.y + (*left_point).second.num_estimation *
			   (x - (*left_point).first) - 0.5 * m * pow(x - (*left_point).first, 2);
	}

	else if (x > (*right_point).second.supposed_x2 &&
			 x <= (*right_point).first) {

		return (*right_point).second.y - (*right_point).second.num_estimation *
			   (x - (*right_point).first) -
			   0.5 * m * pow(x - (*right_point).first, 2);
	}
	else throw - 1;
}

double Minimizer::get_R() {
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

void Minimizer::check_supposed_x() {
	double _x = (*right_point).second.supposed_x2, x_ = (*right_point).second.supposed_x3,
		x1 = (*left_point).first, x2 = (*right_point).first,
		z1 = (*left_point).second.y, z2 = (*left_point).second.y,
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

void Minimizer::recalc_characteristics() {
	do {
		recalc = false;
		while (!pq->empty())
			pq->pop();
		for (reset(); right_point != points->end(); go_Next_interval()) {
			compute_supposed_x();
			check_supposed_x();
			if (recalc) break;
			(*right_point).second.R = get_R();
			pq->push(_interval({ (*left_point).first, (*left_point).second },
				{ (*right_point).first, (*right_point).second }));
		}
	} while (right_point != points->end());
}

void Minimizer::compute_R(double new_point, double new_m) {
	if (new_m != m) {
		m = new_m;
		recalc_characteristics();
	}
	else {
		go_new_left_interval(new_point);
		for (int i = 0;i < 2;++i, go_Next_interval()) {
			compute_supposed_x();
			check_supposed_x();
			if (recalc) break;
			(*right_point).second.R = get_R();
			pq->push(_interval({ (*left_point).first, (*left_point).second },
				{ (*right_point).first, (*right_point).second }));
		}

		if (recalc) recalc_characteristics();
	}
}

double Minimizer::get_m() {
	if (M_Max > 0)
		return r_p * M_Max;
	else if (M_Max == 0)
		return 1;
	else
		throw - 1;
}

void Minimizer::insert_to_map(double _x, double _y, double _R, 
							  double _num_estimation) {
	_characteristics _ch(_y, _R, _num_estimation);
	points->insert({ _x, _ch });
	if (_y < res.y) {
		res.x = _x;
		res.y = _y;
	}
}

void Minimizer::compare_interval_len(double new_point) {
	go_new_left_interval(new_point);
	double _interval_length;
	for (int i = 0;i < 2;++i, go_Next_interval()) {
		_interval_length = (*right_point).first - (*left_point).first;
		if (_interval_length < min_interval_length)
			min_interval_length = _interval_length;
	}
}

void Minimizer::compare_M(double new_point) {
	double M;
	go_new_left_interval(new_point);
	for (int i = 0;i < 2;++i, go_Next_interval()) {
		M = get_M();
		if (M >= M_Max)
			M_Max = M;
	}
}

double Minimizer::get_new_point(_interval i) {
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

void Minimizer::delete_containers() {
	delete points;
	points = nullptr;
	delete pq;
	pq = nullptr;
}

void Minimizer::perform_first_iteration() {
	min_interval_length = b - a;
	res.x = a;
	res.y = (*function)(a);

	if (a != b) {
		insert_to_map(res.x, res.y, 0, 0);
		insert_to_map(b, (*function)(b), 0, 0);
		reset();
		(*left_point).second.num_estimation =
		(*right_point).second.num_estimation = ((*right_point).second.y - (*left_point).second.y) /
											   ((*right_point).first - (*left_point).first);
		M_Max = get_M();
		m = -1;
		pq->push(_interval({ (*left_point).first, (*left_point).second },
			{ (*right_point).first, (*right_point).second }));
	}
}

void Minimizer::set_experiment(double _a, double _b, double(*f)(double x),
	double _eps, double _r_par) {
	a = _a;
	b = _b;
	function = f;
	eps = _eps;
	r_p = _r_par;
	if(points == nullptr)
		points = new std::map<double, _characteristics>;
	if(pq == nullptr)
		pq = new std::priority_queue<_interval, std::vector<_interval>, CompareR_max >;
	recalc = false;
}

_result Minimizer::get_result() {
	return res;
}

_result Minimizer::solve() {
	//eps *= b - a;
	std::pair<double, double> new_point;
	double new_m;
	_interval search_interval;

	perform_first_iteration(); // perform at the boundary points a and b 

	while (!isEnd()) {
		new_m = get_m();
		compute_R(new_point.first, new_m);
		search_interval = pq->top(); pq->pop();
		new_point.first = get_new_point(search_interval);

		new_point.second = (*function)(new_point.first);
		insert_to_map(new_point.first, new_point.second, 0, 0);
		compute_num_estimation();
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	res.k = points->size();
	delete_containers();

	return res;
}