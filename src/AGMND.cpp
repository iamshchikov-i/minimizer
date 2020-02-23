#include "AGMND.h"

characteristics::characteristics() {}

characteristics::characteristics(double _y, double _R, double _num_estimation): 
							     y(_y), R(_R), num_estimation(_num_estimation) {}

characteristics::characteristics(double _R): R(_R) {}

interval::interval() {}

interval::interval(std::pair<double, characteristics> _f_point,
	               std::pair<double, characteristics> _s_point) :
				   first_point(_f_point), second_point(_s_point) {}

Minimizer::Minimizer(double _a, double _b, double(*f)(double x),
	                 double _eps, double _r_par): a(_a), b(_b), function(f),
					 eps(_eps), r_p(_r_par) {
	points = new std::map<double, characteristics>;
	pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
}

double Minimizer::get_num_estimation() {
	return ((*right_point).second.y - (*left_point).second.y) /
		   ((*right_point).first - (*left_point).first);
}

bool Minimizer::isEnd() {
	return min_interval_length <= eps;
}

void Minimizer::go_Next_Interval() {
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

void Minimizer::compute_R(double new_point, double new_m) {
	double supposed_x1, supposed_x2, supposed_x3;
	if (new_m != m) {
		if (!pq->empty()) {
			delete pq;
			pq = new std::priority_queue<interval, std::vector<interval>, 
				                         CompareR>;
		}
		m = new_m;
		for (reset(); right_point != points->end(); go_Next_Interval()) {
			compute_supposed_x();
			(*right_point).second.R = get_R();
			pq->push(interval({ (*left_point).first, (*left_point).second },
				{ (*right_point).first, (*right_point).second }));
		}
	}
	else {
		go_new_left_interval(new_point);
		for (int i = 0;i < 2;++i, go_Next_Interval()) {
			compute_supposed_x();
			(*right_point).second.R = get_R();
			pq->push(interval({ (*left_point).first, (*left_point).second }, 
				{ (*right_point).first, (*right_point).second }));
		}
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
	characteristics _ch(_y, _R, _num_estimation);
	points->insert({ _x, _ch });
	if (_y < res.y) {
		res.x = _x;
		res.y = _y;
	}
}

void Minimizer::compare_interval_len(double new_point) {
	go_new_left_interval(new_point);
	double interval_length;
	for (int i = 0;i < 2;++i, go_Next_Interval()) {
		interval_length = (*right_point).first - (*left_point).first;
		if (interval_length < min_interval_length)
			min_interval_length = interval_length;
	}
}

void Minimizer::compare_M(double new_point) {
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

double Minimizer::get_new_point(interval i) {
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
	double num_estimation;
	min_interval_length = b - a;
	res.x = a; res.y = (*function)(a);
	insert_to_map(res.x, res.y, 0, 0);
	insert_to_map(b, (*function)(b), 0, 0);
	reset();
	(*left_point).second.num_estimation =
	(*right_point).second.num_estimation = get_num_estimation();
	M_Max = get_M();
	m = -1;
	pq->push(interval({ (*left_point).first, (*left_point).second }, 
		              { (*right_point).first, (*right_point).second }));
}

void Minimizer::set_experiment(double _a, double _b, double(*f)(double x),
	double _eps, double _r_par) {
	a = _a;
	b = _b;
	function = f;
	eps = _eps;
	r_p = _r_par;
	if(points == nullptr)
		points = new std::map<double, characteristics>;
	if(pq == nullptr)
		pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
}

result Minimizer::get_result() {
	return res;
}

result Minimizer::solve() {
	std::pair<double, double> new_point;
	double new_m, M;

	perform_first_iteration(); // perform at the boundary points a and b 
	
	while (!isEnd()) {
		new_m = get_m();
		compute_R(new_point.first, new_m);
		new_point.first = get_new_point(pq->top()); pq->pop();
		new_point.second = (*function)(new_point.first);
		insert_to_map(new_point.first, new_point.second, 0, 0);
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	delete_containers();
	return res;
}