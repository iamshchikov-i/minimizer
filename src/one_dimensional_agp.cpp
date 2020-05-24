#include "one_dimensional_agp.h"

One_Dimensional_AGP::One_Dimensional_AGP(double _a, double _b, double _curr_x,
	double(*f)(double x, double y),
	double _eps, double _r_par) : One_Dimensional_Minimizer(_a, _b, _curr_x, f, _eps, _r_par) {
	pq = new std::priority_queue<interval, std::vector<interval>, CompareR_min>;
}
bool One_Dimensional_AGP::isEnd() {
	return min_interval_length <= eps;
}

double One_Dimensional_AGP::get_M() {
    return std::abs(((*right_point).second.z - (*left_point).second.z) /
                    ((*right_point).first - (*left_point).first));
}

double One_Dimensional_AGP::get_m() {
    if (M_Max > 0) {
        return r_p * M_Max;
    } else if (M_Max == 0) {
        return 1;
    } else {
        throw - 1;
    }
}

double One_Dimensional_AGP::get_R() {
    double tmp = m * ((*right_point).first - (*left_point).first);
    return tmp + (std::pow((*right_point).second.z - (*left_point).second.z, 2)
    / tmp) - 2 * ((*right_point).second.z + (*left_point).second.z);
}


double One_Dimensional_AGP::get_new_point(interval i) {
	return 0.5*(i.first_point.first + i.second_point.first) -
		((i.second_point.second.z) - i.first_point.second.z) / (2 * m);
}

void One_Dimensional_AGP::compute_R(double new_point, double new_m) {
    if (new_m != m) {
		while (!pq->empty())
			pq->pop();
        m = new_m;
        reset();
        for (; right_point != points->end(); go_Next_Interval()) {
            (*right_point).second.R = get_R();
            pq->push(interval({ (*left_point).first, (*left_point).second },
                { (*right_point).first, (*right_point).second }));
        }
    } else {
        go_new_left_interval(new_point);
        for (int i = 0; i < 2; ++i, go_Next_Interval()) {
            (*right_point).second.R = get_R();
            pq->push(interval({ (*left_point).first, (*left_point).second },
                { (*right_point).first, (*right_point).second }));
        }
    }
}

void One_Dimensional_AGP::insert_to_map(double _y, double _z,
                                              double _R, double _num_estimation) {
    characteristics _ch(_z, _R, _num_estimation);
    points->insert({ _y, _ch });
    if (_z < res.z) {
		res.x = curr_x;
        res.y = _y;
        res.z = _z;
    }
}


void One_Dimensional_AGP::compare_M(double new_point) {
    go_new_left_interval(new_point);
    double M;
    for (int i = 0; i < 2; ++i, go_Next_Interval()) {
        M = get_M();
        if (M >= M_Max)
            M_Max = M;
    }
}

void One_Dimensional_AGP::perform_first_iteration() {
	min_interval_length = b - a;
	res.x = curr_x;
	res.y = a; res.z = (*function)(curr_x, a);
	if (a != b) {
		insert_to_map(res.y, res.z, 0);
		insert_to_map(b, (*function)(curr_x, b), 0);
		reset();
		M_Max = get_M();
		m = -1;
		pq->push(interval({ (*left_point).first, (*left_point).second },
			{ (*right_point).first, (*right_point).second }));
	}
}

void One_Dimensional_AGP::delete_containers() {
	delete points;
	points = nullptr;
	delete pq;
	pq = nullptr;
}

void One_Dimensional_AGP::set_experiment(double _a, double _b, double _curr_x,
	double(*f)(double x, double y),
	double _eps, double _r_par) {
    if (a > b)
        throw "b is a right border, must be more than a";
	res.k = 0;
    a = _a;
    b = _b;
	eps = _eps;
	r_p = _r_par;
    curr_x = _curr_x;
    function = f;
    if (points == nullptr)
        points = new std::map<double, characteristics>;
    if (pq == nullptr)
        pq = new std::priority_queue<interval, std::vector<interval>,
                                     CompareR_min>;
}

result One_Dimensional_AGP::solve() {
	std::pair<double, double> new_point;
	double new_m;

	perform_first_iteration();

	while (min_interval_length > eps) {
		new_m = get_m();
		compute_R(new_point.first, new_m);
		new_point.first = get_new_point(pq->top()); pq->pop();
		new_point.second = (*function)(curr_x, new_point.first);
		insert_to_map(new_point.first, new_point.second, 0);
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	res.k = points->size();
	delete_containers();

	return res;
}