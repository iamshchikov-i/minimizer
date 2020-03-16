#include "two_dimensional_minimizer.h"


Two_Dimensional_Minimizer::Two_Dimensional_Minimizer(One_Dimensional_Minimizer* _odm, double _a, double _b, double _lower_y, double _upper_y,
	double(*f)(double x, double y), double _eps, double _r_par) :
	One_Dimensional_AGP(_a, _b, 0, f, _eps, _r_par),
	odm(_odm), lower_y(_lower_y), upper_y(_upper_y), function(f) {
	if (a > b)
		throw "b is a right border, must be more than a";
	if (lower_y > upper_y)
		throw "upper_y is a upper border, must be more than lower_y";
}

void Two_Dimensional_Minimizer::perform_first_iteration(result* tmp_res) {
	min_interval_length = b - a;
	odm->set_experiment(lower_y, upper_y, a, function, r_p);

	*tmp_res = odm->solve();
	res.k += tmp_res->k;
	res.k_max_on_y = tmp_res->k;
	res.x = a;
	res.y = tmp_res->y;
	res.z = tmp_res->z;
	insert_to_map(res.x, res.y, res.z, 0);

	if (a != b) {
		odm->set_experiment(lower_y, upper_y, b, function, r_p);
		*tmp_res = odm->solve();
		res.k += tmp_res->k;

		if (tmp_res->k > res.k_max_on_y)
			res.k_max_on_y = tmp_res->k;

		insert_to_map(b, tmp_res->y, tmp_res->z, 0);

		reset();
		M_Max = get_M();
		m = -1;
		pq->push(interval({ (*left_point).first, (*left_point).second },
				 { (*right_point).first, (*right_point).second }));
	}
}

void Two_Dimensional_Minimizer::insert_to_map(double _x, double _y, double _z, double _R) {
	characteristics ch(_z, _R, 0);
	points->insert({ _x, ch });
	if (_z < res.z) {
		res.x = _x;
		res.y = _y;
		res.z = _z;
	}
}

void Two_Dimensional_Minimizer::set_experiment(One_Dimensional_Minimizer* _odm, const double _a, const double _b,
	double _lower_y, double _upper_y, double(*f)(double x, double y),
	const double _eps, const double _r_par) {
	odm = _odm;
	b = _b;
	lower_y = _lower_y;
	upper_y = _upper_y;
	function = f;
	eps = _eps;
	r_p = _r_par;
	if (points == nullptr)
		points = new std::map<double, characteristics>;
	if (pq == nullptr)
		pq = new std::priority_queue<interval, std::vector<interval>,
									CompareR_min>;
}

result Two_Dimensional_Minimizer::solve() {
	std::pair<double, double> new_point;
	double new_m;
	result tmp_res, *ptmp_res = &tmp_res;
	res.k = 0;
	perform_first_iteration(ptmp_res);

	while (min_interval_length > eps) {
		new_m = get_m();
		compute_R(new_point.first, new_m);
		new_point.first = get_new_point(pq->top()); pq->pop();
		odm->set_experiment(lower_y, upper_y, new_point.first, function, r_p);
		tmp_res = odm->solve();
		res.k += tmp_res.k;

		if (tmp_res.k > res.k_max_on_y)
			res.k_max_on_y = tmp_res.k;

		new_point.second = tmp_res.z;
		insert_to_map(new_point.first, tmp_res.y, new_point.second, 0);
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	res.k_on_x = points->size();
	delete_containers();

	return res;
}