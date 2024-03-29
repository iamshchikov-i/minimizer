#include "multi_dimensional_minimizer.h"

Multi_Dimensional_Minimizer::Multi_Dimensional_Minimizer(int _range, std::vector<double>& _lower_bound,
	std::vector<double>& _upper_bound,
	double(*f)(std::vector<double> coords), bool _useMPI, bool _useThreads, int _threadsNum, 
	Upper_method upper_method, double _eps, int _Nmax,
	double _r_par) : range(_range), lower_bound(_lower_bound),
	upper_bound(_upper_bound), function(f), useMPI(_useMPI), eps(_eps), Nmax(_Nmax), r_par(_r_par) {
	One_Dimensional_Minimizer* podm;

	std::vector<std::pair<double, double>> bounds;
	std::vector<double> curr_x;
	for (int i = 0; i < range; ++i) {
		bounds.push_back({ lower_bound[i], upper_bound[i] });
		curr_x.push_back(lower_bound[i]);
		if (i == range - 1) {
			switch (upper_method) {
			case Upper_method::AGP :
				podm = new One_Dimensional_AGP(range, i, odm, bounds, curr_x,
					_useMPI, _useThreads, _threadsNum,
					function, eps, Nmax, r_par);
				odm.push_back(podm);
				break;
			case Upper_method::AGMND :
				podm = new One_Dimensional_AGMND(range, i, odm, bounds, curr_x,
					_useMPI, _useThreads, _threadsNum,
					function, eps, Nmax, r_par);
				odm.push_back(podm);
				break;
			}
		} else {
			One_Dimensional_AGP* podm = new One_Dimensional_AGP(range, i, odm, bounds, curr_x,
				_useMPI, _useThreads, _threadsNum,
				function, eps, Nmax, r_par);
			odm.push_back(podm);
		}
	}
	for (int i = 0; i < range; ++i)
		odm[i]->set_experiment(range, i, odm, bounds, curr_x, 
			_useMPI, _useThreads, _threadsNum, function, eps, Nmax, r_par);
}


void Multi_Dimensional_Minimizer::set_experiment(int _range, std::vector<double>& _lower_bound,
	std::vector<double>& _upper_bound,
	double(*f)(std::vector<double> coords), double _eps,
	double _r_par) {
	
	range = _range;
	lower_bound = _lower_bound;
	upper_bound = _upper_bound;
	function = f;
	eps = _eps;
	r_par = _r_par;
}

result Multi_Dimensional_Minimizer::solve() {
	result res;

	if(useMPI)
		res = odm[0]->solve_mpi();
	else
		res = odm[0]->solve_seq();

	return res;
}