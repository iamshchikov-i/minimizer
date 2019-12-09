#include "minimizer_v4.h"

void Minimizer_v4::prepare_data_for_parallel_job(std::vector<interval>& buf, std::pair<double, double>& new_point) {
	for (int i = 0; i < th_num; ++i) {
		buf[i] = pq->top();
		pq->pop();
	}
}

void Minimizer_v4::do_shared_job(std::vector<interval>& buf, std::vector<double>& data, int th_rank) {
	std::pair<double, double> new_point;
	new_point.first = get_new_point(buf[th_rank]);
	new_point.second = (*function)(new_point.first);
	data[2 * th_rank] = new_point.first;
	data[(2 * th_rank) + 1] = new_point.second;
}

void Minimizer_v4::insert_data_to_map(std::vector<double>& data) {
	for (int i = 0; i < th_num; ++i) {
		insert_to_map(data[2 * i], data[(2 * i) + 1], 0);
		res.k++;
		compare_interval_len(data[2 * i]);
		compare_M(data[2 * i]);
	}
}

int Minimizer_v4::get_procnum() {
	return th_num;
}

void Minimizer_v4::do_first_iteration() {
	min_interval_length = b - a;
	res.x = a; res.y = (*function)(a);
	insert_to_map(res.x, res.y, 0);
	insert_to_map(b, (*function)(b), 0);
	reset();

	M_Max = get_M();
	m = -1;
	pq->push(interval({ (*left_point).first, (*left_point).second.y }, { (*right_point).first, (*right_point).second.y }, 0));
	res.k = 2;
}

void Minimizer_v4::do_first_parallel_iteration(std::thread *pth, std::vector<interval>& buf,
											   std::vector<double>& data, std::pair<double, double>& new_point, int rank) {
	if (rank == 0) {
		double new_m = get_m();
		calculate_R(new_point.first, new_m);
		prepare_data_for_parallel_job(buf, new_point);
		//sem = 1;
	}
	else {
		while (sem == 0) {}
	}
		
	do_shared_job(buf, data, rank);
	/*++sem;
	while (sem != 0 && sem != th_num + 1) { }
	sem = 0;*/
	
	if (rank == 0)
		insert_data_to_map(data);
}


result Minimizer_v4::get_result() {
	return res;
}

Minimizer_v4::Minimizer_v4(double _a, double _b, double(*f)(double x), int _th_num,
	double _eps, int _N_max, double _r_par) : Minimizer_v2(_a, _b, f, _eps, _N_max, _r_par) {
	if (_th_num <= 0)
		throw "number of threads must be more than zero";
	th_num = _th_num;
	sem = 0;
}

void Minimizer_v4::set_experiment(const double _a, const double _b, double(*f)(double x), int _th_num, const double _eps,
	const int _N_max, const double _r_par) {
	a = _a;
	b = _b;
	function = f;
	if (_th_num <= 0)
		throw "number of threads must be more than zero";
	th_num = _th_num;
	eps = _eps;
	N_max = _N_max;
	r_p = _r_par;
	res.k = 0;
	sem = 0;
	if (values == nullptr)
		values = new std::map<double, characteristics>;
	if (pq == nullptr)
		pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
}

void Minimizer_v4::do_parallel_job(std::thread *pth, std::vector<interval>& buf, std::vector<double>& data, 
                                   std::pair<double, double>& new_point, int rank) {
	double new_m;
	do_first_parallel_iteration(pth, buf, data, new_point, rank);

	while (!isEnd()) {
		if (rank == 0) {
			new_m = get_m();

			for (int i = 0; i < th_num; ++i)
				calculate_R(data[2 * i], new_m);

			prepare_data_for_parallel_job(buf, new_point);
			//sem = 1;
		}
		else
			while (sem == 0) {}

		do_shared_job(buf, data, rank);
		/*sem++;
		while (sem != 0 && sem != th_num + 1) {}
		sem = 0;*/

		if(rank == 0)
			insert_data_to_map(data);
	}
	if(rank == 0)
		delete_containers();
}

void Minimizer_v4::solve() {
	std::thread *pth = new std::thread[th_num];
	std::vector<interval> buf(th_num);
	std::vector<double> data(th_num * 2);
	std::pair<double, double> new_point;
	double new_m;
	
	do_first_iteration();
	
	while (!isEnd() && pq->size() < th_num) {
		new_m = get_m();
		calculate_R(new_point.first, new_m);
		new_point.first = get_new_point(pq->top()); pq->pop();
		new_point.second = (*function)(new_point.first);
		insert_to_map(new_point.first, new_point.second, 0);
		res.k++;
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}
	
	//require synchronization! (now it's not works)
	for (int i = 1; i < th_num; ++i)
		pth[i] = std::thread(&Minimizer_v4::do_parallel_job, this, pth, std::ref(buf), std::ref(data), new_point, i);
	do_parallel_job(pth, buf, data, new_point, 0);
	
	for (int i = 1; i < th_num; ++i)
		pth[i].join();
}
