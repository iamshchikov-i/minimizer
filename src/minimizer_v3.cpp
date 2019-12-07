#include "minimizer_v3.h"

void print_arr(std::vector<double> v, int procrank) {
	std::cout <<"procrank is "<<procrank<<"---";
	for (int i = 0; i < v.size(); ++i)
		std::cout << v[i] << " ";
	std::cout << std::endl;
}

int Minimizer_v3::get_procnum() {
	return procnum;
}

void Minimizer_v3::do_first_iteration() {
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

void Minimizer_v3::do_first_parallel_iteration(std::vector<double>& recvbuf, std::vector<double>& data, std::pair<double, double>& new_point) {
	double new_m = get_m();
	calculate_R(new_point.first, new_m);

	for (int i = 0; i < procnum; ++i) {
		if (procrank == i)
			new_point.first = get_new_point(pq->top());
		pq->pop();
	}
	new_point.second = (*function)(new_point.first);
	data[0] = new_point.first; data[1] = new_point.second;
	
	MPI_Allgather(data.data(), 2, MPI_DOUBLE, recvbuf.data(), 2, MPI_DOUBLE, MPI_COMM_WORLD);

	for (int i = 0; i < procnum; ++i) {
		insert_to_map(recvbuf[2 * i], recvbuf[(2 * i) + 1], 0);
		res.k++;
		compare_interval_len(recvbuf[2 * i]);
		compare_M(recvbuf[2 * i]);
	}
}

int Minimizer_v3::get_rank() {
	return procrank;
}

result Minimizer_v3::get_result() {
	return res;
}

Minimizer_v3::Minimizer_v3(double _a, double _b, double(*f)(double x),
	double _eps, int _N_max, double _r_par) :Minimizer_v2(_a, _b, f, _eps, _N_max, _r_par) {
	MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
}

void Minimizer_v3::set_experiment(const double _a, const double _b, double(*f)(double x), const double _eps,
	                              const int _N_max, const double _r_par) {
	a = _a;
	b = _b;
	function = f;
	eps = _eps;
	N_max = _N_max;
	r_p = _r_par;
	res.k = 0;
	if (values == nullptr)
		values = new std::map<double, characteristics>;
	if (pq == nullptr)
		pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
}

void Minimizer_v3::solve() {
	std::vector<double> recvbuf(procnum * 2), data(2); //x,y
	std::pair<double, double> new_point;
	double new_m;

	do_first_iteration();

	while (!isEnd() && pq->size() < procnum) {
		new_m = get_m();
		calculate_R(new_point.first, new_m);
		new_point.first = get_new_point(pq->top()); pq->pop();
		new_point.second = (*function)(new_point.first);
		insert_to_map(new_point.first, new_point.second, 0);
		res.k++;
		compare_interval_len(new_point.first);
		compare_M(new_point.first);
	}

	do_first_parallel_iteration(recvbuf, data, new_point);

	while (!isEnd()) {
		new_m = get_m();
		for (int i = 0; i < procnum; ++i)
			calculate_R(recvbuf[2 * i], new_m);

		for (int i = 0; i < procnum; ++i) {
			if (procrank == i)
				new_point.first = get_new_point(pq->top());
			pq->pop();
		}
		new_point.second = (*function)(new_point.first);
		data[0] = new_point.first; data[1] = new_point.second;

		MPI_Allgather(data.data(), 2, MPI_DOUBLE, recvbuf.data(), 2, MPI_DOUBLE, MPI_COMM_WORLD);

		for (int i = 0; i < procnum; ++i) {
			insert_to_map(recvbuf[2 * i], recvbuf[(2 * i) + 1], 0);
			compare_interval_len(recvbuf[2 * i]);
			compare_M(recvbuf[2 * i]);
			res.k++;
		}
	}
	delete_containers();
}
