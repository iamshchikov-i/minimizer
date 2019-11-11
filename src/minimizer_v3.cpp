#include "minimizer_v3.h"

void print_arr(std::vector<double> v, int procrank) {
	std::cout <<"procrank is "<<procrank<<"---";
	for (int i = 0; i < v.size(); ++i)
		std::cout << v[i] << " ";
	std::cout << std::endl;
}

void Minimizer_v3::calc_involved_procnum(int queue_size) {
	if (procnum <= queue_size)
		curr_involved_procnum = procnum;
	else
		curr_involved_procnum = queue_size;
}

void Minimizer_v3::set_comm() {
	MPI_Group group_world, new_work_group, new_wait_group;
	std::vector<int> member_work, member_wait;

	if (curr_involved_procnum < procnum) {
		MPI_Comm_group(MPI_COMM_WORLD, &group_world);
		member_wait.push_back(0);
		for (int i = 0; i < curr_involved_procnum; ++i) 
			member_work.push_back(i);
		for (int i = 0; i < procnum - curr_involved_procnum; ++i) {
			member_wait.push_back(procnum - 1 - i);
		}
			
		MPI_Group_incl(group_world, member_work.size(), member_work.data(),
			&new_work_group);
		MPI_Group_incl(group_world, member_wait.size(), member_wait.data(),
			&new_wait_group);
		MPI_Comm_create(MPI_COMM_WORLD, new_work_group, &involved_comm);
		MPI_Comm_create(MPI_COMM_WORLD, new_wait_group, &waiting_comm);
	}
	else
		involved_comm = MPI_COMM_WORLD;
}

Minimizer_v3::Minimizer_v3(double _a, double _b, double(*f)(double x),
	double _eps, int _N_max, double _r_par) :Minimizer_v2(_a, _b, f, _eps, _N_max, _r_par) {
	MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	involved_comm = MPI_COMM_NULL;
	waiting_comm = MPI_COMM_NULL;
	curr_involved_procnum;
}

result Minimizer_v3::solve() {
	std::vector<double> recvbuf, data(2); //x,y
	std::pair<double, double> new_point;
	double new_m;

	min_interval_length = b - a;
	insert_to_map(a, (*function)(a), 0);
	res.x = a; res.y = (*function)(a);
	insert_to_map(b, (*function)(b), 0);
	reset();
	M_Max = get_M();
	m = -1;
	pq->push(interval(a, b, 0));
	res.k = 2;
	recvbuf.resize(curr_involved_procnum * 2);
	recvbuf[0] = res.x; recvbuf[1] = res.y;

	while (!isEnd()) {
		new_m = get_m();
		for (int i = 0; i < curr_involved_procnum; ++i)
			calculate_R(recvbuf[2 * i], new_m);
		calc_involved_procnum(pq->size());
		set_comm();
		recvbuf.resize(curr_involved_procnum * 2);
		
		if (procrank < curr_involved_procnum) {
			for (int i = 0; i < curr_involved_procnum; ++i) {
				if (procrank == i)
					new_point.first = get_new_point(pq->top());
				pq->pop();
			}
			new_point.second = (*function)(new_point.first);
			data[0] = new_point.first; data[1] = new_point.second;

			MPI_Allgather(data.data(), 2, MPI_DOUBLE, recvbuf.data(), 2, MPI_DOUBLE, involved_comm);
		}
		
		if (curr_involved_procnum != procnum) {
			MPI_Bcast(recvbuf.data(), curr_involved_procnum * 2, MPI_DOUBLE, 0, waiting_comm);
		}

		for (int i = 0; i < curr_involved_procnum; ++i) {
			insert_to_map(recvbuf[2*i], recvbuf[(2*i) + 1], 0);/**/
			compare_interval_len(recvbuf[2 * i]);
			compare_M(recvbuf[2 * i]);
			res.k++;
		}
	}
	delete_containers();
	return res;
}
