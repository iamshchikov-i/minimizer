#include "one_dimensional_agp.h"

One_Dimensional_AGP::One_Dimensional_AGP(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
	std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
	bool _useMPI, bool _useThreads, int _threadsNum,
	double(*f)(std::vector<double> x),
	double _eps, int _Nmax, double _r_par) :
	One_Dimensional_Minimizer(_range, _curr_dim, _odm, _bounds, _curr_x, 
		_useMPI, _useThreads, _threadsNum, 
		f,
		_eps, _Nmax, _r_par) {
	pq = new std::priority_queue<interval, std::vector<interval>, CompareR_min>;
}
bool One_Dimensional_AGP::isEnd() {
	return min_interval_length <= eps || points->size() > Nmax;
}

double One_Dimensional_AGP::get_M() {
    return std::abs(((*right_point).second.z - (*left_point).second.z) /
                    ((*right_point).first[curr_dim] - (*left_point).first[curr_dim]));
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
    double tmp = m * ((*right_point).first[curr_dim] - (*left_point).first[curr_dim]);
    return tmp + (std::pow((*right_point).second.z - (*left_point).second.z, 2)
    / tmp) - 2 * ((*right_point).second.z + (*left_point).second.z);
}


double One_Dimensional_AGP::get_new_point(interval i) {
	return 0.5*(i.first_point.first[curr_dim] + i.second_point.first[curr_dim]) -
		((i.second_point.second.z) - i.first_point.second.z) / (2 * m);
}

void One_Dimensional_AGP::compute_R(std::vector<double> new_point, double new_m) {
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

void One_Dimensional_AGP::insert_to_map(std::vector<double> coords, double _z, double _R,
	double _num_estimation) {
    characteristics _ch(_z, _R, _num_estimation);
    points->insert({ coords, _ch });
    if (_z < res.z) {
		res.coords = coords;
        res.z = _z;

		//print(curr_x);
		//std::cout << res.z << std::endl;
    }
}


void One_Dimensional_AGP::compare_M(std::vector<double> new_point) {
    go_new_left_interval(new_point);
    double M;
    for (int i = 0; i < 2; ++i, go_Next_Interval()) {
        M = get_M();
        if (M >= M_Max)
            M_Max = M;
    }
}

void One_Dimensional_AGP::perform_first_iteration() {
	double a = bounds[curr_dim].first, b = bounds[curr_dim].second;
	min_interval_length = b - a;
	curr_x[curr_dim] = a;
	res.coords = curr_x;
	result tmp_res;

	if (curr_dim != range - 1) {
		odm[curr_dim + 1]->set_experiment(range, curr_dim + 1, odm, bounds, curr_x,
			useMPI, useThreads, threadsNum,
			function, eps, Nmax, r_p);
		if(useMPI)
			odm[curr_dim + 1]->solve_mpi();
		else
			odm[curr_dim + 1]->solve_seq();
		tmp_res = odm[curr_dim + 1]->get_result();
		res.coords = tmp_res.coords;
		res.z = tmp_res.z;

		for(int i = curr_dim + 1; i < range; ++i)
			res.k[i] += tmp_res.k[i];
	} else {
		tmp_res.coords = curr_x;
		res.z = (*function)(res.coords);
		tmp_res.z = res.z;
	}
	
	if (a != b) {
		insert_to_map(tmp_res.coords, tmp_res.z, 0);

		curr_x[curr_dim] = b;
		if (curr_dim != range - 1) {
			odm[curr_dim + 1]->set_experiment(range, curr_dim + 1, odm, bounds, curr_x,
				useMPI, useThreads, threadsNum, 
				function, eps, Nmax, r_p);
			if(useMPI)
				odm[curr_dim + 1]->solve_mpi();
			else
				odm[curr_dim + 1]->solve_seq();
			tmp_res = odm[curr_dim + 1]->get_result();

			for (int i = curr_dim + 1; i < range; ++i)
				res.k[i] += tmp_res.k[i];
		} else {
			tmp_res.coords = curr_x;
			tmp_res.z = (*function)(curr_x);
		}
		
		insert_to_map(tmp_res.coords, tmp_res.z, 0);
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

void One_Dimensional_AGP::set_experiment(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
	std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
	bool _useMPI, bool _useThreads, int _threadsNum,
	double(*f)(std::vector<double> x),
	double _eps, int _Nmax, double _r_par) {

	res.k = std::vector<int>(range, 0);
	range = _range;
	curr_dim = _curr_dim;
	odm = _odm;
	bounds = _bounds;
	curr_x = _curr_x;
	function = f;
	eps = _eps;
	Nmax = _Nmax;
	r_p = _r_par;

	useMPI = _useMPI;
	useThreads = _useThreads;
	threadsNum = _threadsNum;
    
    if (points == nullptr)
        points = new std::map<std::vector<double>, characteristics>;
    if (pq == nullptr)
        pq = new std::priority_queue<interval, std::vector<interval>,
                                     CompareR_min>;
}

void One_Dimensional_AGP::perform_first_parallel_step(std::pair<double,
	double>& new_point, result& tmp_res, double& new_m, 
	std::map<int, std::pair<double, double>>& processing_intervals_map, 
	std::set<std::pair<double, double>>& processing_intervals_set, MPI_Status& st) {
	new_m = get_m();
	compute_R(curr_x, new_m);

	if (procrank == 0) {
		for (int i = 1; i < procnum; ++i) {
			auto interval = pq->top();
			processing_intervals_map.insert({i, {interval.first_point.first[0], interval.second_point.first[0]}});
			processing_intervals_set.insert({ interval.first_point.first[0], interval.second_point.first[0] });
			pq->pop();
		}

		std::vector<double> data(range + 1);
		std::vector<int> k(range);
		MPI_Recv(data.data(), data.size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
		processing_intervals_set.erase(processing_intervals_map[st.MPI_SOURCE]);
		processing_intervals_map.erase(st.MPI_SOURCE);
		

		tmp_res.z = data.back();
		data.pop_back();
		tmp_res.coords = data;
		curr_x = tmp_res.coords;

		insert_to_map(tmp_res.coords, tmp_res.z, 0);
		compare_interval_len(curr_x);
		compare_M(curr_x);

		MPI_Recv(k.data(), k.size(), MPI_INT, st.MPI_SOURCE, st.MPI_SOURCE, MPI_COMM_WORLD, &st);
		for (int i = curr_dim + 1; i < range; ++i)
			res.k[i] += k[i];
		
	} else {
		for (int i = 1; i < procnum; ++i) {
			if (procrank == i) {
				new_point.first = get_new_point(pq->top());
				break;
			}
			pq->pop();
		}

		curr_x[curr_dim] = new_point.first;
		
		odm[curr_dim + 1]->set_experiment(range, curr_dim + 1, odm, bounds, curr_x, 
			useMPI, useThreads, threadsNum, 
			function, eps, Nmax, r_p);
		odm[curr_dim + 1]->solve_mpi();
		tmp_res = odm[curr_dim + 1]->get_result();

		for (int i = curr_dim + 1; i < range; ++i)
			res.k[i] += tmp_res.k[i];
		

		std::vector<double> data(tmp_res.coords);
		data.push_back(tmp_res.z);
		MPI_Send(data.data(), data.size(), MPI_DOUBLE, 0, procrank, MPI_COMM_WORLD);
		MPI_Send(tmp_res.k.data(), tmp_res.k.size(), MPI_INT, 0, procrank, MPI_COMM_WORLD);
	}
}

void One_Dimensional_AGP::do_parallel_job(double last_coord,
	std::vector<double>& res, int rank) {

	std::vector<double> tmp_coords(curr_x);
	tmp_coords[curr_dim] = last_coord;
	res[rank] = (*function)(tmp_coords);
}

result One_Dimensional_AGP::solve_mpi() {
	std::vector<double> data(range + 1);
	std::vector<int> k(range);
	std::map<int, std::pair<double, double>> processing_intervals_map;
	std::set<std::pair<double, double>> processing_intervals_set;
	std::pair<double, double> tmp_interval;
	result tmp_res;
	MPI_Status st;
	if (curr_dim == 0) {
		for (int i = 0; i < range; ++i)
			curr_x[i] = bounds[i].first;
	}
	std::pair<double, double> new_point;
	double new_m;

	perform_first_iteration();
	if (curr_dim == 0) {
		while (!isEnd() && points->size() < procnum) {
			new_m = get_m();
			compute_R(curr_x, new_m);
			new_point.first = get_new_point(pq->top()); 
			pq->pop();

			curr_x[curr_dim] = new_point.first;
			if (curr_dim != range - 1) {
				odm[curr_dim + 1]->set_experiment(range, curr_dim + 1, odm, bounds, curr_x,
					useMPI, useThreads, threadsNum, 
					function, eps, Nmax, r_p);
				odm[curr_dim + 1]->solve_mpi();
				tmp_res = odm[curr_dim + 1]->get_result();

				for (int i = curr_dim + 1; i < range; ++i)
					res.k[i] += tmp_res.k[i];
			}
			else {
				tmp_res.coords = curr_x;
				tmp_res.z = (*function)(curr_x);
			}

			curr_x = tmp_res.coords;
			insert_to_map(tmp_res.coords, tmp_res.z, 0);
			compare_interval_len(curr_x);
			compare_M(curr_x);
		}

		perform_first_parallel_step(new_point, tmp_res, new_m, processing_intervals_map, 
			processing_intervals_set, st);

		while (!isEnd()) {
			if (procrank == 0) {
				new_m = get_m();
				compute_R(curr_x, new_m);

				do {
					auto interv = pq->top();
					tmp_interval.first = interv.first_point.first[0];
					tmp_interval.second = interv.second_point.first[0];
					new_point.first = get_new_point(interv); 
					pq->pop();
				} while (processing_intervals_set.count(tmp_interval) > 0);

				curr_x[curr_dim] = new_point.first;
				processing_intervals_map.insert({ st.MPI_SOURCE, tmp_interval });
				processing_intervals_set.insert(tmp_interval);
				MPI_Send(&new_point.first, 1, MPI_DOUBLE, st.MPI_SOURCE, st.MPI_SOURCE, MPI_COMM_WORLD);

				MPI_Recv(data.data(), data.size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
				processing_intervals_set.erase(processing_intervals_map[st.MPI_SOURCE]);
				processing_intervals_map.erase(st.MPI_SOURCE);

				std::vector<double>::const_iterator first = data.begin();
				std::vector<double>::const_iterator last = data.begin() + range;

				tmp_res.z = data.back();
				tmp_res.coords = std::vector<double>(first, last);
				curr_x = tmp_res.coords;

				insert_to_map(tmp_res.coords, tmp_res.z, 0);
				compare_interval_len(curr_x);
				compare_M(curr_x);

				MPI_Recv(k.data(), k.size(), MPI_INT, st.MPI_SOURCE, st.MPI_SOURCE, MPI_COMM_WORLD, &st);
				for (int i = curr_dim + 1; i < range; ++i)
					res.k[i] += k[i];
			} else {
				MPI_Recv(&new_point.first, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

				if (st.MPI_TAG == 0)
					break;
					
				curr_x[curr_dim] = new_point.first;

				odm[curr_dim + 1]->set_experiment(range, curr_dim + 1, odm, bounds, curr_x, 
					useMPI, useThreads, threadsNum, 
					function, eps, Nmax, r_p);
				odm[curr_dim + 1]->solve_mpi();
				tmp_res = odm[curr_dim + 1]->get_result();

				for (int i = curr_dim + 1; i < range; ++i)
					res.k[i] += tmp_res.k[i];
				

				std::vector<double> data(tmp_res.coords);
				data.push_back(tmp_res.z);
				MPI_Send(data.data(), data.size(), MPI_DOUBLE, 0, procrank, MPI_COMM_WORLD);
				MPI_Send(tmp_res.k.data(), tmp_res.k.size(), MPI_INT, 0, procrank, MPI_COMM_WORLD);
			}
		}

		if (procrank == 0) {
			MPI_Send(&new_point.first, 1, MPI_DOUBLE, st.MPI_SOURCE, 0, MPI_COMM_WORLD);
			while (processing_intervals_set.size() > 0) {
				MPI_Recv(data.data(), data.size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
				processing_intervals_set.erase(processing_intervals_map[st.MPI_SOURCE]);
				processing_intervals_map.erase(st.MPI_SOURCE);
				MPI_Recv(k.data(), k.size(), MPI_INT, st.MPI_SOURCE, st.MPI_SOURCE, MPI_COMM_WORLD, &st);

				MPI_Send(&new_point.first, 1, MPI_DOUBLE, st.MPI_SOURCE, 0, MPI_COMM_WORLD);
			}
		}
	} else if (curr_dim == range - 1 
		&& procrank == 0 
		&& points->size() >= threadsNum - 1
		&& useThreads) {
		while (!isEnd()) {
			std::thread *pth = new std::thread[threadsNum];
			std::vector<double> last_coord(threadsNum);
			std::vector<double> result(threadsNum);

			new_m = get_m();
			compute_R(curr_x, new_m);
			for (int i = 0; i < threadsNum; ++i) {
				last_coord[i] = get_new_point(pq->top());
				pq->pop();
			}
				

			for (int i = 1; i < threadsNum; ++i)
				pth[i] = std::thread(&One_Dimensional_AGP::do_parallel_job, this,
					last_coord[i], std::ref(result), i);
			do_parallel_job(last_coord[0], result, 0);

			curr_x[curr_dim] = last_coord[0];
			tmp_res.coords = curr_x;
			tmp_res.z = result[0];

			insert_to_map(tmp_res.coords, tmp_res.z, 0);
			compare_interval_len(curr_x);
			compare_M(curr_x);

			for (int i = 1; i < threadsNum; ++i) {
				pth[i].join();

				curr_x[curr_dim] = last_coord[i];
				tmp_res.coords = curr_x;
				tmp_res.z = result[i];

				insert_to_map(tmp_res.coords, tmp_res.z, 0);
				compare_interval_len(curr_x);
				compare_M(curr_x);
			}
		}
	} else {
		while (!isEnd()) {
			new_m = get_m();
			compute_R(curr_x, new_m);
			new_point.first = get_new_point(pq->top()); 
			pq->pop();

			curr_x[curr_dim] = new_point.first;
			if (curr_dim != range - 1) {
				odm[curr_dim + 1]->set_experiment(range, curr_dim + 1, odm, bounds, curr_x, 
					useMPI, useThreads, threadsNum, 
					function, eps, Nmax, r_p);
				odm[curr_dim + 1]->solve_mpi();
				tmp_res = odm[curr_dim + 1]->get_result();

				for (int i = curr_dim + 1; i < range; ++i)
					res.k[i] += tmp_res.k[i];
			}
			else {
				tmp_res.coords = curr_x;
				tmp_res.z = (*function)(curr_x);
			}

			curr_x = tmp_res.coords;
			insert_to_map(tmp_res.coords, tmp_res.z, 0);
			compare_interval_len(curr_x);
			compare_M(curr_x);
		}
	}
	
	res.k[curr_dim] += points->size();
	delete_containers();
		
	return res;
}

result One_Dimensional_AGP::solve_seq() {
	result tmp_res;
	if (curr_dim == 0) {
		for (int i = 0; i < range; ++i)
			curr_x[i] = bounds[i].first;
	}
	std::pair<double, double> new_point;
	double new_m;

	std::thread *pth;
	if(curr_dim == range - 1
	   && useThreads)
		pth = new std::thread[threadsNum];

	perform_first_iteration();
	while (!isEnd()) {
		if (curr_dim == range - 1
			&& points->size() >= threadsNum + 1
			&& useThreads) {
				std::vector<double> last_coord(threadsNum);
				std::vector<double> result(threadsNum);

				new_m = get_m();
				compute_R(curr_x, new_m);
				for (int i = 0; i < threadsNum; ++i) {
					last_coord[i] = get_new_point(pq->top());
					pq->pop();
				}
					

				for (int i = 1; i < threadsNum; ++i)
					pth[i] = std::thread(&One_Dimensional_AGP::do_parallel_job, this,
						last_coord[i], std::ref(result), i);
				do_parallel_job(last_coord[0], result, 0);

				curr_x[curr_dim] = last_coord[0];
				insert_to_map(curr_x, result[0], 0);
				compare_interval_len(curr_x);
				compare_M(curr_x);

				for (int i = 1; i < threadsNum; ++i) {
					pth[i].join();

					curr_x[curr_dim] = last_coord[i];
					insert_to_map(curr_x, result[i], 0);
					compare_interval_len(curr_x);
					compare_M(curr_x);
				}
		} else {
			new_m = get_m();
			compute_R(curr_x, new_m);
			new_point.first = get_new_point(pq->top()); 
			pq->pop();
			curr_x[curr_dim] = new_point.first;
			if (curr_dim != range - 1) {
				odm[curr_dim + 1]->set_experiment(range, curr_dim + 1, odm, bounds, curr_x,
					useMPI, useThreads, threadsNum,
					function, eps, Nmax, r_p);
				odm[curr_dim + 1]->solve_seq();
				tmp_res = odm[curr_dim + 1]->get_result();

				for (int i = curr_dim + 1; i < range; ++i)
					res.k[i] += tmp_res.k[i];
			} else {
				tmp_res.coords = curr_x;
				tmp_res.z = (*function)(curr_x);
			}

			curr_x = tmp_res.coords;
			insert_to_map(tmp_res.coords, tmp_res.z, 0);
			compare_interval_len(curr_x);
			compare_M(curr_x);
		}
	}
	res.k[curr_dim] += points->size();
	delete_containers();

	return res;
}