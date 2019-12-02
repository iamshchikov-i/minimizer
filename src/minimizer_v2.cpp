#include "minimizer_v2.h"

characteristics::characteristics(double _y, double _R): y(_y), R(_R) {}

characteristics::characteristics(double _R): R(_R) {}

interval::interval(std::pair<double, double> _f_point, std::pair<double, double> _s_point, double _R) : first_point(_f_point), second_point(_s_point), _ch(_R) {}

Minimizer_v2::Minimizer_v2(double _a, double _b, double(*f)(double x), double _eps, int _N_max, double _r_par): a(_a), b(_b), function(f),
																										eps(_eps), N_max(_N_max), r_p(_r_par) {
	res.k = 0;
	values = new std::map<double, characteristics>;
	pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
}

bool Minimizer_v2::stop1() {
	return res.k >= N_max;
}

bool Minimizer_v2::stop2() {
	return min_interval_length <= eps;
}

bool Minimizer_v2::isEnd() {
	return stop1() || stop2();
}

void Minimizer_v2::go_Next_Interval() {
	left_point++; right_point++;
}

void Minimizer_v2::go_new_left_interval(double new_point) {
	right_point = values->find(new_point);
	left_point = right_point; left_point--;
}

void Minimizer_v2::reset() {
	left_point = values->begin();
	right_point = left_point; right_point++;
}

double Minimizer_v2::get_M() {
	return abs(((*right_point).second.y - (*left_point).second.y) / ((*right_point).first - (*left_point).first));
}

double Minimizer_v2::get_R() {
	double tmp = m * ((*right_point).first - (*left_point).first);
	return tmp + (pow((*right_point).second.y - (*left_point).second.y, 2) / tmp) - 2 * ((*right_point).second.y + (*left_point).second.y);
}

void Minimizer_v2::calculate_R(double new_point, double new_m) {
	if (new_m != m) {
		if (!pq->empty()) {
			delete pq;
			pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
		}
		m = new_m;
		reset();
		for (; right_point != values->end(); go_Next_Interval()) {
			(*left_point).second.R = get_R();
			pq->push(interval({ (*left_point).first, (*left_point).second.y },
				{ (*right_point).first, (*right_point).second.y }, (*left_point).second.R) );
		}
	}
	else {
		go_new_left_interval(new_point);
		for (int i = 0;i < 2;++i, go_Next_Interval()) {
			(*left_point).second.R = get_R();
			pq->push(interval({ (*left_point).first, (*left_point).second.y }, 
				{ (*right_point).first, (*right_point).second.y }, (*left_point).second.R));
		}
	}
}

double Minimizer_v2::get_m() {
	if (M_Max > 0)
		return r_p * M_Max;
	else if (M_Max == 0)
		return 1;
	else
		throw - 1;
}

void Minimizer_v2::insert_to_map(double _x, double _y, double _R) {
	characteristics _ch(_y, _R);
	values->insert({ _x, _ch });
	if (_y < res.y) {
		res.x = _x;
		res.y = _y;
	}
}

void Minimizer_v2::compare_interval_len(double new_point) {
	go_new_left_interval(new_point);
	double interval_length;
	for (int i = 0;i < 2;++i, go_Next_Interval()) {
		interval_length = (*right_point).first - (*left_point).first;
		if (interval_length < min_interval_length)
			min_interval_length = interval_length;
	}
}

void Minimizer_v2::compare_M(double new_point) {
	go_new_left_interval(new_point);
	double M;
	for (int i = 0;i < 2;++i, go_Next_Interval()) {
		M = get_M();
		if (M >= M_Max)
			M_Max = M;
	}
}

double Minimizer_v2::get_new_point(interval i) {
	return 0.5*(i.first_point.first + i.second_point.first) - ((i.second_point.second) - i.first_point.second) / (2 * m);
}

void Minimizer_v2::delete_containers() {
	delete values;
	values = nullptr;
	delete pq;
	pq = nullptr;
}

void Minimizer_v2::set_experiment(double _a, double _b, double(*f)(double x),
	double _eps, int _N_max, double _r_par) {
	a = _a;
	b = _b;
	function = f;
	eps = _eps;
	N_max = _N_max;
	r_p = _r_par;
	res.k = 0;
	if(values == nullptr)
		values = new std::map<double, characteristics>;
	if(pq == nullptr)
		pq = new std::priority_queue<interval, std::vector<interval>, CompareR>;
}

result Minimizer_v2::get_result() {
	return res;
}

result Minimizer_v2::solve() {
	std::chrono::milliseconds elapsed_ms, t[7];
	for (int i = 0; i < 7;i++)
		t[i] = { 0 };
	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	std::pair<double, double> new_point;
	double new_m;

	
	//1я итерация
	

	min_interval_length = b - a;
	res.x = a; res.y = (*function)(a);
	insert_to_map(res.x, res.y, 0);
	insert_to_map(b, (*function)(b), 0);
	reset();
	M_Max = get_M();
	m = -1;
	pq->push(interval({(*left_point).first, (*left_point).second.y }, { (*right_point).first, (*right_point).second.y }, 0));
	res.k = 2;

	
	while (!isEnd()) {

		begin = std::chrono::steady_clock::now();
		
		new_m = get_m();
		
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		t[0] += elapsed_ms;

		

		begin = std::chrono::steady_clock::now();
		calculate_R(new_point.first, new_m);
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		t[1] += elapsed_ms;

		begin = std::chrono::steady_clock::now();
		new_point.first = get_new_point(pq->top()); pq->pop();
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		t[2] += elapsed_ms;


		begin = std::chrono::steady_clock::now();
		new_point.second = (*function)(new_point.first);
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		t[3] += elapsed_ms;

		begin = std::chrono::steady_clock::now();
		insert_to_map(new_point.first, new_point.second, 0);
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		t[4] += elapsed_ms;

		res.k++; 


		begin = std::chrono::steady_clock::now();
		compare_interval_len(new_point.first);
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		t[5] += elapsed_ms;

		begin = std::chrono::steady_clock::now();
		compare_M(new_point.first);
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		t[6] += elapsed_ms;
	}

	std::cout << "Time of get m = " << t[0].count() << std::endl;
	std::cout << "Time of calculate R = " << t[1].count() << std::endl;
	std::cout << "Time of get new point = " << t[2].count() << std::endl;
	std::cout << "Time of get value of new point = " << t[3].count() << std::endl;
	std::cout << "Time of insert to map = " << t[4].count() << std::endl;
	std::cout << "Time of get compare interval len = " << t[5].count() << std::endl;
	std::cout << "Time of get compare M = " << t[6].count() << std::endl;
	
	return res;
}