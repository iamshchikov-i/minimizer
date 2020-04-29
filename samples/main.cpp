#include <thread>
#include <chrono>
//#include "minimizer_v1.h"
#include "minimizer_v2.h"
//#include "minimizer_v3.h"
#include "functions.h"
#include "AGMND.h"
#include "hansen_functions.h"

void set_bords(std::vector<double>&);

void set_actual_values(std::vector<double>&);

void time();

void measure_time_of_function(double(*f)(double x));

void measure_time_v2_v3();

int main(int argc, char** argv) {
	bool flag = false;
	int n = 20;
	double eps = 0.01, r1 = 1.1;
	_result res;

	Minimizer agmnd(intervals[0][0], intervals[0][1], pfn[0]);
	for (int j = 0; j < 20; j++) {

		agmnd.set_experiment(intervals[j][0], intervals[j][1], pfn[j], 0.0001, r1);
		agmnd.solve();
		res = agmnd.get_result();
		
		std::cout << j <<"/"<< res.k <<"\n";
	} 
		

	return 0;
}

