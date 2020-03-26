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
	double eps = 0.01, r1 = 1.1, r2 = 1.1;
	double delta1, delta2, delta = 0.1, k1, k2, res1, res2;
	Minimizer_v2 m1(intervals[0][0], intervals[0][1], pfn[0]);
	Minimizer m2(intervals[0][0], intervals[0][1], pfn[0]);
	
	
	for (int j = 0; j < n; j++) {
		r1 = 1.1, r2 = 1.1;
		do {
			flag = false;
			m1.set_experiment(intervals[j][0], intervals[j][1], pfn[j], 0.001, r1);
			m1.solve();
			res1 = m1.get_result().x;
			k1 = m1.get_result().k;
			for (int k = 0; k < res[j].size();k++) {
				delta1 = fabs(res1 - res[j][k]);
				if (delta1 <= eps) {
					flag = true;
					break;
				}	
			}
			r1 += delta;
		} while (!flag);
		
		
		do {
			flag = false;
			m2.set_experiment(intervals[j][0], intervals[j][1], pfn[j], 0.001, r2);
			m2.solve();
			res2 = m2.get_result().x;
			k2 = m2.get_result().k;
			for (int k = 0; k < res[j].size();k++) {
				delta2 = fabs(res2 - res[j][k]);
				if (delta2 <= eps) {
					flag = true;
					break;
				}
			}
			r2 += delta;
			
		} while (!flag);

		std::cout<<"| " << j << " | " << r1 << " | " << k1 << " | " << r2 << " | " << k2 << " | " << (double)k1 / (double)k2 << " |\n";
	}
	


	return 0;
}

