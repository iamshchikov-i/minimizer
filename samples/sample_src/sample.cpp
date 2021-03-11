#include <iostream>

#include "Grishagin/grishagin_function.hpp"
#include "Grishagin/GrishaginProblemFamily.hpp"
#include "Grishagin/GrishaginConstrainedProblem.hpp"
#include "Grishagin/GrishaginConstrainedProblemFamily.hpp"

#include "GKLS/GKLSProblem.hpp"
#include "GKLS/GKLSProblemFamily.hpp"
#include "GKLS/GKLSConstrainedProblem.hpp"
#include "GKLS/GKLSConstrainedProblemFamily.hpp"

//#include "one_dimensional_agp_discontinuous.h"
//#include "two_dimensional_minimizer.h"
#include "multi_dimensional_minimizer.h"
using std::cout;

#include <chrono>

int i;
int range = 4;
TGrishaginProblemFamily grishFam;
TGKLSProblemFamily gklsFam(range);

double f_grish(std::vector<double> coords) {
	return grishFam[i]->ComputeFunction({ coords[0], coords[1] });
}

double f_gkls(std::vector<double> coords) {
	return gklsFam[i]->ComputeFunction(coords);
}

void get_odm_values();

int main()
{
	std::vector<int> task_nums = { 0, 1, 2, 3, 4 };
	std::string status_agp, status_agmnd;
	result res_agp, res_agmnd;
	double eps = 0.01;
	std::vector<double> lower_bound;
	std::vector<double> upper_bound;

	vector<double> lb, ub;

	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	std::chrono::milliseconds elapsed_ms_agp, elapsed_ms_agmnd;

	for (int j : task_nums) {
		i = j;
		gklsFam[i]->GetBounds(lb, ub);
		for (int k = 0; k < range; ++k) {
			lower_bound.push_back(lb[k]);
			upper_bound.push_back(ub[k]);
		}

		Multi_Dimensional_Minimizer mdm_agp(range, lower_bound, upper_bound, f_gkls, Upper_method::AGP);

		begin = std::chrono::steady_clock::now();
		res_agp = mdm_agp.solve();
		end = std::chrono::steady_clock::now();

		elapsed_ms_agp =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

		Multi_Dimensional_Minimizer mdm_agmnd(range, lower_bound, upper_bound, f_gkls, Upper_method::AGMND);

		begin = std::chrono::steady_clock::now();
		res_agmnd = mdm_agmnd.solve();
		end = std::chrono::steady_clock::now();

		elapsed_ms_agmnd =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

		double actual_res = gklsFam[i]->GetOptimumValue();
		if (abs(abs(res_agp.z) - abs(actual_res)) <= eps)
			status_agp = "OK";
		else
			status_agp = "Fail";

		if (abs(abs(res_agmnd.z) - abs(actual_res)) <= eps)
			status_agmnd = "OK";
		else
			status_agmnd = "Fail";

		std::cout << "AGP " << j << " " << status_agp << ", eps =  " << abs(abs(res_agp.z) - abs(actual_res)) << ", total time = " << elapsed_ms_agp.count() << std::endl;
		//print(res_agp.k);

		std::cout << "AMND " << j << " " << status_agmnd << ", eps =  " << abs(abs(res_agmnd.z) - abs(actual_res)) << ", total time = " << elapsed_ms_agmnd.count() << std::endl;
		//print(res_agmnd.k);
		//std::cout << std::endl;
	}

	return 0;
}

void get_odm_values() {
	double h = 0.001;
	double z, min_z, tmp, min_x, min_y;
	vector<double> lb(2), ub(2);
	
	for (int j = 1; j < 2; j++) {
		i = j;
		gklsFam[i]->GetBounds(lb, ub);
		min_x = lb[0];
		min_y = lb[1];
		
		for (double x = lb[0]; x < ub[0]; x += h) {
			min_z = f_gkls({ x, lb[1] });
			for (double y = lb[1]; y < ub[1]; y += h) {
				tmp = f_gkls({ x, y });
				if (tmp < min_z) {
					min_z = tmp;
					min_x = x;
					min_y = y;
				}
			}

			std::cout << x << ":" << min_z << std::endl;
		}
		//std::cout << min_x << " " << min_y << " : " << min_z << std::endl;
	}
}