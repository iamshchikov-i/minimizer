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
#include <assert.h>

struct ResultInfo {
	int successCount;
	int failCount;
	double averageTime;
	double totalTime;

ResultInfo(int _successCount, int _failCount, double _averageTime, double _totalTime);
};
ResultInfo::ResultInfo(int _successCount, int _failCount, double _averageTime, double _totalTime) :
	successCount(_successCount), failCount(_failCount), 
	averageTime(_averageTime), totalTime(_totalTime) {}

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
bool check_result_coords(std::vector<double>& expected,
						 std::vector<double>& actual, double eps);

int main()
{
	//std::vector<int> task_nums = { 0, 1, 2, 3, 4 };
	std::string status_agp, status_agmnd;
	result res_agp, res_agmnd;
	double eps = 0.01;
	double eps_par = 0.001, r_par = 2.5;
	int Nmax = 100000000;
	std::vector<double> lower_bound;
	std::vector<double> upper_bound;

	vector<double> lb, ub;

	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	std::chrono::milliseconds elapsed_ms_agp, elapsed_ms_agmnd;

	int taskNumber = gklsFam.GetFamilySize();
	ResultInfo resInfoAgp(0, 0, 0.0, 0.0), resInfoAgmnd(0, 0, 0.0, 0.0);

	for (int j = 0; j < taskNumber; ++j) {
		i = j;
		gklsFam[i]->GetBounds(lb, ub);
		for (int k = 0; k < range; ++k) {
			lower_bound.push_back(lb[k]);
			upper_bound.push_back(ub[k]);
		}

		Multi_Dimensional_Minimizer mdm_agp(range, lower_bound, upper_bound, f_gkls, Upper_method::AGP, eps_par, Nmax, r_par);

		begin = std::chrono::steady_clock::now();
		res_agp = mdm_agp.solve();
		end = std::chrono::steady_clock::now();

		elapsed_ms_agp =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

		vector<double> actual_res = gklsFam[i]->GetOptimumPoint();
		bool res = check_result_coords(res_agp.coords, actual_res, eps);
		if (res) {
			status_agp = "OK";
			resInfoAgp.successCount++;
		} else {
			status_agp = "Fail";
			resInfoAgp.failCount++;
		}
		resInfoAgp.totalTime += elapsed_ms_agp.count();
			
		std::cout << "AGP " << j << " " << status_agp << ", total time = " << elapsed_ms_agp.count() << std::endl;
		std::cout << "number of processed points: ";
		print(res_agp.k);
		std::cout << std::endl;

		Multi_Dimensional_Minimizer mdm_agmnd(range, lower_bound, upper_bound, f_gkls, Upper_method::AGMND, eps_par, Nmax, r_par);

		begin = std::chrono::steady_clock::now();
		res_agmnd = mdm_agmnd.solve();
		end = std::chrono::steady_clock::now();

		elapsed_ms_agmnd =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

		res = check_result_coords(res_agmnd.coords, actual_res, eps);
		if (res) {
			status_agmnd = "OK";
			resInfoAgmnd.successCount++;
		}
		else {
			status_agmnd = "Fail";
			resInfoAgmnd.failCount++;
		}
		resInfoAgmnd.totalTime += elapsed_ms_agmnd.count();

		std::cout << "AGMND " << j << " " << status_agmnd << ", total time = " << elapsed_ms_agmnd.count() << std::endl;
		std::cout << "number of processed points: ";
		print(res_agmnd.k);
		std::cout << std::endl << std::endl;
	}
	resInfoAgp.averageTime = resInfoAgp.totalTime / taskNumber;
	resInfoAgmnd.averageTime = resInfoAgmnd.totalTime / taskNumber;

	std::cout << "AGP\n" << "success-count fail-count average-time\n";
	std::cout << resInfoAgp.successCount << " " << resInfoAgp.failCount << " " << resInfoAgp.averageTime << std::endl;

	std::cout << "AGMND\n" << "success-count fail-count average-time\n";
	std::cout << resInfoAgmnd.successCount << " " << resInfoAgmnd.failCount << " " << resInfoAgmnd.averageTime << std::endl;

	return 0;
}

bool check_result_coords(std::vector<double>& expected,
						 std::vector<double>& actual, double eps) {
	bool status = true;

	assert(expected.size() == actual.size());
	for (int i = 0; i < actual.size(); ++i) {
		if (abs(expected[i] - actual[i]) > eps) {
			std::cout << "coords " << i << " are too much different " << expected[i] << " " << actual[i] << std::endl;
			status = false;
		}
	}

	return status;
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