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

#include <chrono>
#include <assert.h>
#include <map>

using std::cout;

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
int procrank;
int range = 2;
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
void execExperiment(int task_num, double eps_par, double r_par,
	int Nmax, Upper_method um, double eps = 0.01);

ResultInfo resInfoAgp(0, 0, 0.0, 0.0), resInfoAgmnd(0, 0, 0.0, 0.0);
std::map<std::string, ResultInfo*> resultsInfo;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	//std::vector<int> task_nums = { 0, 1, 2, 3, 4 };
	double eps = 0.01;
	double eps_par = 0.001, r_par = 2.5;
	int Nmax = 100000000;

	int taskNumber = gklsFam.GetFamilySize();
	resultsInfo.insert({ "AGP", &resInfoAgp });
	resultsInfo.insert({ "AGMND", &resInfoAgmnd });

	for (int j = 0; j < taskNumber; ++j) {
		execExperiment(j, eps_par, r_par, Nmax, Upper_method::AGP, eps);
		execExperiment(j, eps_par, r_par, Nmax, Upper_method::AGMND, eps);
	}
	
	if (procrank == 0) {
		resInfoAgp.averageTime = resInfoAgp.totalTime / taskNumber;
		resInfoAgmnd.averageTime = resInfoAgmnd.totalTime / taskNumber;

		std::cout << "AGP\n" << "success-count fail-count average-time\n";
		std::cout << resInfoAgp.successCount << " " << resInfoAgp.failCount << " " << resInfoAgp.averageTime << std::endl;

		std::cout << "AGMND\n" << "success-count fail-count average-time\n";
		std::cout << resInfoAgmnd.successCount << " " << resInfoAgmnd.failCount << " " << resInfoAgmnd.averageTime << std::endl;

	}

	MPI_Finalize();

	return 0;
}

void execExperiment(int task_num, double eps_par, double r_par,
	int Nmax, Upper_method um, double eps) {
	std::string status;
	result result;
	std::vector<double> lower_bound;
	std::vector<double> upper_bound;
	std::vector<double> actual_res;
	bool res;

	vector<double> lb, ub;
	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	std::chrono::milliseconds elapsed_ms;

	int j = task_num;
	i = j;
	gklsFam[i]->GetBounds(lb, ub);
	for (int k = 0; k < range; ++k) {
		lower_bound.push_back(lb[k]);
		upper_bound.push_back(ub[k]);
	}

	int procnum;

	MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	MPI_Barrier(MPI_COMM_WORLD);


	Multi_Dimensional_Minimizer mdm(range, lower_bound, upper_bound, f_gkls, um, eps_par, Nmax, r_par);

	if (procrank == 0)
		begin = std::chrono::steady_clock::now();
	result = mdm.solve();
	MPI_Barrier(MPI_COMM_WORLD);
	if (procrank == 0) {
		end = std::chrono::steady_clock::now();

		elapsed_ms =
			std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	}

	std::string type;
	if (um == Upper_method::AGP)
		type = "AGP";
	else if (um == Upper_method::AGMND)
		type = "AGMND";
	else throw - 1;

	if (procrank == 0) {
		actual_res = gklsFam[i]->GetOptimumPoint();
		res = check_result_coords(result.coords, actual_res, eps);
		if (res) {
			status = "OK";
			resultsInfo[type]->successCount++;
		} else {
			status = "Fail";
			resultsInfo[type]->failCount++;
		}
		resultsInfo[type]->totalTime += elapsed_ms.count();
			
		std::cout << type << j << " " << status << ", total time = " << elapsed_ms.count() << std::endl;
		std::cout << "number of processed points: ";
		print(result.k);
		std::cout << std::endl;
	}
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