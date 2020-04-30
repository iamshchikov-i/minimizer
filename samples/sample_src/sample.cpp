#include <iostream>

#include "Grishagin/grishagin_function.hpp"
#include "Grishagin/GrishaginProblemFamily.hpp"
#include "Grishagin/GrishaginConstrainedProblem.hpp"
#include "Grishagin/GrishaginConstrainedProblemFamily.hpp"

#include "two_dimensional_minimizer.h"
#include "functions.h"

using std::cout;

int i;
TGrishaginProblemFamily grishFam;

double f(double x, double y) {
	return grishFam[i]->ComputeFunction({ x, y });
}

void grish_fam();

int main()
{
	grish_fam();

	return 0;
}


void grish_fam() {
	double eps_par = 0.001, r1 = 1.1, r2 = 1.1;
	vector<double> lb(2), ub(2);
	vector<double> delta1(2), delta2(2);
	const double eps = 0.01;
	double delta = 0.1;
	double delta_x, delta_y, actual_x, actual_y;

	double(*fptr)(double, double) = f;
	result res1;
	result res2;

	One_Dimensional_AGMND odm_agmnd(0, 0, 0, nullptr, eps_par / 10, r1), *p_odm_agmnd = &odm_agmnd;
	One_Dimensional_AGP odm_agp(0, 0, 0, nullptr, eps_par, r1), *p_odm_agp = &odm_agp;
	Two_Dimensional_Minimizer m1(p_odm_agp, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r2);
	Two_Dimensional_Minimizer m2(p_odm_agmnd, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r2);

	for (int j = 0; j < grishFam.GetFamilySize(); j++) {
		r1 = 2.0, r2 = 2.0;
		i = j;
		actual_x = grishFam[j]->GetOptimumPoint()[0], actual_y = grishFam[j]->GetOptimumPoint()[1];
		grishFam[j]->GetBounds(lb, ub);

		do {
			p_odm_agp->set_experiment(0, 0, 0, nullptr, r1);
			m1.set_experiment(p_odm_agp, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r2);
			m1.solve();
			res1 = m1.get_result();
			delta1[0] = std::abs(actual_x - res1.x);
			delta1[1] = std::abs(actual_y - res1.y);
			if (r1 < 3.6) r1 += delta;
			else {
				r1 = 1.1;
				r2 += delta;
			}
		} while (delta1[0] > eps || delta1[1] > eps);
		printf("r1 = %f, r2 = %f, k_x = %d max_k_y = %d k = %d\n", r1, r2, res1.k_on_x, res1.k_max_on_y, res1.k);

		r1 = 1.1, r2 = 1.1;
		do {
			p_odm_agmnd->set_experiment(0, 0, 0, nullptr, r1);
			m2.set_experiment(p_odm_agmnd, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r2);
			m2.solve();
			res2 = m2.get_result();
			delta2[0] = std::abs(actual_x - res2.x);
			delta2[1] = std::abs(actual_y - res2.y);
			if (r1 < 3.6) r1 += delta;
			else {
				r1 = 1.1;
				r2 += delta;
			}
		} while (delta2[0] > eps || delta2[1] > eps);
		printf("r1 = %f, r2 = %f, k_x = %d max_k_y = %d k = %d\n\n", r1, r2, res2.k_on_x, res2.k_max_on_y, res2.k);

	}

}