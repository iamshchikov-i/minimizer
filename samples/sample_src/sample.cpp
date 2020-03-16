#include <iostream>

#include "Grishagin/grishagin_function.hpp"
#include "Grishagin/GrishaginProblemFamily.hpp"
#include "Grishagin/GrishaginConstrainedProblem.hpp"
#include "Grishagin/GrishaginConstrainedProblemFamily.hpp"

#include "two_dimensional_minimizer.h"

using std::cout;

int i;
TGrishaginProblemFamily grishFam;

double f(double x, double y) {
	return grishFam[i]->ComputeFunction({ x, y });
}

int main()
{
	double eps_par = 0.001, r = 1.1;
	vector<double> lb(2), ub(2);
	vector<double> delta1(2), delta2(2);
	const double eps = 0.01;
	double delta = 0.1;
	double delta_x, delta_y, actual_x, actual_y;

	double(*fptr)(double, double) = f;
	result res1;
	result res2;

	One_Dimensional_AGMND odm_agmnd(0, 0, 0, nullptr, eps_par, r), *p_odm_agmnd = &odm_agmnd;
	One_Dimensional_AGP odm_agp(0, 0, 0, nullptr, eps_par, r), *p_odm_agp = &odm_agp;
	Two_Dimensional_Minimizer m1(p_odm_agp, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r);
	Two_Dimensional_Minimizer m2(p_odm_agmnd, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r);

	for (int j = 28; j < 29; j++) {
		r = 1.1;
		i = j;
		actual_x = grishFam[j]->GetOptimumPoint()[0], actual_y = grishFam[j]->GetOptimumPoint()[1];
		grishFam[j]->GetBounds(lb, ub);

		do {
			m1.set_experiment(p_odm_agp, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r);
			m1.solve();
			res1 = m1.get_result();
			delta1[0] = std::abs(actual_x - res1.x);
			delta1[1] = std::abs(actual_y - res1.y);
			r += 0.1;
		}

		while (delta1[0] > eps || delta1[1] > eps);
	
		/*m2.set_experiment(p_odm_agmnd, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r);
		m2.solve();
		res2 = m2.get_result();
		delta2[0] = std::abs(actual_x - res2.x);
		delta2[1] = std::abs(actual_y - res2.y);*/
	}

	return 0;
}
