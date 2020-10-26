#include <iostream>

#include "Grishagin/grishagin_function.hpp"
#include "Grishagin/GrishaginProblemFamily.hpp"
#include "Grishagin/GrishaginConstrainedProblem.hpp"
#include "Grishagin/GrishaginConstrainedProblemFamily.hpp"

#include "one_dimensional_agp_discontinuous.h"
#include "two_dimensional_minimizer.h"

using std::cout;

int i;
TGrishaginProblemFamily grishFam;

double f(double x, double y) {
	return grishFam[i]->ComputeFunction({ x, y });
}

void grish_fam();
void get_odm_values();

int main()
{
	grish_fam();

	return 0;
}


void grish_fam() {
	double eps_par = 0.001, r = 2.0;
	vector<double> lb(2), ub(2);
	const double eps = 0.01;
	double delta_x, delta_y, actual_x, actual_y;

	double(*fptr)(double, double) = f;
	result res;

	One_Dimensional_AGP_D odm_agp_d(0, 0, 0, nullptr, eps_par, r), *p_odm_agp_d = &odm_agp_d;
	Two_Dimensional_Minimizer m(p_odm_agp_d, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r);


	for (int j = 0; j < grishFam.GetFamilySize(); j++) {
		r = 2.0;
		i = j;
		actual_x = grishFam[j]->GetOptimumPoint()[0], actual_y = grishFam[j]->GetOptimumPoint()[1];
		grishFam[j]->GetBounds(lb, ub);

		p_odm_agp_d->set_experiment_d(0, 0, 0, nullptr);
		m.set_experiment(p_odm_agp_d, lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r);
		m.solve();
		res = m.get_result();
		//std::cout << res.k << std::endl;
		//std::cout <<res.k_on_x<<"/"<<res.k_max_on_y << std::endl;

	}

}

void get_odm_values() {
	double h = 0.01;
	double z, min_z, tmp;
	vector<double> lb(2), ub(2);
	
	for (int j = 0; j < grishFam.GetFamilySize(); j++) {
		i = j;
		grishFam[i]->GetBounds(lb, ub);

		for (double x = lb[0]; x < ub[0]; x += h) {
			min_z = f(x, lb[1]);
			for (double y = lb[1]; y < ub[1]; y += h) {
				tmp = f(x, y);
				if (tmp < min_z) {
					min_z = tmp;
				}
			}

			std::cout << x << ":" << min_z << std::endl;
		}
		std::cout << std::endl;
	}


}