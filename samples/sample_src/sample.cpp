#include <iostream>

#include "Grishagin/grishagin_function.hpp"
#include "Grishagin/GrishaginProblemFamily.hpp"
#include "Grishagin/GrishaginConstrainedProblem.hpp"
#include "Grishagin/GrishaginConstrainedProblemFamily.hpp"

#include"minimizer_seq.h"
#include "AGMND.h"
#include "functions.h"
#include "assert.h"

using std::cout;

int i;
TGrishaginProblemFamily grishFam;

double f(double x, double y) {
	return grishFam[i]->ComputeFunction({ x, y });
}

int main()
{
	double eps_par = 0.001, r = 3.3;
	vector<double> lb(2), ub(2);
	vector<double> delta1(2), delta2(2);
	const double eps = 0.01;
	double delta = 0.1;

	for (int j = 0; j < grishFam.GetFamilySize(); j++) {
		i = j;
		eps_par = 0.001;
		double(*fptr)(double, double) = f;
		result res1;
		_result res2;
		double delta_x, delta_y, actual_x = grishFam[j]->GetOptimumPoint()[0], actual_y = grishFam[j]->GetOptimumPoint()[1];
		grishFam[j]->GetBounds(lb, ub);

		Minimizer m(lb[0], ub[0], lb[1], ub[1], fptr, eps_par, 500, r);
		m.solve();
		res1 = m.get_result();
		delta1[0] = std::abs(actual_x - res1.x);
		delta1[1] = std::abs(actual_y - res1.y);

		AGMND agmnd(lb[0], ub[0], lb[1], ub[1], fptr, eps_par, r);
		agmnd.solve();
		res2 = agmnd.get_result();
		delta2[0] = std::abs(actual_x - res2.x);
		delta2[1] = std::abs(actual_y - res2.y);
		
		cout << "| " << j << " | " << res1.k << " | " << res2.k << " | " << (double)res1.k / (double)res2.k << " |" << std::endl;
	}

	return 0;
}
