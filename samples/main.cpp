#include <chrono>
#include "hansen_functions.h"
#include "one_dimensional_agp_discontinuous.h"

int main(int argc, char** argv) {
	int n = 20;
	int k;
	One_Dimensional_AGP_D agp_d(intervals[0][0], intervals[0][1], pfn[0]);

	for (int j = 0; j < n; j++) {
		agp_d.set_experiment(intervals[j][0], intervals[j][1], pfn[j]);
		agp_d.solve();
		k = agp_d.get_result().k;
		
		std::cout << k << std::endl;
	}

	return 0;
}

