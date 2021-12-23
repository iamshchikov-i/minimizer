#include <chrono>
#include "minimizer_v1.h"
#include "minimizer_v2.h"
#include "hansen_functions.h"

//void measure_time_v1_v2() {
//	std::vector<std::pair<double, double> > bord(6);
//	bord[0] = { -15.0, 8.0 };  bord[1] = { -10.0, 5.0 }; bord[2] = { 0.12, 10.0 };
//	bord[3] = { -10.2, 10.2 }; bord[4] = { 0.0, 0.1 };   bord[5] = { 0.1, 10.57 };
//	std::chrono::milliseconds r1{ 0 }, r2{ 0 }, elapsed_ms;
//	std::chrono::time_point<std::chrono::steady_clock> begin, end;
//	double(*fptr[6])(double) = { f1,f2,f3,f4,f5,f6 };
//
//	std::cout << "Via Minimizer_v1: " << std::endl;
//	for (int i = 0;i < 6;++i) {
//		Minimizer_v1 m(bord[i].first, bord[i].second, fptr[i]);
//		begin = std::chrono::steady_clock::now();
//		std::cout << "x* of f" << i << " = " << m.find_point() << std::endl;
//		end = std::chrono::steady_clock::now();
//		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
//		r1 += elapsed_ms;
//	}
//	std::cout << "Total time 1: " << r1.count() << " ms\n";
//
//	std::cout << std::endl << "Via Minimizer_v2: " << std::endl;
//	Minimizer_v2 m(bord[0].first, bord[0].second, fptr[0]);
//	for (int i = 0;i < 6;++i) {
//		m.set_experiment(bord[i].first, bord[i].second, fptr[i]);
//		begin = std::chrono::steady_clock::now();
//		std::cout << "x* of f" << i << " = " << m.solve().x << std::endl;
//		end = std::chrono::steady_clock::now();
//		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
//		r2 += elapsed_ms;
//	}
//	std::cout << "Total time 2: " << r2.count() << " ms\n";
//}

int main(int argc, char** argv) {
	int n = 20;
	int k;

	std::vector<int> taskList = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };

	std::chrono::milliseconds r1{ 0 }, r2{ 0 }, elapsed_ms;
	std::chrono::time_point<std::chrono::steady_clock> begin, end;

	Minimizer_v1 agp(intervals[0][0], intervals[0][1], pfn[0]);

	for (auto j : taskList) {
		Minimizer_v1 agp(intervals[j][0], intervals[j][1], pfn[j], 0.0001, 1000, 5.0);
		begin = std::chrono::steady_clock::now();
		agp.find_point();
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		std::cout << elapsed_ms.count() << std::endl;
	}

	/*Minimizer_v2 agp(intervals[0][0], intervals[0][1], pfn[0]);

	for (auto j : taskList) {
		Minimizer_v2 agp(intervals[j][0], intervals[j][1], pfn[j], 0.0001, 1000, 5.0);
		begin = std::chrono::steady_clock::now();
		agp.solve();
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		std::cout << elapsed_ms.count() << std::endl;
	}*/

	return 0;
}
