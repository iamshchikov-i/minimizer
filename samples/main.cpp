#include <chrono>
#include "minimizer_v1.h"
#include "minimizer_v2.h"
#include "functions.h"

int main() {
	std::vector<std::pair<double, double> > bord(6);
	bord[0] = { -15.0, 8.0 };  bord[1] = { -10.0, 5.0 }; bord[2] = { 0.12, 10.0 };
	bord[3] = { -10.2, 10.2 }; bord[4] = { 0.0, 0.1 };   bord[5] = { 0.1, 10.57 };
	std::chrono::milliseconds r1{ 0 }, r2{ 0 }, elapsed_ms;
	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	double(*fptr[6])(double) = { f1,f2,f3,f4,f5,f6 };
	
	std::cout << "Via Minimizer_v1: " << std::endl;
	for (int i = 0;i < 6;++i) {
		Minimizer_v1 m(bord[i].first, bord[i].second, fptr[i]);
		begin = std::chrono::steady_clock::now();
		std::cout<<"x* of f"<<i<< " = "<< m.find_point() << std::endl;
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		r1 += elapsed_ms;
	}
	std::cout << "Total time 1: " << r1.count() << " ms\n";

	std::cout << std::endl << "Via Minimizer_v2: " << std::endl;
	Minimizer_v2 m(bord[0].first, bord[0].second, fptr[0]);
	for (int i = 0;i < 6;++i) {
		m.set_experiment(bord[i].first, bord[i].second, fptr[i]);
		begin = std::chrono::steady_clock::now();
		std::cout << "x* of f" << i << " = " << m.solve().x << std::endl;
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		r2 += elapsed_ms;
	}
	std::cout << "Total time 2: " << r2.count() << " ms\n";

	return 0;
}
