#include <chrono>
#include "minimizer_v1.h"
#include "minimizer_v2.h"
#include "minimizer_v3.h"
#include "functions.h"

void measure_time_v1_v2() {
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
		std::cout << "x* of f" << i << " = " << m.find_point() << std::endl;
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
}

void measure_time_of_function(double(*f)(double x)) {
	double y;
	std::chrono::milliseconds elapsed_ms;
	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	begin = std::chrono::steady_clock::now();
	y = (*f)(5);
	end = std::chrono::steady_clock::now();
	elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << elapsed_ms.count() << " ms\n";
}

void measure_time_v2_v3() {
	int rank;
	std::pair<double, double> bord;
	std::chrono::milliseconds elapsed_ms;
	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	bord = { -15.0, 8.0 };
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		/*std::cout << "Via Minimizer_v1: " << std::endl;

		Minimizer_v1 m1(bord.first, bord.second, f1);
		begin = std::chrono::steady_clock::now();
		m1.find_point();
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		std::cout << "Total time 1: " << elapsed_ms.count() << " ms\n";*/

		std::cout << std::endl << "Via Minimizer_v2: " << std::endl;
		Minimizer_v2 m2(bord.first, bord.second, f1);
		begin = std::chrono::steady_clock::now();
		m2.solve();
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		std::cout << "Total time 2: " << elapsed_ms.count() << " ms\n";
	}

	Minimizer_v3 m3(bord.first, bord.second, f1);
	MPI_Barrier(MPI_COMM_WORLD);
	if (m3.get_rank() == 0) {
		std::cout << std::endl << "Via Minimizer_v3: " << std::endl;
		begin = std::chrono::steady_clock::now();
		m3.solve();
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		std::cout << "Total time 3: " << elapsed_ms.count() << " ms\n";
	}
	else
		m3.solve();
}

int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);

	measure_time_v2_v3();

	MPI_Finalize();
	return 0;
}
