#include <chrono>
#include "windows.h"
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

void set_bords(std::vector<double>& bord) {
	bord[0] = -15.0; bord[1] = 8.0; bord[2] = 6.36; bord[3] = 20.0; bord[4] = -7.56; bord[5] = -1.92;
	bord[6] = -2.0; bord[7] = 5.0; bord[8] = -0.21; bord[9] = 10.0; bord[10] = -7.56; bord[11] = -1.92;
	bord[12] = 0.12; bord[13] = 10.0; bord[14] = 4.23; bord[15] = 10.0; bord[16] = -5.78; bord[17] = -1.92;
	bord[18] = -10.2; bord[19] = 10.2; bord[20] = 0.12; bord[21] = 0.98; bord[22] = -19.2; bord[23] = -5.87;
	bord[24] = 0.0; bord[25] = 1.0; bord[26] = 0.923; bord[27] = 1.338; bord[28] = -1.806; bord[29] = -1.448;
	bord[30] = 0.1; bord[31] = 10.57; bord[32] = 0.1; bord[33] = 1.298; bord[34] = 10.53; bord[35] = 19.23;
}

void set_actual_values(std::vector<double>& values) {
	values[0] = 5.46; values[1] = 6.36; values[2] = 1.92;
	values[3] = 0.5; values[4] = 0.21; values[5] = 1.92;
	values[6] = 3.0; values[7] = 4.23; values[8] = 5.78;
	values[9] = 1.0; values[10] = 0.12; values[11] = 5.87;
	values[12] = 0.8028; values[13] = 0.923; values[14] = 1.448;
	values[15] = 1.353; values[16] = 1.298; values[17] = 10.53;
}

void test_parallel_version() {
	int test_num = 18;
	std::vector<double> bord(test_num * 2);
	std::vector<double> actual_value(test_num);
	double(*fptr[6])(double) = {f1, f2, f3, f4, f5, f6};
	const double eps = 0.01;
	double delta = 0.0;
	
	
	set_bords(bord);
	set_actual_values(actual_value);
	
	Minimizer_v3 m(0, 0, nullptr);
	
	for (int i = 0; i < test_num; i++) {
		m.set_experiment(bord[2 * i], bord[(2 * i) + 1], fptr[i / 3]);
		
		m.solve();
		
		delta = abs(abs(m.get_result().x) - actual_value[i]);
		
		if (m.get_rank() == 0) {
			std::cout << "Test number " << i << " - ";
			
			if (delta <= eps) {
				std::cout << "OK" << std::endl;	
			}
			
			else {
				std::cout << "Error" << std::endl;
				std::cout << "Expected value = " << m.get_result().x << std::endl;
				std::cout << "Actual value = " << actual_value[i] << std::endl;
			}
			std::cout<<std::endl;
		}
	}
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
	
	test_parallel_version();

	MPI_Finalize();
	return 0;
}
