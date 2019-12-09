#include <thread>
#include <chrono>
#include "minimizer_v1.h"
#include "minimizer_v2.h"
#include "minimizer_v3.h"
#include "functions.h"

void set_bords(std::vector<double>&);

void set_actual_values(std::vector<double>&);

void time();

void measure_time_of_function(double(*f)(double x));

void measure_time_v2_v3();

void pr(int argc, char** argv) {
	for (int i = 0; i < argc;i++)
		std::cout << argv[i] << std::endl;
}

int main(int argc, char** argv) {
	/*MPI_Init(&argc, &argv);
	MPI_Finalize();*/

	pr(argc, argv);
	
	return 0;
}

void time() {
	int test_num = 18;
	
	std::chrono::milliseconds r1 = std::chrono::milliseconds::zero(),
							  r2 = std::chrono::milliseconds::zero(),
						      r3 = std::chrono::milliseconds::zero(),
	elapsed_ms;
	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	std::vector<double> bord(test_num * 2);
	std::vector<double> actual_value(test_num);
	double(*fptr[6])(double) = { f1, f2, f3, f4, f5, f6 };
	const double eps = 0.01;
	double delta = 0.0;
	double res;
	set_bords(bord);
	set_actual_values(actual_value);
	
	Minimizer_v3 m3(0, 0, nullptr);

	if (m3.get_rank() == 0) {

		std::cout << "Via Minimizer_v1" << std::endl;
		
		for (int i = 0;i < test_num;++i) {
			Minimizer_v1 m1(bord[2 * i], bord[(2 * i) + 1], fptr[i / 3]);
			begin = std::chrono::steady_clock::now();
			res = m1.find_point();
			end = std::chrono::steady_clock::now();
			elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
			std::cout << "    Test number " << i << std::endl;
			std::cout << "		time =  " << elapsed_ms.count() << std::endl;
			std::cout << "		number of point =  " << m1.get_k() << std::endl;
			std::cout << "		x =  " << res << std::endl;
			r1 += elapsed_ms;
		}
		std::cout << "Total time 1: " << r1.count() << " ms\n";
	}

	if (m3.get_rank() == 0) {
		
		std::cout << "Via Minimizer_v2" << std::endl;
		Minimizer_v2 m(0, 0, nullptr);
		for (int i = 0;i < test_num;++i) {
			m.set_experiment(bord[2 * i], bord[(2 * i) + 1], fptr[i / 3]);
			begin = std::chrono::steady_clock::now();
			m.solve();
			end = std::chrono::steady_clock::now();
			elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
			std::cout << "    Test number " << i << std::endl;
			std::cout << "		time = " << elapsed_ms.count() << std::endl;
			std::cout << "		number of point =  " << m.get_result().k << std::endl;
			std::cout << "		x =  " << m.get_result().x << std::endl;
			r2 += elapsed_ms;
		}
		std::cout << "Total time 2: " << r2.count() << " ms\n";
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (m3.get_rank() == 0)
		std::cout << "Via Minimizer_v3" << std::endl;
	
	for (int i = 0;i < test_num;++i) {
		m3.set_experiment(bord[2 * i], bord[(2 * i) + 1], fptr[i / 3]);

		MPI_Barrier(MPI_COMM_WORLD);
		begin = std::chrono::steady_clock::now();
		m3.solve();
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

		if (m3.get_rank() == 0) {
			std::cout << "    Test number " << i << std::endl;
			std::cout << "		time " << elapsed_ms.count() << std::endl;
			std::cout << "		number of point =  " << m3.get_result().k << std::endl;
			std::cout << "		x =  " << m3.get_result().x << std::endl;
			r3 += elapsed_ms;
		}
	}
	if (m3.get_rank() == 0)
		std::cout << "Total time 3: " << r3.count() << " ms\n";
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
	result res;
	std::pair<double, double> bord;
	std::chrono::milliseconds elapsed_ms;
	std::chrono::time_point<std::chrono::steady_clock> begin, end;
	bord = { -15.0, 8.0 };
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (rank == 0) {
		std::cout << std::endl << "Via Minimizer_v2: " << std::endl;
		Minimizer_v2 m2(bord.first, bord.second, f1);
		begin = std::chrono::steady_clock::now();
		res = m2.solve();
		end = std::chrono::steady_clock::now();
		elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
		std::cout << "Number of points = " << res.k << std::endl;
		std::cout << "Total time 2: " << elapsed_ms.count() << " ms\n";
	}

	MPI_Barrier(MPI_COMM_WORLD);
	Minimizer_v3 m3(bord.first, bord.second, f1);
	begin = std::chrono::steady_clock::now();
	m3.solve();
	end = std::chrono::steady_clock::now();
	elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	if (m3.get_rank() == 0) {
		std::cout << std::endl << "Via Minimizer_v3: " << std::endl;
		std::cout << "Number of points = " << m3.get_result().k << std::endl;
		std::cout << "Total time 3: " << elapsed_ms.count() << " ms\n";
	}
}