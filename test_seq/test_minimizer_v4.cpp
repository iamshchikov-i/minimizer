#include "minimizer_v4.h"
#include "functions.h"
#include "gtest.h"

const double eps = 0.01;
const int th_num = 2;

TEST(minimizer_v4, can_create_minimizer) {
	double(*fptr)(double) = f1;

	ASSERT_NO_THROW(Minimizer_v4 m(0, 10, fptr, th_num));
}

/*TEST(minimizer_v4, test_f1_1) {
	double delta;
	double(*fptr)(double) = f1;
	Minimizer_v4 m(-15, 8, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x - 5.46));

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f1_2) {
	double delta;
	double(*fptr)(double) = f1;
	Minimizer_v4 m(6.36, 20, fptr, th_num);
	m.solve();
	delta = abs(m.get_result().x - 6.36);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f1_3) {
	double delta;
	double(*fptr)(double) = f1;
	Minimizer_v4 m(-7.56, -1.92, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 1.92);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f2_1) {
	double delta;
	double(*fptr)(double) = f2;
	Minimizer_v4 m(-2, 5, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 0.5);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f2_2) {
	double delta;
	double(*fptr)(double) = f2;
	Minimizer_v4 m(-0.21, 10, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 0.21);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f2_3) {
	double delta;
	double(*fptr)(double) = f2;
	Minimizer_v4 m(-7.56, -1.92, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 1.92);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f3_1) {
	double delta;
	double(*fptr)(double) = f3;
	Minimizer_v4 m(0.12, 10.0, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 3.0);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f3_2) {
	double delta;
	double(*fptr)(double) = f3;
	Minimizer_v4 m(4.23, 10.0, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 4.23);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f3_3) {
	double delta;
	double(*fptr)(double) = f3;
	Minimizer_v4 m(-5.78, -1.92, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 5.78);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f3_4) {
	double delta;
	double(*fptr)(double) = f3;
	Minimizer_v4 m(-5.78, 5.78, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 5.78);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f4_1) {
	double delta;
	double(*fptr)(double) = f4;
	Minimizer_v4 m(-10.2, 10.2, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 1.0);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f4_2) {
	double delta;
	double(*fptr)(double) = f4;
	Minimizer_v4 m(0.12, 0.98, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 0.12);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f4_3) {
	double delta;
	double(*fptr)(double) = f4;
	Minimizer_v4 m(-19.2, -5.87, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 5.87);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f5_1) {
	double delta;
	double(*fptr)(double) = f5;
	Minimizer_v4 m(0.0, 1.0, fptr, th_num, 0.0001);
	m.solve();
	delta = abs(abs(m.get_result().x) - 0.8028);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f5_2) {
	double delta;
	double(*fptr)(double) = f5;
	Minimizer_v4 m(0.923, 1.338, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 0.923);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f5_3) {
	double delta;
	double(*fptr)(double) = f5;
	Minimizer_v4 m(-1.806, -1.448, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 1.448);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f6_1) {
	double delta;
	double(*fptr)(double) = f6;
	Minimizer_v4 m(0.1, 10.57, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 1.353);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f6_2) {
	double delta;
	double(*fptr)(double) = f6;
	Minimizer_v4 m(0.1, 1.298, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 1.298);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer_v4, test_f6_3) {
	double delta;
	double(*fptr)(double) = f6;
	Minimizer_v4 m(10.53, 19.23, fptr, th_num);
	m.solve();
	delta = abs(abs(m.get_result().x) - 10.53);

	EXPECT_EQ(1, delta <= eps);
}
*/