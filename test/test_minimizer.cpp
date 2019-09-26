#include "minimizer.h"

#include "gtest.h"

const double eps = 0.001;

double f1(double x) {
	return pow(x - 5.46, 2) + 9;
}

double f2(double x) {
	return 5 + x + pow(x, 2);
}

double f3(double x) {
	
	return pow(x, 3) - 6*pow(x, 2) + 9*x -4;
}

double f4(double x) {
	return (16.4*x) / (1 + pow(x, 2));
}

double f5(double x) {
	return sin(5.87*x);
}

double f6(double x) {
	return sqrt(x)*log(0.1*x);
}

TEST(minimizer, can_create_minimizer) {
	double(*fptr)(double) = f1;

	ASSERT_NO_THROW(Minimizer m(0, 10, fptr));
}

TEST(minimizer, test_f1_1) {
	double delta;
	double(*fptr)(double) = f1;
	Minimizer m(-15, 8, fptr);
	delta = abs(m.find_point() - 5.46);
	
	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f1_2) {
	double delta;
	double(*fptr)(double) = f1;
	Minimizer m(6.36, 20, fptr);
	delta = abs(m.find_point() - 6.36);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f1_3) {
	double delta;
	double(*fptr)(double) = f1;
	Minimizer m(-7.56, -1.92, fptr);
	delta = abs(abs(m.find_point()) - 1.92);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f2_1) {
	double delta;
	double(*fptr)(double) = f2;
	Minimizer m(-10, 5, fptr);
	delta = abs(abs(m.find_point()) - 0.5);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f2_2) {
	double delta;
	double(*fptr)(double) = f2;
	Minimizer m(-0.21, 10, fptr);
	delta = abs(abs(m.find_point()) - 0.21);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f2_3) {
	double delta;
	double(*fptr)(double) = f2;
	Minimizer m(-7.56, -1.92, fptr);
	delta = abs(abs(m.find_point()) - 1.92);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f3_1) {
	double delta;
	double(*fptr)(double) = f3;
	Minimizer m(0.12, 10.0, fptr);
	delta = abs(abs(m.find_point()) - 3.0);
	
	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f3_2) {
	double delta;
	double(*fptr)(double) = f3;
	Minimizer m(4.23, 10.0, fptr);
	delta = abs(abs(m.find_point()) - 4.23);
	
	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f3_3) {
	double delta;
	double(*fptr)(double) = f3;
	Minimizer m(-5.78, -1.92, fptr);
	delta = abs(abs(m.find_point()) - 5.78);
	
	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f3_4) {
	double delta;
	double(*fptr)(double) = f3;
	Minimizer m(-5.78, 5.78, fptr);
	delta = abs(abs(m.find_point()) - 5.78);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f4_1) {
	double delta;
	double(*fptr)(double) = f4;
	Minimizer m(-10.2, 10.2, fptr);
	delta = abs(abs(m.find_point()) - 1.0);
	
	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f4_2) {
	double delta;
	double(*fptr)(double) = f4;
	Minimizer m(0.12, 0.98, fptr);
	delta = abs(abs(m.find_point()) - 0.12);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f4_3) {
	double delta;
	double(*fptr)(double) = f4;
	Minimizer m(-19.2, -5.87, fptr);
	delta = abs(abs(m.find_point()) - 5.87);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f5_1) {
	double delta;
	double(*fptr)(double) = f5;
	Minimizer m(0.0, 1.0, fptr, 0.0001);
	delta = abs(abs(m.find_point()) - 0.8028);
	
	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f5_2) {
	double delta;
	double(*fptr)(double) = f5;
	Minimizer m(0.923, 1.338, fptr);
	delta = abs(abs(m.find_point()) - 0.923);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f5_3) {
	double delta;
	double(*fptr)(double) = f5;
	Minimizer m(-1.806, -1.448, fptr);
	delta = abs(abs(m.find_point()) - 1.448);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f6_1) {
	double delta;
	double(*fptr)(double) = f6;
	Minimizer m(0.1, 10.57, fptr);
	delta = abs(abs(m.find_point()) - 1.353);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f6_2) {
	double delta;
	double(*fptr)(double) = f6;
	Minimizer m(0.1, 1.298, fptr);
	delta = abs(abs(m.find_point()) - 1.298);

	EXPECT_EQ(1, delta <= eps);
}

TEST(minimizer, test_f6_3) {
	double delta;
	double(*fptr)(double) = f6;
	Minimizer m(10.53, 19.23, fptr);
	delta = abs(abs(m.find_point()) - 10.53);

	EXPECT_EQ(1, delta <= eps);
}