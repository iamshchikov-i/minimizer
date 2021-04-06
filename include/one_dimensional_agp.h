#ifndef __ONE_DIMENSIONAL_AGP_H__
#define __ONE_DIMENSIONAL_AGP_H__

#include "one_dimensional_minimizer.h"

class One_Dimensional_AGP : public One_Dimensional_Minimizer {
protected:
	std::priority_queue<interval, std::vector<interval>, CompareR_min>* pq;
	virtual bool isEnd();
	virtual double get_M(), get_m(), get_R(), get_new_point(interval i);
	virtual void compute_R(std::vector<double> new_point, double new_m),
		insert_to_map(std::vector<double> res_coords, double _z, double _R,
			double _num_estimation = 0.0),
		compare_M(std::vector<double> new_point), perform_first_iteration(), delete_containers();
	void perform_first_parallel_step(std::pair<double, double>& new_point, result& tmp_res, double& new_m,
		std::map<int, std::pair<double, double>>& processing_intervals_map,
		std::set<std::pair<double, double>>& processing_intervals_set, MPI_Status& st);
	void do_parallel_job(double last_coord,
		std::vector<double>& res, int rank);
public:
	 void set_experiment(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
		 std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
		 bool useMPI, bool useThreads, int threadsNum,
		 double(*f)(std::vector<double> x),
		 double _eps = 0.001, int _Nmax = 1000, double _r_par = 2.0);
	 One_Dimensional_AGP(int _range, int _curr_dim, std::vector<One_Dimensional_Minimizer*> _odm,
		 std::vector<std::pair<double, double>> _bounds, std::vector<double> _curr_x,
		 bool useMPI, bool useThreads, int threadsNum,
		double(*f)(std::vector<double> x),
		double _eps = 0.001, int _Nmax = 1000, double _r_par = 2.0);
	virtual result solve_mpi();
	virtual result solve_seq();
};

#endif  
