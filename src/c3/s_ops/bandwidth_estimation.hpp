
///
/// Bandwidth Estimation Operation
///
/// This class attempts to estimate sigma for a specific subject. It will
/// produce output to a file if specified. -1 num_points denotes the use of
/// all the points.
///
/// Optimizations:
///   -Moving the num_harmonics variable into a constant Subject variable
///   makes loop unrolling very effective (~10x speed up).
///   -I also moved the allocation and de-alocation of two arrays to class
///   memory. This saves something like 1% run time, but hey, it works.

#ifndef C3_Bandwidth_Estimation_HPP
#define C3_Bandwidth_Estimation_HPP

#include<iomanip>
#include "c3/subject.hpp"

namespace c3 { 

class Bandwidth_Estimation : public c3::Subject {
public:
  Bandwidth_Estimation(); //blank
  ~Bandwidth_Estimation(); //blank

  int set(std::unordered_map<std::string, std::string> params);

  int subj_specific_load(std::string subj){
    return 0;
  }

  int run();

  int save_file(std::string filename);

  int load_file(std::string filename);
 
  bool check_req();

protected:
  //-1 denotes all points
  bool monte_carlo_flag;
  bool rand_selection;

  int num_test_points;
  int num_test_fibers;
  double range_min;
  double range_max;
  int range_size;
  double range_step;

  double* sigma_objective;

private:
  void monte_carlo_integrator();
  void harmonic_helper_func(
    int* idxA, float* coordsA, //test
    int* idxB, float* coordsB, //data
    double* px_term,
    double** test_fiber_cross_harm);

  //this calculates
  // $\sigma (i + i^2 + j + j^2)\frac{(2i + 1)(2j + 1)}{16\pi^2}$
  //for each value of sigma in the range
  //indices are [sigma_idx][i][j]
  void calc_exp_lookup();
  double *** exp_lookup;

  double *temp_harmonics_0;
  double *temp_harmonics_1;

  //this checks whether or not we should filp the coordinates,
};

}//end of namespace c3

#endif

