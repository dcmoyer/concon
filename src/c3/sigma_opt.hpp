//Copyright 2017 Daniel Moyer
//
//Permission is hereby granted, free of charge, to any person obtaining a
//copy of this software and associated documentation files (the "Software"),
//to deal in the Software without restriction, including without limitation
//the rights to use, copy, modify, merge, publish, distribute, sublicense,
//and/or sell copies of the Software, and to permit persons to whom the
//Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in
//all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
//OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
//OTHER DEALINGS IN THE SOFTWARE.
//



#ifndef SIGMA_OPT_HPP
#define SIGMA_OPT_HPP

#include <fstream>
#include <float.h>
#include <string>
#include <math.h>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <limits>

#include "c3/check_flip.h"
#include "sh/spherical_harmonics.h"
#include <memory>
#include <array>
#include <functional>
#include "MeshLib/Solid.h"
//#include <Eigen/Sparse>


namespace sigma_opt {

double calc_half_kernel(const double dot_prod, const int& L_resolution,
  const int& num_harmonics, double** harmonic_lookup,
  const double& sigma, double * exp_lookup=NULL);

double calc_kernel(const int& num_harmonics, double** harmonics,
  const double& sigma, double ** exp_lookup=NULL);
double calc_kernel_dsigma(const int& num_harmonics, double** harmonics,
  const double& sigma);


double grad_descent_for_sigma_with_null_hyp(
  const int& num_test_points, const int& num_test_fibers,
  const int& num_harmonics, double*** tf_harmonics, double *** tp_harmonics,
  double epsilon, int max_iters=10000,double weight=0.1);

double grad_descent_for_sigma(const int& num_test_points,
  const int& num_harmonics, double*** harmonics, double epsilon,
  int max_iters=10000);

double grad_descent_for_sigma_ISE(const int& num_test_points,
  const int& num_harmonics, double*** tp_harmonics, double ***tf_harmonics, double epsilon,
  int max_iters=10000);

double exp_grid_search_sigma_ISE(const int& num_test_points, const int& num_test_fibers,
  const int& num_harmonics, double*** tp_harmonics, double ***tf_harmonics,
  double interval_min=0.01, double interval_max=0.5, double delta=0.01,
  bool verbose=false);
double grid_search_sigma_ISE(const int& num_test_points, const int& num_test_fibers,
  const int& num_harmonics, double*** tp_harmonics, double ***tf_harmonics,
  double interval_min=0.01, double interval_max=0.5, double delta=0.01,
  bool verbose=false);

}

#endif

