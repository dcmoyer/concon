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


#include "sigma_opt.hpp"
namespace sigma_opt{
double calc_half_kernel(const double dot_prod, const int& L_resolution,
  const int& num_harmonics, double** harmonic_lookup,
  const double& sigma, double * exp_lookup){
  double ret_val = 0;
  if(exp_lookup == NULL){
    for(int i = 0; i < num_harmonics; ++i){
      //WE NEED A LOOKUP TABLE!
      ret_val += exp(-sigma * i * (i+1)) * (2.0 * i + 1.0) *
        harmonic_lookup[i][(int) (dot_prod * L_resolution + L_resolution)];
      //4 * M_PI factor removed due to unit sphere normalization
    }
  } else { 
    for(int i = 0; i < num_harmonics; ++i){
      ret_val += exp_lookup[i] * 
        harmonic_lookup[i][(int) (dot_prod * L_resolution + L_resolution)];
    }
  }
  return fmax(ret_val,1e-20);
  //return fmax(ret_val,std::numeric_limits<double>::denorm_min());
}

double calc_kernel(const int& num_harmonics, double** harmonics,
  const double& sigma, double ** exp_lookup){
  double ret_val = 0;
  if(exp_lookup == NULL){
    for(int i = 0; i < num_harmonics; ++i){
      for(int j = 0; j < num_harmonics; ++j){
      //WE NEED A LOOKUP TABLE!
        ret_val += exp(-sigma *( i * (i+1) + j * (j + 1))) *
          (2 * i + 1) * (2 * j + 1) *
          harmonics[i][j];
        //4 * M_PI factor removed due to unit sphere normalization
      }
    }
  } else { 
    for(int i = 0; i < num_harmonics; ++i){
      for(int j = 0; j < num_harmonics; ++j){
        ret_val += exp_lookup[i][j] * harmonics[i][j];
      }
    }
  }
  return fmax(ret_val,1e-10);
  //return fmax(ret_val,std::numeric_limits<double>::denorm_min());
}

double calc_kernel_dsigma(const int& num_harmonics, double** harmonics, const double& sigma){

  double ret_val = 0;
  for(int i = 0; i < num_harmonics; ++i){
    for(int j = 0; j < num_harmonics; ++j){
      //WE NEED A LOOKUP TABLE!
      ret_val += -sigma *( i * (i+1) + j * (j + 1)) *
        exp(-sigma *( i * (i+1) + j * (j + 1))) *
        (2 * i + 1) * (2 * j + 1) *
        harmonics[i][j];
      //4 * M_PI factor removed due to unit sphere normalization
    }
  }
  return ret_val;

}

double grad_descent_for_sigma_with_null_hyp(
  const int& num_test_points,const int& num_test_fibers,
  const int& num_harmonics, double*** tf_harmonics, double *** tp_harmonics,
  double epsilon, int max_iters,double weight){

  //Switched Optimal calcuation to Jones Marron Sheather
  // using the Integrated Squared Error
  // e.g. \int \hat{f}^2 - \sum_{i}^N LeaveOneOut.

  double sigma = 0.01;
  double alpha = 1.0/64.0 / 2.0;
  double** saved_calcs_points = new double*[2];
  saved_calcs_points[0] = new double[num_test_points];
  saved_calcs_points[1] = new double[num_test_points];

  double** saved_calcs_fibers = new double*[2];
  saved_calcs_fibers[0] = new double[num_test_fibers];
  saved_calcs_fibers[1] = new double[num_test_fibers];


  double prev_value = 0;
  double test_value = 0;
  double gradient = 0;
  for(int i = 0; i < num_test_points; ++i){
    saved_calcs_points[0][i] += -weight*
      log(calc_kernel(num_harmonics, tp_harmonics[i], sigma));
    prev_value += saved_calcs_points[0][i];
  }

  for(int i = 0; i < num_test_fibers; ++i){
    saved_calcs_fibers[0][i] += (1-weight)*
      log(calc_kernel(num_harmonics, tf_harmonics[i], sigma));
    prev_value += saved_calcs_fibers[0][i];
  }

  //the descent
  int iter_count = 0;
  do{
    alpha = 2*alpha;
    gradient = 0;

    for(int i = 0; i < num_test_points; ++i){
      gradient += -weight*
        calc_kernel_dsigma(num_harmonics, tp_harmonics[i], sigma) /
        calc_kernel(num_harmonics, tp_harmonics[i], sigma);
    }
    for(int i = 0; i < num_test_fibers; ++i){
      gradient += (1-weight)*
        calc_kernel_dsigma(num_harmonics, tf_harmonics[i], sigma) /
        calc_kernel(num_harmonics, tf_harmonics[i], sigma);
    }
    gradient = gradient / (num_test_points + num_test_fibers);

    //find new step size
    do{
      test_value = 0;

      //check for lower bound
      if(sigma + alpha * gradient < 0){
        alpha = alpha/2.0;
        continue;
      }
      for(int i = 0; i < num_test_points; ++i){
        saved_calcs_points[1][i] = -weight*
          log(calc_kernel(num_harmonics,tp_harmonics[i],sigma + alpha * gradient));
        test_value += saved_calcs_points[1][i];
      }

      for(int i = 0; i < num_test_fibers; ++i){
        saved_calcs_fibers[1][i] += (1-weight)*log(
          calc_kernel(num_harmonics,tf_harmonics[i],sigma + alpha * gradient));
        test_value += saved_calcs_fibers[1][i];
      }

      //if test_value better than prev_value
      if(test_value > prev_value){
        break;
      } else {
        alpha = alpha/2.0;
      }
    }
    while(alpha > epsilon);

    sigma = sigma + alpha*gradient;

    //swap the calcs
    double* temp = saved_calcs_points[0];
    saved_calcs_points[0] = saved_calcs_points[1];
    saved_calcs_points[1] = temp;

    temp = saved_calcs_fibers[0];
    saved_calcs_fibers[0] = saved_calcs_fibers[1];
    saved_calcs_fibers[1] = temp;

    prev_value = test_value;
    if(iter_count > max_iters){
      std::cout << "Warning! Sigma optimization passed max iters." << std::endl;
      break;
    }
    iter_count += 1;
    //std::cout << sigma << " : " << test_value << std::endl;
  }
  while(alpha > epsilon );

  std::cout << "sigma optimization done in " << iter_count << " iterations"
    << std::endl;

  delete[] saved_calcs_points[0];
  delete[] saved_calcs_points[1];
  delete[] saved_calcs_points;
  delete[] saved_calcs_fibers[0];
  delete[] saved_calcs_fibers[1];
  delete[] saved_calcs_fibers;

  return sigma;

}



double grad_descent_for_sigma(const int& num_test_points,
  const int& num_harmonics, double*** harmonics, double epsilon,
  int max_iters){

  //Switched Optimal calcuation to Jones Marron Sheather
  // using the Integrated Squared Error
  // e.g. \int \hat{f}^2 - \sum_{i}^N LeaveOneOut.

  double sigma = 0.1;
  double alpha = 1.0/64.0 / 2.0;
  double** saved_calcs = new double*[2];
  saved_calcs[0] = new double[num_test_points];
  saved_calcs[1] = new double[num_test_points];

  double prev_value = 0;
  double test_value = 0;
  double gradient = 0;
  for(int i = 0; i < num_test_points; ++i){
    saved_calcs[0][i] = calc_kernel(num_harmonics, harmonics[i], sigma);
    if(saved_calcs[0][i] < epsilon)
      std::cout << saved_calcs[0][i] << std::endl;
    prev_value += log(saved_calcs[0][i]);
  }

  //the descent
  int iter_count = 0;
  do{
    alpha = 2*alpha;
    gradient = 0;

    for(int i = 0; i < num_test_points; ++i){
      gradient += calc_kernel_dsigma(num_harmonics, harmonics[i], sigma) /
        saved_calcs[0][i];
        //calc_kernel(num_harmonics, harmonics[i], sigma);
    }
    //gradient = gradient / num_test_points;
    gradient = gradient / fabs(gradient);

    std::cout << "value is " << prev_value << "\tgradient is " << gradient
      << std::endl;
    //find new step size
    do{
      test_value = 0;

      //check for lower bound
      if(sigma + alpha * gradient < 0){
        alpha = alpha/2.0;
        continue;
      }
      for(int i = 0; i < num_test_points; ++i){
        saved_calcs[1][i] = calc_kernel(num_harmonics,harmonics[i],sigma + alpha * gradient);
        test_value += log(saved_calcs[1][i]);
      }

      //if test_value better than prev_value
      if(test_value > prev_value){
        break;
      } else {
        alpha = alpha/2.0;
      }
    }
    while(alpha > epsilon);

    sigma = sigma + alpha*gradient;

    //swap the calcs
    double* temp = saved_calcs[0];
    saved_calcs[0] = saved_calcs[1];
    saved_calcs[1] = temp;

    prev_value = test_value;
    if(iter_count > max_iters){
      std::cout << "Warning! Sigma optimization passed max iters." << std::endl;
      break;
    }
    iter_count += 1;
    //std::cout << sigma << " : " << test_value << std::endl;
  }
  while(alpha > epsilon );

  std::cout << "sigma optimization done in " << iter_count << " iterations"
    << std::endl;

  delete[] saved_calcs[0];
  delete[] saved_calcs[1];
  delete[] saved_calcs;

  return sigma;

}



double grad_descent_for_sigma_ISE(const int& num_test_points,
  const int& num_harmonics, double*** tp_harmonics, double ***tf_harmonics, double epsilon,
  int max_iters){

  //Switched Optimal calcuation to Jones Marron Sheather
  // using the Integrated Squared Error
  // e.g. \int \hat{f}^2 - \sum_{i}^N LeaveOneOut.

  double sigma = 1;
  double alpha = 1.0/64.0 / 2.0;
  double** saved_calcs = new double*[2];
  saved_calcs[0] = new double[num_test_points];
  saved_calcs[1] = new double[num_test_points];

  double prev_value = 0;
  double test_value = 0;
  double temp_kern = 0;
  double gradient = 0;
  //for(int i = 0; i < num_test_points; ++i){
  //  saved_calcs[0][i] = calc_kernel(num_harmonics, tp_harmonics[i], sigma);
  //  prev_value += log(saved_calcs[0][i]);
  //}

  //the descent
  int iter_count = 0;
  do{
    alpha = 2*alpha;
    gradient = 0;

    for(int i = 0; i < num_test_points; ++i){
      gradient += 2*calc_kernel_dsigma(num_harmonics, tp_harmonics[i],sigma)-
        calc_kernel_dsigma(num_harmonics, tf_harmonics[i],sigma);
      //gradient += calc_kernel_dsigma(num_harmonics, harmonics[i], sigma) /
      //  saved_calcs[0][i];
        //calc_kernel(num_harmonics, harmonics[i], sigma);
    }
    gradient = gradient / num_test_points;

    //find new step size
    do{
      test_value = 0;

      //check for lower bound
      if(sigma - alpha * gradient < 0){
        alpha = alpha/2.0;
        continue;
      }
      for(int i = 0; i < num_test_points; ++i){
        temp_kern = calc_kernel(num_harmonics, tp_harmonics[i],sigma - alpha * gradient);
        saved_calcs[1][i] = temp_kern * temp_kern -
          calc_kernel(num_harmonics, tf_harmonics[i],sigma - alpha * gradient);
        test_value += saved_calcs[1][i];
      }

      //if test_value better than prev_value
      if(test_value < prev_value){
        break;
      } else {
        alpha = alpha/2.0;
      }
    }
    while(alpha > epsilon);

    sigma = sigma - alpha*gradient;

    //swap the calcs
    double* temp = saved_calcs[0];
    saved_calcs[0] = saved_calcs[1];
    saved_calcs[1] = temp;

    prev_value = test_value;
    if(iter_count > max_iters){
      std::cout << "Warning! Sigma optimization passed max iters." << std::endl;
      break;
    }
    iter_count += 1;
    //std::cout << sigma << " : " << test_value << std::endl;
  }
  while(alpha > epsilon );

  std::cout << "sigma optimization done in " << iter_count << " iterations"
    << std::endl;

  delete[] saved_calcs[0];
  delete[] saved_calcs[1];
  delete[] saved_calcs;

  return sigma;

}

double exp_grid_search_sigma_ISE(const int& num_test_points, const int& num_test_fibers,
  const int& num_harmonics, double*** tp_harmonics, double ***tf_harmonics,
  double interval_min, double interval_max, double delta, bool verbose){

  double min_fvalue = 10000000;
  double min_x = -1;
  double test_fvalue;
  double temp_kernel;
  double sigma;

  for(double x = interval_min; x < interval_max; x += delta){
    sigma = exp(x);
    test_fvalue = 0;
    for(int i = 0; i < num_test_points; ++i ){
      temp_kernel = calc_kernel(num_harmonics, tp_harmonics[i],sigma);
      //test_fvalue += (4.0 * M_PI/num_test_points) * temp_kernel * temp_kernel;
      test_fvalue += (1.0/num_test_points) * temp_kernel * temp_kernel;
    }
    for(int i = 0; i < num_test_fibers; ++i ){
      test_fvalue -= 2*(1.0/num_test_fibers)*
        log(calc_kernel(num_harmonics, tf_harmonics[i],sigma));
    }

    if(test_fvalue < min_fvalue){
      min_fvalue = test_fvalue;
      min_x = sigma;
      printf("new best: %f : %f\n", sigma, test_fvalue);
    } else if (verbose == true ){
      printf("sigma val: %f : %f\n",sigma, test_fvalue);
    }

  }

  return min_x;

}

double grid_search_sigma_ISE(const int& num_test_points, const int& num_test_fibers,
  const int& num_harmonics, double*** tp_harmonics, double ***tf_harmonics,
  double interval_min, double interval_max, double delta, bool verbose){

  double min_fvalue = 10000000;
  double min_x = -1;
  double test_fvalue;
  double temp_kernel;

  for(double x = interval_min; x < interval_max; x += delta){
    test_fvalue = 0;
    for(int i = 0; i < num_test_points; ++i ){
      temp_kernel = calc_kernel(num_harmonics, tp_harmonics[i],x);
      test_fvalue += (1.0/num_test_points) * temp_kernel * temp_kernel;
    }
    for(int i = 0; i < num_test_fibers; ++i ){
      test_fvalue -= 2*(1.0/num_test_fibers)*
        log(calc_kernel(num_harmonics, tf_harmonics[i],x));
    }

    if(test_fvalue < min_fvalue){
      min_fvalue = test_fvalue;
      min_x = x;
      printf("new best: %f : %f\n", x, test_fvalue);
    } else if (verbose == true ){
      printf("sigma val: %f : %f\n",x, test_fvalue);
    }

  }

  return min_x;

}

}//end of namespace


