
#include "c3/s_ops/bandwidth_estimation.hpp"

namespace c3 {

Bandwidth_Estimation::Bandwidth_Estimation() : c3::Subject() {

  r_trk_endpoints = true;
  r_grid = true;
  r_harm_lookup_table = true;

  sigma_objective = NULL;
  exp_lookup = NULL;

  temp_harmonics_0 = new double[num_harmonics];
  temp_harmonics_1 = new double[num_harmonics];
}

Bandwidth_Estimation::~Bandwidth_Estimation(){
  rand_selection = false;
  monte_carlo_flag = false;

  delete[] temp_harmonics_0;
  delete[] temp_harmonics_1;

  if(sigma_objective != NULL)
    delete[] sigma_objective;

  if(exp_lookup != NULL){
    for(int i = 0; i < range_size; ++i){
      for(int j = 0; j < num_harmonics; ++j)
        delete[] exp_lookup[i][j];
      delete[] exp_lookup[i];
    }
    delete[] exp_lookup;
  }

}

int Bandwidth_Estimation::set(
  std::unordered_map<std::string, std::string> params){

  bool num_points_used = false;
  bool num_fibers_used = false;
  bool range_min_used = false;
  bool range_max_used = false;
  bool range_size_used = false;
  monte_carlo_flag = false;  

  for( std::pair<std::string,std::string> p : params ){

    if(p.first == "num_test_points"){
      num_test_points = c3::string_to_int(p.second);
      num_points_used = true;
    }

    if(p.first == "num_test_fibers"){
      num_test_fibers = c3::string_to_int(p.second);
      num_fibers_used = true;
    }

    if(p.first == "rand_selection"){
      rand_selection = true;
    }

    if(p.first == "verbose"){
      verbose = true;
    }

    if(p.first == "monte_carlo_flag"){
      monte_carlo_flag = true;
    }

    if(p.first == "range_min"){
      range_min = c3::string_to_double(p.second);
      range_min_used = true;
    }

    if(p.first == "range_max"){
      range_max = c3::string_to_double(p.second);
      range_max_used = true;
    }

    if(p.first == "range_size"){
      range_size = c3::string_to_int(p.second);
      range_size_used = true;
    }

  }

  if((!num_points_used && monte_carlo_flag) || !num_fibers_used ||
    !range_min_used || !range_max_used || !range_size_used)
    return 1;
  else
    return 0;

}

bool Bandwidth_Estimation::check_req(){

  if(trk_endpoints_idx == NULL){
    std::cerr << "Error: trk endpoint file not loaded. Aborting"
    << std::endl;
    return false;
  }

  if(!monte_carlo_flag && (lh_grid == NULL || rh_grid == NULL)){
    std::cerr << "Error: integration grid not set. Aborting." << std::endl;
    return false;
  }

  if(harm_lookup_table == NULL){
    std::cerr << "Error: Harmonic Lookup Table not set. Aborting."
      << std::endl;
    return false;
  }

  return true;
}

int Bandwidth_Estimation::run(){

  if(!check_req()){
    return 1;
  }

  //
  //  Computes the Leave-One-Out prob as well as the integrated squared error.
  //

  //double sum_of_errors = 0;
  //double max_error = 0;
  //double temp_harmonics_0[num_harmonics];
  //double temp_harmonics_1[num_harmonics];
  int test_idx, local_N;

  sigma_objective = new double[range_size];
  for(int i = 0; i < range_size; ++i){
    sigma_objective[i] = 0;
  }
  range_step = ( range_max - range_min )/((double) range_size-1.0);

  double ** test_fiber_cross_harm = new double*[num_harmonics];
  for(int i = 0; i < num_harmonics; ++i){
    test_fiber_cross_harm[i] = new double[num_harmonics];
  }

  if(verbose)
    std::cout << "computing coefs for test points" << std::endl;

  time_t start_kernel_calc = time(0);

  //THIS IS WHERE WE'D OPENMP
  //for each track

  if(num_test_fibers == -1){
    rand_selection = true;
    local_N = num_tracks;
  } else {
    rand_selection = true;
    local_N = num_test_fibers;
  }

  if(verbose)
    std::cout << "Calculating exp lookup" << std::endl;
  calc_exp_lookup();

  if(verbose)
    std::cout << "Switched to Penalized Log Like" << std::endl;
  double px_term[2];
  for(int cv_idx = 0; cv_idx < local_N; ++cv_idx){

    if(rand_selection)
      test_idx = rand() % num_tracks;
    else
      test_idx = cv_idx;

    //if invalid point
    if(trk_endpoints_idx[test_idx][0] == -1 ||
      trk_endpoints_idx[test_idx][2] == -1){

      if(rand_selection)
        cv_idx--;
      continue;
    }

    //zero out the test_fiber_cross_harm
    for(int i = 0; i < num_harmonics; ++i)
      for(int j = 0; j < num_harmonics; ++j)
        test_fiber_cross_harm[i][j] = 0;

    for(int data_idx = 0; data_idx < num_tracks; ++data_idx) {

      if(data_idx == test_idx)
        continue;

      harmonic_helper_func(
        trk_endpoints_idx[test_idx], trk_endpoints_coords[test_idx],
        trk_endpoints_idx[data_idx], trk_endpoints_coords[data_idx],
        px_term,
        test_fiber_cross_harm);

    }

    //if this is slow
    //make the exp_lookup
    //This is
    //$\mathcal{L}_{D-p}(p) = \frac{1}{|D - p|}\sum_{q \in D-p} \kappa(p,q)$
    //added log term
    //TODO: make the log into an injectable form
    for(int sig_idx = 0; sig_idx < range_size; ++sig_idx)
      sigma_objective[sig_idx] += (1.0/(num_nonnull_tracks - 1)) * 
        log(sigma_opt::calc_kernel(num_harmonics,
        test_fiber_cross_harm,
        range_min + range_step*sig_idx,
        exp_lookup[sig_idx]));

    //timing
    if(verbose){
      std::cout << std::setw(40) << "\raverage time per fiber: "
        << std::setw(10) << std::setprecision(6) << std::fixed
        << ((double) (time(0) - start_kernel_calc))/(cv_idx + 1)
        << std::flush;
    }
  }

  if(verbose){
    std::cout << std::setw(40) << "\raverage time per fiber: "
      << std::setw(10) << std::setprecision(6) << std::fixed
      << ((double) (time(0) - start_kernel_calc))/ (local_N)
      << std::endl;
  }
  std::cout << std::setw(40) << "leave one out coef calc runtime: "
    << std::setw(10) << time(0) - start_kernel_calc << std::endl;

  //$(2.0/|D|) \sum_{p\in D} \mathcal{L}_{D - p}(p)$
  for(int sig_idx = 0; sig_idx < range_size; ++sig_idx)
    sigma_objective[sig_idx] *= -(2.0/num_nonnull_tracks);

  //
  //  Integration for $\int \hat{f}(x)^2 dx$
  //
  int mc_hemi[2];
  float mc_position[6];
  double temp_kernel = 0;
  time_t start_integ_calc;
  int num_point_pairs;
  int grid_size1, grid_size2;
  double integrator_block_size;

  if(monte_carlo_flag){
    monte_carlo_integrator();
  } else {

    //rh to rh
    //lh to lh
    //rh to lh
    //lh to rh (skipped)
    MeshLib::Solid* grid1, *grid2;
    MeshLib::Vertex *v1,*v2;
    MeshLib::Point p1,p2;
    for(int t = 0; t < 4; ++t){ 

      if(t == 0){
        grid1 = grid2 = rh_grid;
        mc_hemi[0] = mc_hemi[1] = 0;
      } else if (t == 1) {
        grid1 = grid2 = lh_grid;
        mc_hemi[0] = mc_hemi[1] = 0;
      } else if (t == 2) {
        grid1 = rh_grid;
        grid2 = lh_grid;
        mc_hemi[0] = 0;
        mc_hemi[1] = 1;
      } else {
        continue;
      }

      grid_size1 = grid1->numVertices();
      grid_size2 = grid2->numVertices();

      num_point_pairs = 0;
      start_integ_calc = time(0);
      for(MeshLib::SolidVertexIterator viter1(grid1);
        !viter1.end(); ++viter1) {

        v1 = *viter1;
        p1 = v1->point();

        MeshLib::SolidVertexIterator viter2(grid2);
        if(t < 2){
          viter2 = viter1;
          ++viter2;
          integrator_block_size = (2.0/(grid_size1 * grid_size2));
        } else {
          viter2 = MeshLib::SolidVertexIterator(grid2);
          integrator_block_size = (1.0/(grid_size1 * grid_size2));
        }

        for(; !viter2.end(); ++viter2){
          v2 = *viter2;
          p2 = v2->point();

          //zero out the test_fiber_cross_harm
          for(int i = 0; i < num_harmonics; ++i)
            for(int j = 0; j < num_harmonics; ++j)
              test_fiber_cross_harm[i][j] = 0;

          for(int pos = 0; pos < 3; ++pos){
            //someone uses 1 indexing.
            //someone's in trouble yo.
            mc_position[pos] = p1[pos];
            mc_position[pos + 3] = p2[pos];
          }

          for(int data_idx = 0; data_idx < num_tracks; ++data_idx) {
            
            //if invalid point
            if(trk_endpoints_idx[data_idx][0] == -1 ||
              trk_endpoints_idx[data_idx][2] == -1){
              continue;
            }

            harmonic_helper_func(
              mc_hemi, mc_position,
              trk_endpoints_idx[data_idx], trk_endpoints_coords[data_idx],
              px_term,
              test_fiber_cross_harm);

          }//end data loop
 
          //if this is slow
          //make the exp_lookup
          //This is
          //$\mathcal{L}_{D-p}(p) = \frac{1}{|D - p|}\sum_{q \in D-p} \kappa(p,q)$
          for(int sig_idx = 0; sig_idx < range_size; ++sig_idx){
              temp_kernel = sigma_opt::calc_kernel(num_harmonics,
                test_fiber_cross_harm,
                range_min + range_step*sig_idx,
                exp_lookup[sig_idx]);

              //note: we only integrated half the domain because of its
              //  symmetry. Thus, x2. Also, full surface area = 1.
              //  However, further recall that we need to divide by a further
              //  \frac{1}{n} for the cross validation window
              //  See Hall and Marron 87
              sigma_objective[sig_idx] +=
                integrator_block_size * temp_kernel * temp_kernel;
          }
 
          //timing
          if(verbose){
            num_point_pairs += 1;
            std::cout << std::setw(40) << "\raverage time per point: "
              << std::setw(10) << std::setprecision(6) << std::fixed
              << ((double) (time(0) - start_integ_calc))/(num_point_pairs)
              << std::flush;
          }
        }//end v2 loop



      }//end v1 loop
    
      if(verbose){
        std::cout << std::setw(40) << "\raverage time per fiber: "
          << std::setw(10) << std::setprecision(6) << std::fixed
          << ((double) (time(0) - start_integ_calc))/ num_point_pairs
          << std::endl;
      }
      std::cout << std::setw(40) << "integrated square func calc runtime: "
        << std::setw(10) << time(0) - start_integ_calc << std::endl;

    }//end hemi loop
  }

  for(int i = 0; i < num_harmonics; ++i){
    delete[] test_fiber_cross_harm[i];
  }
  delete[] test_fiber_cross_harm;

  return 0;
}

void Bandwidth_Estimation::monte_carlo_integrator(){

  double temp_harmonics_0[num_harmonics];
  double temp_harmonics_1[num_harmonics];

  double ** test_point_cross_harm = new double*[num_harmonics];
  for(int i = 0; i < num_harmonics; ++i){
    test_point_cross_harm[i] = new double[num_harmonics];
  }

  double temp_kernel = 0;
  double px_term[2];
  int mc_hemi[2];
  float mc_position[6];
  for(int cv_idx = 0; cv_idx < num_test_points; ++cv_idx){

    //zero out the test_point_cross_harm
    for(int i = 0; i < num_harmonics; ++i)
      for(int j = 0; j < num_harmonics; ++j)
        test_point_cross_harm[i][j] = 0;

    //MAKE RANDOM POINT
    mc_hemi[0] = rand() % 2;
    mc_hemi[1] = rand() % 2;
    for(int i = 0; i < 6; ++i){
      mc_position[i] = (((double) rand())/RAND_MAX) - 0.5;
    }
    float mc_norm = sqrt(mc_position[0]*mc_position[0] +
      mc_position[1] * mc_position[1] + mc_position[2] * mc_position[2]);
    mc_position[0] = mc_position[0]/mc_norm;
    mc_position[1] = mc_position[1]/mc_norm;
    mc_position[2] = mc_position[2]/mc_norm;

    mc_norm = sqrt(mc_position[3]*mc_position[3] +
      mc_position[4] * mc_position[4] + mc_position[5] * mc_position[5]);
    mc_position[3] = mc_position[3]/mc_norm;
    mc_position[4] = mc_position[4]/mc_norm;
    mc_position[5] = mc_position[5]/mc_norm;

    for(int data_idx = 0; data_idx < num_tracks; ++data_idx) {

      px_term[0] = 0;
      px_term[1] = 0;

      //if it IS in the test set, skip
      if(trk_endpoints_idx[data_idx][0] == -1 ||
        trk_endpoints_idx[data_idx][2] == -1)
        continue;

      int check_flip_output = check_flip::check_flip(
        mc_hemi[0],&mc_position[0],
        mc_hemi[1],&mc_position[3],
        trk_endpoints_idx[data_idx][1],trk_endpoints_coords[data_idx],
        trk_endpoints_idx[data_idx][3],trk_endpoints_coords[data_idx] + 3);

      if(check_flip_output == 0){
        continue;
      } else if(check_flip_output == 1){
        //q
        //dot product
        for(int d = 0; d < 3; ++d){
          //x_1 * y_1
          px_term[0] += mc_position[d] * trk_endpoints_coords[data_idx][d];
          //x_2 * y_2
          px_term[1] += mc_position[3 + d] * trk_endpoints_coords[data_idx][3 + d];
        }
      } else { //otherwise we need to swap y_1 and y_2
        //reversed dot product
        for(int d = 0; d < 3; ++d){
          //x_1 * y_2
          px_term[0] += mc_position[d] * trk_endpoints_coords[data_idx][3 + d];
          //x_2 * y_1
          px_term[1] += mc_position[d + 3] * trk_endpoints_coords[data_idx][d];
        }
      }

      //round off errors? (they should be conservative but...just in case...)
      px_term[0] = std::min(px_term[0],1.0);
      px_term[1] = std::min(px_term[1],1.0);

      for(int L = 0; L < num_harmonics; ++L){
        temp_harmonics_0[L] =
          harm_lookup_table[L][(int)
            ( px_term[0]*num_harm_samps + num_harm_samps)];
        temp_harmonics_1[L] =
          harm_lookup_table[L][(int)
            ( px_term[1]*num_harm_samps + num_harm_samps)];
      }
      for(int L1 = 0; L1 < num_harmonics; ++L1)
        for(int L2 = 0; L2 < num_harmonics; ++L2)
          test_point_cross_harm[L1][L2] +=
            temp_harmonics_0[L1] * temp_harmonics_1[L2];
    }

    for(int sig_idx = 0; sig_idx < range_size; ++sig_idx){
      temp_kernel = (1.0 / num_nonnull_tracks) * 
        sigma_opt::calc_kernel(num_harmonics,
        test_point_cross_harm,
        range_min + range_step*sig_idx,
        exp_lookup[sig_idx]);
      sigma_objective[sig_idx] += (1.0/(num_test_points)) *
        temp_kernel * temp_kernel;
    }

  }
  for(int i = 0; i < num_harmonics; ++i){
    delete[] test_point_cross_harm[i];
  }
  delete[] test_point_cross_harm;

 
}

void Bandwidth_Estimation::harmonic_helper_func(
  int* idxA, float* coordsA, //test
  int* idxB, float* coordsB, //data
  double* px_term,
  double** test_fiber_cross_harm){

  //see the paper
  int check_flip_output = check_flip::check_flip(
    idxA[1],coordsA,
    idxA[3],coordsA + 3,
    idxB[1],coordsB,
    idxB[3],coordsB + 3);

  px_term[0] = 0;
  px_term[1] = 0;

  if(check_flip_output == 0){
    //if they're not in the same hemispheres...
    return;
  } else if(check_flip_output == 1){
    //otherwise if no flip

    //dot product
    for(int d = 0; d < 3; ++d){
      //x_1 * y_double1
      px_term[0] += coordsA[d] *
        coordsB[d];
      //x_2 * y_2
      px_term[1] += coordsA[3 + d] *
        coordsB[3 + d];
    }
  } else { //otherwise we need to swap y_1 and y_2
    //reversed dot product
    for(int d = 0; d < 3; ++d){
      //x_1 * y_2
      px_term[0] += coordsA[d] *
        coordsB[3 + d];
      //x_2 * y_1
      px_term[1] += coordsA[d + 3] *
        coordsB[d];
    }
  }

  //round off errors? (they should be conservative but...just in case...)
  px_term[0] = std::min(px_term[0],1.0);
  px_term[1] = std::min(px_term[1],1.0);

  for(int L = 0; L < num_harmonics; ++L){
    temp_harmonics_0[L] =
      harm_lookup_table[L][(int)
        ( px_term[0]*num_harm_samps + num_harm_samps)];
    temp_harmonics_1[L] =
      harm_lookup_table[L][(int)
        ( px_term[1]*num_harm_samps + num_harm_samps)];

  }

  for(int L1 = 0; L1 < num_harmonics; ++L1)
    for(int L2 = 0; L2 < num_harmonics; ++L2)
      test_fiber_cross_harm[L1][L2] +=
        temp_harmonics_0[L1] * temp_harmonics_1[L2];
}



//[range_size][num_harmonics][num_harmonics]
void Bandwidth_Estimation::calc_exp_lookup(){
  exp_lookup = new double**[range_size];

  double s = range_min;
  for(int i = 0; i < range_size; ++i){
    exp_lookup[i] = new double*[num_harmonics];
    for(int j = 0; j < num_harmonics; ++j){
      exp_lookup[i][j] = new double[num_harmonics];
      for(int k = 0; k < num_harmonics; ++k)
        exp_lookup[i][j][k] = exp(-s*(j*(j+1.0) + k*(k+1.0)))
          * (2.0*j + 1.0) * (2.0*k + 1.0);
          //4 * M_PI factor removed due to unit sphere normalization
    }        
    s += range_step;
  }

}
 

int Bandwidth_Estimation::save_file(std::string filename){

  std::ofstream output( filename.c_str(), std::ios::binary );

  for(int i = 0; i < range_size; ++i){
    output << sigma_objective[i] << std::endl;
  }
  output.close();

  return 0;
}

int Bandwidth_Estimation::load_file(std::string filename){

  std::cerr << "Warning: Probably can't ever use load file..." << std::endl;

  std::ifstream input( filename.c_str(), std::ios::binary );

  for(int i = 0; i < range_size; ++i){
    input >> sigma_objective[i];
  }
    
  return 0;
}
 





}//end of namespace c3



/*
      px_term[0] = 0;
      px_term[1] = 0;

      //if it IS in the test set, skip
      if(data_idx == test_idx)
        continue;
      if(trk_endpoints_idx[data_idx][0] == -1 ||
          trk_endpoints_idx[data_idx][2] == -1)
        continue;

      //this checks whether or not we should filp the coordinates,
      //see the paper
      int check_flip_output = check_flip::check_flip(
        trk_endpoints_idx[test_idx][1],trk_endpoints_coords[test_idx],
        trk_endpoints_idx[test_idx][3],trk_endpoints_coords[test_idx] + 3,
        trk_endpoints_idx[data_idx][1],trk_endpoints_coords[data_idx],
        trk_endpoints_idx[data_idx][3],trk_endpoints_coords[data_idx] + 3);

      if(check_flip_output == 0){
        //if they're not in the same hemispheres...
        continue;
      } else if(check_flip_output == 1){
        //otherwise if no flip

        //dot product
        for(int d = 0; d < 3; ++d){
          //x_1 * y_1
          px_term[0] += trk_endpoints_coords[test_idx][d] *
            trk_endpoints_coords[data_idx][d];
          //x_2 * y_2
          px_term[1] += trk_endpoints_coords[test_idx][3 + d] *
            trk_endpoints_coords[data_idx][3 + d];
        }
      } else { //otherwise we need to swap y_1 and y_2
        //reversed dot product
        for(int d = 0; d < 3; ++d){
          //x_1 * y_2
          px_term[0] += trk_endpoints_coords[test_idx][d] *
            trk_endpoints_coords[data_idx][3 + d];
          //x_2 * y_1
          px_term[1] += trk_endpoints_coords[test_idx][d + 3] *
            trk_endpoints_coords[data_idx][d];
        }
      }

      //round off errors? (they should be conservative but...just in case...)
      px_term[0] = std::min(px_term[0],1.0);
      px_term[1] = std::min(px_term[1],1.0);

      for(int L = 0; L < num_harmonics; ++L){
        temp_harmonics_0[L] =
          harm_lookup_table[L][(int)
            ( px_term[0]*num_harm_samps + num_harm_samps)];
        temp_harmonics_1[L] =
          harm_lookup_table[L][(int)
            ( px_term[1]*num_harm_samps + num_harm_samps)];

        //std::cout << test_fiber_cross_harm[cv_idx][L][0] << "\t";
        //TODO: implement an error-tracking flag
        //double error = fabs(test_fiber_cross_harm[cv_idx][L][0] -
        //  sh::EvalSH(L,0,0,acos(px_term[0])));
        //if( error > max_error){
        //  max_error = error;
        //  std::cout << "new max: " << error << std::endl;
        //  std::cout << "Harmonic: " << L << std::endl;
        //  std::cout << "px_term[0]: " << px_term[0] << std::endl;
        //  std::cout << "index: " <<
        //    (int) (px_term[0]*num_harm_samps + num_harm_samps) << std::endl;
        //}
        //sum_of_errors += error;
        //test_fiber_cross_harm[cv_idx][L][1] = sh::EvalSH(L,0,0,acos(px_term[1]));
      }

      for(int L1 = 0; L1 < num_harmonics; ++L1)
        for(int L2 = 0; L2 < num_harmonics; ++L2)
          test_fiber_cross_harm[L1][L2] +=
            temp_harmonics_0[L1] * temp_harmonics_1[L2];

 */

