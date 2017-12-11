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



#include "c3/s_ops/effective_connectivity.hpp"

namespace c3{

Effective_Connectivity::Effective_Connectivity() : c3::Subject(){

  r_trk_endpoints = true;
  r_grid = true;
  r_harm_lookup_table = true;
  r_kern_lookup_table = true;

  output_numerator = NULL;
  output_denominator = NULL;
  N_output = -1;

  temp_harmonics_0 = new double[num_harmonics];
  temp_harmonics_1 = new double[num_harmonics];

  weights_set = false;
  trk_weights = NULL;

}

Effective_Connectivity::~Effective_Connectivity(){

  delete[] temp_harmonics_0;
  delete[] temp_harmonics_1;

  if(output_numerator != NULL){
    for(int i = 0; i < N_output; ++i)
      delete[] output_numerator[i];
    delete[] output_numerator;
  }

  if(output_denominator != NULL){
    for(int i = 0; i < N_output; ++i)
      delete[] output_denominator[i];
    delete[] output_denominator;

  }

  if(trk_weights != NULL)
    delete trk_weights;

}

int Effective_Connectivity::set(
  std::unordered_map<std::string, std::string> params){

  bool sigma_used = false;
  bool epsilon_used = false;
  bool final_thold_used = false;

  //set vars here
  //
  for( std::pair<std::string,std::string> p : params ){
    if(p.first == "sigma"){
      sigma = c3::string_to_double(p.second);
      sigma_used = true;
    }

    if(p.first == "epsilon"){
      epsilon = c3::string_to_double(p.second);
      epsilon_used = true;
    }

    if(p.first == "final_thold"){
      final_thold = c3::string_to_double(p.second);
      final_thold_used = true;
    }

    if(p.first == "verbose"){
      verbose = true;
    }

    if(p.first == "trk_weight_prefix"){
      trk_weight_prefix = p.second;
    }

    if(p.first == "trk_weight_postfix"){
      trk_weight_postfix = p.second;
    }

  }

  if(!sigma_used && !epsilon_used && !final_thold_used){
    return 1;
  }
  return 0;

}

int Effective_Connectivity::subj_specific_load(std::string subj){


  if(trk_weight_prefix.length() != 0 || trk_weight_postfix.length() != 0){
    weights_set = true;
    trk_weights = new double[num_tracks];

    //load in weights
    //TODO: FIX THIS HACK
    std::string filename =
      trk_weight_prefix + subj + "/" + subj + trk_weight_postfix;
    std::ifstream weight_in(filename.c_str());

    //check this against file config
    for(int i = 0; i < num_tracks; ++i){
      weight_in >> trk_weights[i];
    }
    weight_in.close();

  }

  return 0;
}

bool Effective_Connectivity::check_req(){

  if(trk_endpoints_idx == NULL){
    std::cerr << "Error: trk endpoint file not loaded. Aborting"
     << std::endl;
    return false;
  }
  //also checks trk_endpoints_coords

  if(kern_lookup_table == NULL){
    std::cerr << "Error: Kernel Lookup Table not set. Aborting"
     << std::endl;
    return false;
  }

  if(harm_lookup_table == NULL){
    std::cerr << "Error: Harmonic Lookup Table not set. Aborting"
     << std::endl;
    return false;
  }

  if(lh_grid == NULL || rh_grid == NULL){
    std::cerr << "Error: Grid not set. Aborting" << std::endl;
    return false;
  }

  if(lh_vertex_lookup_QT == NULL || rh_vertex_lookup_QT == NULL){
    std::cerr << "Error: Grid vertex lookup not set. Aborting."
     << std::endl;
    return false;
  }

  return true;
}
  
int Effective_Connectivity::run(){

  if(!check_req()){
    return 1;
  }

  if(verbose)
    std::cout << "Computing Cutoff Distance" << std::endl;

  double cutoff_distance = 0;
  double * cross_harm[num_harmonics];
  double temp_kernel = 0;
  for(int i = 0; i < num_harmonics; ++i)
    cross_harm[i] = new double[num_harmonics];
  for(double x = 1.0; x > -1.0; x -= 0.0001){

    for(int L = 0; L < num_harmonics; ++L){
      temp_harmonics_0[L] =
        harm_lookup_table[L][(int) ( x*num_harm_samps + num_harm_samps)];
      temp_harmonics_1[L] =
        harm_lookup_table[L][(int) ( 1*num_harm_samps + num_harm_samps)];
    }

    for(int L1 = 0; L1 < num_harmonics; ++L1)
      for(int L2 = 0; L2 < num_harmonics; ++L2)
        cross_harm[L1][L2] = temp_harmonics_0[L1] * temp_harmonics_1[L2];

    temp_kernel = sigma_opt::calc_kernel(num_harmonics, cross_harm, sigma);

    if(temp_kernel < epsilon){
      cutoff_distance = x;
      break;
    }
  }
  for(int i = 0; i < num_harmonics; ++i)
    delete cross_harm[i];

  if(verbose)
    std::cout << "cutoff distance : " << cutoff_distance << std::endl
      << "Beginning Actual Kernel Computation" << std::endl;

  double cutoff_angle = acos(cutoff_distance);

  //
  // Actual kernel computation
  //

  N_output = rh_grid->numVertices() + lh_grid->numVertices();

  output_numerator = new double * [N_output];
  for(int i = 0; i < N_output; ++i){
    output_numerator[i] = new double[N_output];
    for(int j = 0; j < N_output; ++j){
      output_numerator[i][j] = 0;
    }
  }

  output_denominator = new double * [N_output];
  for(int i = 0; i < N_output; ++i){
    output_denominator[i] = new double[N_output];
    for(int j = 0; j < N_output; ++j){
      output_denominator[i][j] = 0;
    }
  }

  int N_mesh_1 = rh_grid->numVertices();
  time_t start_kernel_calc = time(0);
  double dot_prod1, dot_prod2, scale_factor;
  double y_min,y_max;
  int hemi_1, hemi_2, s1, s2;
  double x,y;
  double trk_w;
  int num_true_mesh = 0;
  std::vector<quad_tree::Quad_Tree_Data> qtd_vector_1,qtd_vector_2;
  MeshLib::Point p1, p2, fiber_p1, fiber_p2;

  quad_tree::Quad_Tree* vertex_lookup_QT_1, * vertex_lookup_QT_2;
  MeshLib::Solid* mesh_1,* mesh_2;

  for(int fiber_idx = 0; fiber_idx < num_tracks; ++fiber_idx){

    qtd_vector_1.clear();
    qtd_vector_2.clear();

    hemi_1 = trk_endpoints_idx[fiber_idx][1]; 
    hemi_2 = trk_endpoints_idx[fiber_idx][3];

    if(hemi_1 == -1 || hemi_2 == -1)
      continue;

    //get specific track weight
    if(weights_set){
      trk_w = trk_weights[fiber_idx];
    } else {
      trk_w = 1.0;
    }

    //these vars are from the Subject class
    if(hemi_1 == 0){
      mesh_1 = rh_grid;
      vertex_lookup_QT_1 = rh_vertex_lookup_QT;
    } else { //someday, subcortical structures? cerebellum?
      mesh_1 = lh_grid;
      vertex_lookup_QT_1 = lh_vertex_lookup_QT;
    }

    if(hemi_2 == 0){
      mesh_2 = rh_grid;
      vertex_lookup_QT_2 = rh_vertex_lookup_QT;
    } else { //someday, subcortical structures? cerebellum?
      mesh_2 = lh_grid;
      vertex_lookup_QT_2 = lh_vertex_lookup_QT;
    }

    fiber_p1 = MeshLib::Point(trk_endpoints_coords[fiber_idx][0],
     trk_endpoints_coords[fiber_idx][1], trk_endpoints_coords[fiber_idx][2]);
    fiber_p2 = MeshLib::Point(trk_endpoints_coords[fiber_idx][3],
     trk_endpoints_coords[fiber_idx][4], trk_endpoints_coords[fiber_idx][5]);

    //Finds lists of mesh vertices within the cutoff boundaries
    //unrolled in the loop unrolling...
    for(int fiber_endpoint = 0; fiber_endpoint < 2; ++fiber_endpoint){

      if(fiber_endpoint == 0){
        x = fiber_p1.theta();
        y = fiber_p1.phi();
      } else {
        x = fiber_p2.theta();
        y = fiber_p2.phi();
      }
  
      //scale cylindrical projection
      if(x - cutoff_angle < 0 || x + cutoff_angle > M_PI){
        y_min = 0;
        y_max = 2*M_PI;
      } else {
        scale_factor = fabs(1.0/cos(fmax(
          fabs(M_PI/2.0 - x - cutoff_angle),
          fabs(M_PI/2.0 - x + cutoff_angle)
        )));
        y_min = y - (cutoff_angle * scale_factor);
        y_max = y + (cutoff_angle * scale_factor);
        if(y_max - y_min > 2.0*M_PI){
          y_min = 0;
          y_max = 2*M_PI;
        }
      }

      if(fiber_endpoint == 0){
        vertex_lookup_QT_1->get_range_periodic(
          fmax(x - cutoff_angle,0),
          y_min,
          fmin(x + cutoff_angle,M_PI),
          y_max,
          qtd_vector_1,
          false,true); //y is periodic, x is not
      } else {
        vertex_lookup_QT_2->get_range_periodic(
          fmax(x - cutoff_angle,0),
          y_min,
          fmin(x + cutoff_angle,M_PI),
          y_max,
          qtd_vector_2,
          false,true); //y is periodic, x is not
      }
    }

    s1 = qtd_vector_1.size();
    s2 = qtd_vector_2.size();

    for(int i = 0; i < s1; ++i){

      p1 = mesh_1->idVertex(qtd_vector_1[i].idx)->point();
      dot_prod1 = fiber_p1 * p1;
      if(dot_prod1 < cutoff_distance)
        continue;

      for(int j = 0; j < s2; ++j){

        //no self comparisions
        if(hemi_1 == hemi_2 && qtd_vector_1[i].idx == qtd_vector_2[j].idx)
          continue;

        p2 = mesh_2->idVertex(qtd_vector_2[j].idx)->point();
        dot_prod2 = fiber_p2 * p2;
        if(dot_prod2 < cutoff_distance)
          continue;

        output_numerator[qtd_vector_1[i].idx + (hemi_1 * N_mesh_1) - 1] //.m files are 1 indexed
          [qtd_vector_2[j].idx + (hemi_2 * N_mesh_1) - 1] //.m files are 1 indexed
          += 
          trk_w;
          //trk_w *
          //kern_lookup_table[(int)(dot_prod1 * num_kern_samps + num_kern_samps)] *
          //kern_lookup_table[(int)(dot_prod2 * num_kern_samps + num_kern_samps)];

        output_denominator[qtd_vector_1[i].idx + (hemi_1 * N_mesh_1) - 1] //.m files are 1 indexed
          [qtd_vector_2[j].idx + (hemi_2 * N_mesh_1) - 1] //.m files are 1 indexed
          +=
          1.0;
          //kern_lookup_table[(int)(dot_prod1 * num_kern_samps + num_kern_samps)] *
          //kern_lookup_table[(int)(dot_prod2 * num_kern_samps + num_kern_samps)];
	
      }
    }

    if(verbose){
      num_true_mesh++;
      std::cout << std::setw(40) << "\raverage time per fiber: "
        << std::setw(10) << std::setprecision(6) << std::fixed
        << ((double) (time(0) - start_kernel_calc))/ ((double) num_true_mesh)
        << std::flush;
    }

  }
  if(verbose){
    std::cout << std::setw(40) << "\raverage time per fiber: "
      << std::setw(10) << std::setprecision(6) << std::fixed
      << ((double)(time(0) - start_kernel_calc))/
           ((double) num_nonnull_tracks)
      << std::endl;
  }

  return 0;
}

int Effective_Connectivity::save_file(std::string filename){

  //FILE * output_file = fopen( argv[3] , "wb");
  std::cout << "printing to " << filename << std::endl;
  std::ofstream ofs( filename.c_str(), std::ios::binary );
  ofs.write( reinterpret_cast<char*>( &N_output ), sizeof(int));

  double temp;
  for(int i = 0; i < N_output; ++i){
    for(int j = 0; j < N_output; ++j){
      if(j==i)
        continue;

      if((output_denominator[i][j] + output_denominator[i][j]) < final_thold)
        continue;

      temp = (output_numerator[i][j] + output_numerator[j][i]) /
        (output_denominator[i][j] + output_denominator[i][j]);

      if(temp > final_thold){
        ofs.write( reinterpret_cast<char*>( &i ), sizeof(int));
        ofs.write( reinterpret_cast<char*>( &j ), sizeof(int));
        ofs.write( reinterpret_cast<char*>( &(temp) ), sizeof(double));
      }
        //fprintf(output_file,"%d%d%f",i,j,output[i][j]);
    }
  }
  //fprintf(output_file,"%d%d%f",-1,-1,-1.0);
  //fclose(output_file);
  ofs.close();
  std::cout << "finished printing to " << filename << std::endl;

  return 0;
}

int Effective_Connectivity::load_file(std::string filename){
  return 1;
}



}//end of namespace c3




