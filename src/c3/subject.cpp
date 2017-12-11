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


#include "c3/subject.hpp"

//damnit
//TODO: reimplement everything using either c files or streams.
int fpeek(FILE *stream)
{
    int c;

    c = fgetc(stream);
    ungetc(c, stream);

    return c;
}

namespace c3 {

double string_to_double(std::string in){
  double out;
  std::stringstream ss;
  ss << in;
  ss >> out;
  return out;
}

int string_to_int(std::string in){
  int out;
  std::stringstream ss;
  ss << in;
  ss >> out;
  return out;
}

bool Subject::requires_trk_endpoints(){
  return r_trk_endpoints;
}

bool Subject::requires_grid(){
  return r_grid;
}

bool Subject::requires_harm_lookup_table(){
  return r_harm_lookup_table;
}

bool Subject::requires_kern_lookup_table(){
  return r_kern_lookup_table;
}

bool Subject::requires_kernel(){
  return r_kernel;
}

bool Subject::requires_mask(){
  return r_mask;
}



Subject::Subject(){

  r_trk_endpoints = false;
  r_grid = false;
  r_harm_lookup_table = false;
  r_kern_lookup_table = false;
  r_kernel = false;
  r_mask = false;

  trk_endpoints_idx = NULL;
  trk_endpoints_coords = NULL;

  lh = rh = NULL;
  lh_grid = rh_grid = NULL;

  //num_harmonics now global constant
  //num_harmonics = num_harm_samps = sigma = num_kern_samps = -1;
  num_harm_samps = sigma = num_kern_samps = -1;
  harm_lookup_owner = false;
  kern_lookup_owner = false;
  harm_lookup_table = NULL;
  kern_lookup_table = NULL;

  lh_vertex_lookup_QT = rh_vertex_lookup_QT = NULL;
  grid_are_same = false;

  kernel = NULL;

  mask = NULL;

  num_tracks = -1;
  num_nonnull_tracks = -1;

  verbose = false;

}

Subject::~Subject(){

  if(trk_endpoints_idx != NULL)
    free_trk_endpoints();
  if(lh != NULL)
    delete lh;
  if(rh != NULL)
    delete rh;
  if(lh_grid != NULL)
    delete lh_grid;
  if( !grid_are_same && rh_grid != NULL)
    delete rh_grid;
  if(kernel != NULL)
    free_kernel();
  if(mask != NULL)
    free_mask();

  if(harm_lookup_owner){
    free_harm_table();
  }
  if(kern_lookup_owner){
    free_kern_table();
  }

}



int Subject::load_trk_endpoints(std::string filename){

  FILE * sphere_coord_file = fopen( filename.c_str() ,"r");
  char buffer[500];

  //read file header
  fscanf(sphere_coord_file,"#%i\n", &num_tracks);
  num_nonnull_tracks = num_tracks;

  trk_endpoints_idx = new int*[num_tracks];
  trk_endpoints_coords = new float*[num_tracks];
  for(int i = 0; i < num_tracks; ++i){
    trk_endpoints_idx[i] = new int[4];
    trk_endpoints_coords[i] = new float[6];
  }

  int err;
  //line burn comments
  while(fpeek(sphere_coord_file) == '#'){
    fgets(buffer, 500, sphere_coord_file);
  }
  for(int i = 0; i < num_tracks; ++i){
    
    err = fscanf(sphere_coord_file,"%i\t%i\t%f\t%f\t%f\t%i\t%i\t%f\t%f\t%f",
      &(trk_endpoints_idx[i][0]), //triangle (useless due to sphere conversion)
      &(trk_endpoints_idx[i][1]), //surface idx
      &(trk_endpoints_coords[i][0]), //spherical coord
      &(trk_endpoints_coords[i][1]), //spherical coord
      &(trk_endpoints_coords[i][2]), //spherical coord
      &(trk_endpoints_idx[i][2]), //triangle (useless due to sphere conversion)
      &(trk_endpoints_idx[i][3]), //surface idx
      &(trk_endpoints_coords[i][3]), //spherical_coord
      &(trk_endpoints_coords[i][4]), //spherical_coord
      &(trk_endpoints_coords[i][5])); //spherical_coord

    if(err < 10){
      std::cout << "trk read error" << std::endl;
      exit(1);
    }
    if(trk_endpoints_idx[i][0] == -1 ||
      trk_endpoints_idx[i][2] == -1){
      num_nonnull_tracks--;
    }
  }
  fclose(sphere_coord_file);

  return 0;
}

void Subject::free_trk_endpoints(){
  if(trk_endpoints_idx == NULL)
    return;

  for(int i = 0; i < num_tracks; ++i){
    delete[] trk_endpoints_idx[i];
    delete[] trk_endpoints_coords[i];
  }
  delete[] trk_endpoints_idx;
  delete[] trk_endpoints_coords;
  trk_endpoints_idx = NULL;
  trk_endpoints_coords = NULL;

}

int Subject::load_rh_lh(std::string rh_filename, std::string lh_filename){

  free_rh_lh();
  rh = new MeshLib::Solid();
  rh->read(rh_filename.c_str());

  lh = new MeshLib::Solid();
  lh->read(lh_filename.c_str());  

  return 0;
}

void Subject::free_rh_lh(){
  if(rh != NULL)
    delete rh;
  if(lh != NULL)
    delete lh;
}

int Subject::load_grid(std::string filename, bool lookup, bool same,
  std::string rh_filename){

  lh_grid = new MeshLib::Solid;
  lh_grid->read(filename.c_str());
  if(same){
    rh_grid = lh_grid;
    grid_are_same = true;
  } else {
    rh_grid = new MeshLib::Solid;
    rh_grid->read(filename.c_str());
  }

  if(lookup){
    //TODO: Explain these args.
    lh_vertex_lookup_QT = new quad_tree::Quad_Tree(2.0 * M_PI, M_PI, 0, 0);
    helper_construct_QT(lh_vertex_lookup_QT, lh_grid);
    if(same){
      rh_vertex_lookup_QT = new quad_tree::Quad_Tree(2.0 * M_PI, M_PI, 0, 0);
      helper_construct_QT(rh_vertex_lookup_QT, rh_grid);
    } else {
      rh_vertex_lookup_QT = lh_vertex_lookup_QT;
    }

  }

  return 0;
}

void Subject::free_grid(){

  delete lh_grid;
  if(!grid_are_same)
    delete rh_grid;
  rh_grid = lh_grid = NULL;

  if(lh_vertex_lookup_QT != NULL){
    delete lh_vertex_lookup_QT;
    lh_vertex_lookup_QT = NULL;
  }
  if(rh_vertex_lookup_QT != NULL){
    delete rh_vertex_lookup_QT;
    rh_vertex_lookup_QT = NULL;
  }

}

double** Subject::calc_harm_lookup_table(
  int exp_num_harm_samps, int num_harm){

  num_harm_samps = pow(10, exp_num_harm_samps);
  //num_harmonics = num_harm;
  if(num_harm != 32)
    std::cout << "Warning: num_harmonics now global constant (for speed)"
     << std::endl;

  harm_lookup_table = new double*[num_harmonics];

  for(int i = 0; i < num_harmonics; ++i){
    harm_lookup_table[i] = new double[2 * num_harm_samps + 1];

    for(int j = -num_harm_samps; j < num_harm_samps + 1; ++j){
      harm_lookup_table[i][j + num_harm_samps] =
        sh::EvalSH(i,0,0,acos(((double)j)/num_harm_samps));
    }
  }
  harm_lookup_owner = true;

  return harm_lookup_table;
}

void Subject::set_harm_lookup_table(double** external_lookup,
  int exp_num_harm_samps, int num_harm){

  if(harm_lookup_owner)
    free_harm_table();
  num_harm_samps = pow(10, exp_num_harm_samps);
  //num_harmonics = num_harm;
  if(num_harm != 32)
    std::cout << "Warning: num_harmonics now global constant (for speed)"
     << std::endl;
  harm_lookup_table = external_lookup;

}


double* Subject::calc_kern_lookup_table(double sig,
  int exp_num_kern_samps){

  if(harm_lookup_table == NULL){
    std::cerr
      << "Cannot set kernel lookup without setting harmonic lookup. Aborting"
      << std::endl;
    exit(1);
  }

  sigma = sig;
  num_kern_samps = (int) pow(10, exp_num_kern_samps);
  kern_lookup_table = new double[2 * num_kern_samps + 1];
  for(int i = -num_kern_samps; i < num_kern_samps + 1; ++i){
    kern_lookup_table[num_kern_samps + i] = sigma_opt::calc_half_kernel(
      ((double) i / num_kern_samps), num_harm_samps,
      num_harmonics, harm_lookup_table, sigma);
  }
  kern_lookup_owner = true;

  return kern_lookup_table;
}


void Subject::set_kern_lookup_table(double* external_lookup,
  int exp_num_kern_samps){

  if(kern_lookup_owner)
    free_kern_table();
  num_kern_samps = pow(10, exp_num_kern_samps);
  kern_lookup_table = external_lookup;

}

void Subject::free_harm_table(){

  //dont think this can happen, but safety first?
  if(harm_lookup_table != NULL){
    for(int i = 0; i < num_harmonics; ++i)
      delete[] harm_lookup_table[i];
    delete[] harm_lookup_table;
  }
  harm_lookup_table = NULL;
  harm_lookup_owner = false;
  num_harm_samps = -1;
  //num_harmonics now a global constant
  //num_harm_samps = num_harmonics = -1;

}


void Subject::free_kern_table(){

  //dont think this can happen, but safety first?
  if(kern_lookup_table != NULL){
    delete[] kern_lookup_table;
  }
  kern_lookup_table = NULL;
  kern_lookup_owner = false;
  sigma = num_kern_samps = -1;

}

//TODO: Safety checks?
int Subject::load_kernel(std::string filename){

  kernel = new CC_Kernel;
  return kernel->load_file(filename);

}

void Subject::free_kernel(){

  if(kernel != NULL){
    delete kernel;
    kernel = NULL;
  }

}


int Subject::load_mask(std::string filename){
  //TODO this...
  return 0;
}

void Subject::free_mask(){
  return;//TODO this...
}

void Subject::helper_construct_QT(quad_tree::Quad_Tree* QT,
  MeshLib::Solid* mesh){
  
  double x,y;
  quad_tree::Quad_Tree_Data QTD;
  MeshLib::Point p;
  
  for(MeshLib::SolidVertexIterator viter(mesh); !viter.end(); ++viter)
	{
		MeshLib::Vertex *v = *viter;
    p = v->point();
    x = p.theta();
    y = p.phi();

    QTD.x = x;
    QTD.y = y;
    QTD.idx = v->id();
    QT->insert(x,y, QTD );
  }

}

}//end of namespace c3

