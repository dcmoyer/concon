
#ifndef C3_SUBJ_HPP
#define C3_SUBJ_HPP

#include <fstream>
#include <float.h>
#include <string>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <memory>
#include <array>
#include <functional>
#include <unordered_map>

#include "c3/check_flip.h"
#include "sh/spherical_harmonics.h"
#include "c3/sigma_opt.hpp"
#include "MeshLib/Solid.h"
#include "MeshLib/iterators.h"
#include "utils/quad_tree.hpp"
#include "c3/cc_kernel.hpp"

//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#ifdef WIN32
#include <crtdbg.h>
#endif

//TODO:move this or delete this.
int fpeek(FILE *stream);

namespace c3 {

double string_to_double(std::string in);

int string_to_int(std::string in);


//
//  Subject Class
//
//  An interface for loading, saving, and accessing ConCon kernels.
//
//
class Subject {
public:
  const int num_harmonics = 33;
  Subject();
  virtual ~Subject();

  ///
  /// This method will set params using a dictionary like pass.
  /// It needs to be overloaded for every actual operation.
  ///
  virtual int set(
    std::unordered_map<std::string,std::string> params) = 0;

  ///
  /// This method will run the specific subject child class' operation.
  ///
  virtual int run() = 0;

  virtual int subj_specific_load(std::string subj) = 0;

  virtual bool check_req() = 0;

  virtual int save_file(std::string filename) = 0;

  virtual int load_file(std::string filename) = 0;

  int load_rh_lh(std::string rh_filename, std::string lh_filename);
  void free_rh_lh();

  bool requires_trk_endpoints();
  int load_trk_endpoints(std::string filename);
  void free_trk_endpoints();

  bool requires_grid();
  int load_grid(std::string filename, bool lookup=true, bool same=true,
    std::string rh_filename="default");
  void free_grid();

  //we're not using std::shared_ptr since it doesn't play well with the
  //2d array deconstructor

  bool requires_harm_lookup_table();
  double** calc_harm_lookup_table(int exp_num_harm_samps=6, int num_harm=32);
  void set_harm_lookup_table( double ** external_lookup,
    int exp_num_harm_samps=6, int num_harm=32);

  bool requires_kern_lookup_table();
  double* calc_kern_lookup_table(double sig,
    int exp_num_kern_samps=6);
  void set_kern_lookup_table( double * external_lookup,
    int exp_num_kern_samps=6);

  bool requires_kernel();
  int load_kernel(std::string filename);
  void free_kernel();

  bool requires_mask();
  int load_mask(std::string filename);
  void free_mask();

  //TODO: remove this!!
  virtual bool no_main_save(){ return false; }

//because we implemented our own shared pointers.
private:
  void free_harm_table();
  void free_kern_table();

protected:

  bool r_trk_endpoints;
  bool r_grid;
  bool r_harm_lookup_table;
  bool r_kern_lookup_table;
  bool r_kernel;
  bool r_mask;

  void helper_construct_QT(quad_tree::Quad_Tree* QT, MeshLib::Solid* mesh);

  bool verbose;

  int** trk_endpoints_idx;
  //TODO: Change this to double...
  float** trk_endpoints_coords;

  MeshLib::Solid* lh;
  MeshLib::Solid* rh;

  MeshLib::Solid* lh_grid;
  MeshLib::Solid* rh_grid;

//  int num_harmonics;
  int num_harm_samps;
  bool harm_lookup_owner;
  double** harm_lookup_table;

  double sigma;
  int num_kern_samps;
  bool kern_lookup_owner;
  double* kern_lookup_table;

  quad_tree::Quad_Tree* lh_vertex_lookup_QT;
  quad_tree::Quad_Tree* rh_vertex_lookup_QT;
  bool grid_are_same;

  CC_Kernel* kernel;

  //for only looking at certain regions
  int ** mask;

  int num_tracks;
  int num_nonnull_tracks;

};

}//end of namespace c3



#endif

