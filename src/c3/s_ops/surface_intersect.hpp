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


//
//  Surface Intersect Operation
//
//  This class attempts to find surface intersections between tracts and
//  meshes. It assumes that meshes are in the same space as the tracts...
//
//
//

#ifndef C3_SURFACE_INTERSECT_HPP
#define C3_SURFACE_INTERSECT_HPP

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <unordered_map>
#include "MeshLib/Point.h"
#include "c3/subject.hpp"
#include "utils/oct_tree.hpp"
#include "trk/trackset.hpp"

namespace c3 { 

struct trk_surf_intersect {
  int face_id_1;
  int face_id_2;
  int mesh_id_1; //OR1L
  int mesh_id_2;
  MeshLib::Barry bary_coord_1;
  MeshLib::Barry bary_coord_2;

  trk_surf_intersect(){
    face_id_1 = -1;
    face_id_2 = -1;
    mesh_id_1 = -1;
    mesh_id_2 = -1;
  }
};

const float SURF_INT_RANGE_CONST = 2.0;
const int NUM_STEPS_CONST = 10;
const float MIN_BIN_SIZE_CONST = 1.0;

class Surface_Intersect : public c3::Subject {
public:
  Surface_Intersect() : c3::Subject(){;} //usually blank?
  ~Surface_Intersect(){;} //usually blank as well? Unless you defined vars.

  int set(std::unordered_map<std::string, std::string> params);

  int subj_specific_load(std::string subj_);

  bool check_req();  

  int save_file(std::string filename);

  int load_file(std::string filename);

  int run();

  //TODO:remove this!!
  bool no_main_save(){ return true; }

protected:

  std::string subj;

  float ex_dist;

  //add new data here
  std::string original_mesh_prefix;
  std::string original_mesh_lh_postfix;
  std::string original_mesh_rh_postfix;

  std::string trk_prefix;
  std::string trk_postfix;

  std::string output_prefix;
  std::string output_postfix;

  MeshLib::Solid rh_orig;
  MeshLib::Solid lh_orig;

  //mesh props
  double max_x;
  double max_y;
  double max_z;

  double min_x;
  double min_y;
  double min_z;

  //POSITIVE RH
  //NEGATIVE LH
  oct_tree::Oct_Tree OT;

  trk::Trackset tracks;

  int check_faces(MeshLib::Point& p1, MeshLib::Point& p2,
    MeshLib::Vertex* v, std::unordered_map<int,bool>& idx_lookup,
    MeshLib::Barry& inter_bary
  );

  void lookup_helper(
    MeshLib::Point& p1,
    MeshLib::Point& p2,
    std::vector<oct_tree::Oct_Tree::OT_Data>& lookup_results,
    std::unordered_map<int,bool>& idx_lookup,
    std::vector<int>& intersect_status,
    std::vector<int>& inter_idx_vec,
    std::vector<int>& inter_surf_id_vec,
    std::vector<MeshLib::Barry>& inter_bary_vec,
    int status
  );

  std::vector<trk_surf_intersect> output_data;

};

}//end of namespace c3


#endif


