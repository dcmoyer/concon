
#include "c3/s_ops/surface_intersect.hpp"

namespace c3{

int Surface_Intersect::set(
  std::unordered_map<std::string, std::string> params){

  //set vars here
  //
  for( std::pair<std::string,std::string> p : params ){

    if(p.first == "original_mesh_prefix"){
      original_mesh_prefix = p.second;
    }

    if(p.first == "original_mesh_lh_postfix"){
      original_mesh_lh_postfix = p.second;
    }

    if(p.first == "original_mesh_rh_postfix"){
      original_mesh_rh_postfix = p.second;
    }

    if(p.first == "trk_prefix"){
      trk_prefix = p.second;
    }

    if(p.first == "trk_postfix"){
      trk_postfix = p.second;
    }

    if(p.first == "extension_distance"){
      ex_dist = c3::string_to_double(p.second);
    }

    if(p.first == "output_prefix"){
      output_prefix = p.second;
    }

    if(p.first == "output_postfix"){
      output_postfix = p.second;
    }

    if(p.first == "verbose"){
      verbose = true;
    }

  }

  //TODO:checks

  return 0;
}

int Surface_Intersect::subj_specific_load( std::string subj_ ){
  subj = subj_;

  //load in mesh

  if(verbose)
    std::cout << "[c3/s_ops/surface_intersect] " << 
      "Starting Subj Specific Load" << std::endl;

  std::string filename = original_mesh_prefix + subj
    + original_mesh_rh_postfix;
  rh_orig.read(filename.c_str());

  filename = original_mesh_prefix + subj
    + original_mesh_lh_postfix;
  lh_orig.read(filename.c_str());

  max_x = -std::numeric_limits<float>::max();
  max_y = -std::numeric_limits<float>::max();
  max_z = -std::numeric_limits<float>::max();

  min_x = std::numeric_limits<float>::max();
  min_y = std::numeric_limits<float>::max();
  min_z = std::numeric_limits<float>::max();

  //find max and min range

  //TODO: multiple meshes
  MeshLib::Point p;
  for(int i = 0; i < 2; ++i){
    MeshLib::Solid* target = NULL;
    if(i == 0)
      target = &rh_orig;
    else
      target = &lh_orig;

    for(MeshLib::SolidVertexIterator viter(target); !viter.end(); ++viter){
      p = (*viter)->point();

      if(p[0] < min_x){
        min_x = p[0];
      }
      if ( p[0] > max_x ){
        max_x = p[0];
      }

      if(p[1] < min_y){
        min_y = p[1];
      }
      if ( p[1] > max_y ){
        max_y = p[1];
      }

      if(p[2] < min_z){
        min_z = p[2];
      }
      if ( p[2] > max_z ){
        max_z = p[2];
      }
    }
  }

  if(verbose)
    std::cout << "[c3/s_ops/surface_intersect] "
      << "Starting Oct Tree Load In." << std::endl;

  //construct oct tree
  OT.init(
    min_x - 1.0,max_x + 1.0,
    min_y - 1.0,max_y + 1.0,
    min_z - 1.0,max_z + 1.0,
    MIN_BIN_SIZE_CONST
  );

  int i = 0;
  //TODO: rework to templated OT
  for(MeshLib::SolidVertexIterator viter(&rh_orig); !viter.end(); ++viter){

    p = (*viter)->point();
    OT.insert(p[0],p[1],p[2],(*viter)->id());
    i += 1;
  }
  for(MeshLib::SolidVertexIterator viter(&lh_orig); !viter.end(); ++viter){

    p = (*viter)->point();
    OT.insert(p[0],p[1],p[2],-(*viter)->id());
    i += 1;
  }

  if(verbose)
    std::cout  << "[c3/s_ops/surface_intersect] "
      << "Starting Trk Load In." << std::endl;

  tracks = trk::Trackset();
  filename = trk_prefix + subj + trk_postfix;
  tracks.load_trk( filename.c_str() );

  if(verbose)
    std::cout << "[c3/s_ops/surface_intersect] "
      << "Leaving Subj Specific Load In." << std::endl;

  return 0;
}

bool Surface_Intersect::check_req(){
  return true;  //happy!
}

//returns the (absolute) index of the face with which p1->p2 intersects
//returns -1 otherwise
int Surface_Intersect::check_faces(MeshLib::Point& p1, MeshLib::Point& p2,
  MeshLib::Vertex* v, std::unordered_map<int,bool>& idx_lookup,
  MeshLib::Barry& inter_bary){

  bool crossed = false;
  int ret_val = -1;
  MeshLib::Point inter_point;
  for(MeshLib::VertexFaceIterator vf_it = MeshLib::VertexFaceIterator(v);
    !vf_it.end(); ++vf_it){

    //Removed 171114 1629
    //if( idx_lookup.find((*vf_it)->id()) == idx_lookup.end() ){
    //  continue;
    //}
    //idx_lookup[(*vf_it)->id()] = true;

    //check if p1->p2 intersects
    crossed = (*vf_it)->segment_cross2(p1,p2,inter_point);

    if(crossed){

      ret_val = (*vf_it)->id();
      inter_bary =
        (*vf_it)->find_barrycentric_coords(inter_point);

      break;
    }
  }

  return ret_val;
}

//
//  lookup helper
//    This function performs the majority of the work.
//    It collects the Oct_Tree range, searches for p1->p2 hits
//    to faces adj to vertices in the range, and then add any hits
//    to the intersections lists (inter_*)
//
void Surface_Intersect::lookup_helper(
  MeshLib::Point& p1,
  MeshLib::Point& p2,
  std::vector<oct_tree::Oct_Tree::OT_Data>& lookup_results,
  std::unordered_map<int,bool>& idx_lookup,
  std::vector<int>& inter_status,
  std::vector<int>& inter_idx_vec,
  std::vector<int>& inter_surf_id_vec,
  std::vector<MeshLib::Barry>& inter_bary_vec,
  int status
){
  int inter_idx;
  int surf_id;
  MeshLib::Barry inter_bary;
  MeshLib::Vertex* v;

  OT.get_range(
    p2[0] - SURF_INT_RANGE_CONST,p2[0] + SURF_INT_RANGE_CONST,
    p2[1] - SURF_INT_RANGE_CONST,p2[1] + SURF_INT_RANGE_CONST,
    p2[2] - SURF_INT_RANGE_CONST,p2[2] + SURF_INT_RANGE_CONST,
    lookup_results
  );

  //for v in vertices nearby
  for(oct_tree::Oct_Tree::OT_Data& data : lookup_results){
    if( data.idx >= 0 ){ //rh
      v = rh_orig.idVertex( data.idx );
      surf_id = 0;
    } else { //lh
      v = lh_orig.idVertex( -data.idx );
      surf_id = 1;
    }

    inter_idx = check_faces(p1,p2,v,idx_lookup,inter_bary);

    if(inter_idx != -1){
      inter_idx_vec.push_back(inter_idx);
      inter_surf_id_vec.push_back(surf_id);
      inter_bary_vec.push_back(inter_bary);
    }
  }

  return;
}


int Surface_Intersect::run(){

  if(!check_req()){
    //condition?
    return 1;
  }

  trk::TrackIterator end_of_the_track = tracks.get_end();
  bool first = true;
  MeshLib::Point p1,p2;
  std::unordered_map<int,bool> idx_lookup;

  //  single instance variables

  std::vector<oct_tree::Oct_Tree::OT_Data> lookup_results;
  int trk_size;

  //0 is regular, -1 is start ext, 1 is end ext, maybe we'll have others?
  std::vector<int> inter_status;
  std::vector<int> inter_idx_vec;
  std::vector<int> inter_surf_id_vec;
  std::vector<MeshLib::Barry> inter_bary_vec;

  float step_size = ex_dist / NUM_STEPS_CONST;

  trk_surf_intersect endpoints;

  //  output data
  int ti_idx = 0;

  for(trk::TrackIterator ti = tracks.get_iterator();
    ti != end_of_the_track; ++ti){

    if(verbose){
      std::cout << "\r[c3/s_ops/surface_intersect:run] " << ti_idx << std::flush;
      ti_idx += 1;
    }

    first = true;
    trk::StreamlineIterator si_end = (*ti).get_end();
    trk_size = (*ti).size();

    if(trk_size < 2 && verbose){
      std::cerr << "[c3/s_ops/surface_intersect:run] " <<
        "Warning, track has size < 2." << std::endl;
      output_data.push_back(trk_surf_intersect());
      continue;
    }

    //
    // extension stuff (Off the front)
    //

    MeshLib::Point start_pt =
      MeshLib::Point((*ti)[0].x,(*ti)[0].y,(*ti)[0].z);
    MeshLib::Point start_plus_one =
      MeshLib::Point((*ti)[1].x,(*ti)[1].y,(*ti)[1].z);

    for(int i = 1; i <= NUM_STEPS_CONST; ++i){
      p1 = start_pt + ((start_pt - start_plus_one) * step_size * i);
      p2 = start_pt + ((start_pt - start_plus_one) * step_size * (i - 1));

      lookup_helper(p1,p2,lookup_results,idx_lookup,
        inter_status, inter_idx_vec, inter_surf_id_vec,inter_bary_vec,
        -1);

      p1 = p2;
      lookup_results.clear();
      idx_lookup.clear();
    }

    //
    // regular stuff (In the middle)
    //

    for(trk::StreamlineIterator si = (*ti).get_iterator();
      si != si_end; ++si){

      if(first){
        p1 = MeshLib::Point((*si).x,(*si).y,(*si).z);
        first = false;
        continue;
      }

      p2 = MeshLib::Point((*si).x,(*si).y,(*si).z);

      lookup_helper(p1,p2,lookup_results,idx_lookup,
        inter_status, inter_idx_vec, inter_surf_id_vec,inter_bary_vec,
        0);

      p1 = p2;
      lookup_results.clear();
      idx_lookup.clear();

    }

    //
    // extension stuff (Off the end)
    //

    MeshLib::Point end_pt =
      MeshLib::Point(
        (*ti)[trk_size - 1].x,
        (*ti)[trk_size - 1].y,
        (*ti)[trk_size - 1].z);
    MeshLib::Point end_minus_one =
      MeshLib::Point(
        (*ti)[trk_size - 2].x,
        (*ti)[trk_size - 2].y,
        (*ti)[trk_size - 2].z);

    for(int i = 1; i <= NUM_STEPS_CONST; ++i){
      p1 = end_pt + (end_pt - end_minus_one) * step_size * i;
      p2 = end_pt + (end_pt - end_minus_one) * step_size * (i - 1);

      lookup_helper(p1,p2,lookup_results,idx_lookup,
        inter_status, inter_idx_vec, inter_surf_id_vec,inter_bary_vec,
        1);

      p1 = p2;
      lookup_results.clear();
      idx_lookup.clear();
    }


    //process results
    //TODO: allow endpoint filtering here (masking etc)
    if(inter_idx_vec.size() > 1){
      int iiv_size = inter_idx_vec.size() - 1;
      endpoints.face_id_1 = inter_idx_vec[0];
      endpoints.face_id_2 = inter_idx_vec[iiv_size];

      endpoints.mesh_id_1 = inter_surf_id_vec[0];
      endpoints.mesh_id_2 = inter_surf_id_vec[iiv_size];

      endpoints.bary_coord_1 = inter_bary_vec[0];
      endpoints.bary_coord_2 = inter_bary_vec[iiv_size];

      output_data.push_back(endpoints);
    } else {
      output_data.push_back(trk_surf_intersect());
    }

    inter_idx_vec.clear();
    inter_surf_id_vec.clear();
    inter_bary_vec.clear();
  }

  if(verbose)
    std::cout << std::endl;

  //std::cerr << "[c3/s_ops/surface_intersect:run] " <<
  //  "WARNING: save_file not called in run func, please verify file output."
  //  << std::endl;

  save_file(output_prefix + subj + output_postfix);

  return 0;
}

int Surface_Intersect::save_file(std::string filename){

  if(verbose)
    std::cout << "[c3/s_ops/surface_intersect:save_file] "
      << "Start of save." << std::endl;

  std::ofstream out(filename.c_str());
  std::string delim = "\t";

  out << "#" << output_data.size() << std::endl;
  out << "#face_id\tmesh_id0R1L\tbcoord\tbcoord\tbcoord\t"
    << "face_id\tmesh_id0R1L\tbcoord\tbcoord\tbcoord"
    << std::endl;

  for(trk_surf_intersect& tsi : output_data){
    out << 
      tsi.face_id_1 << delim <<
      tsi.mesh_id_1 << delim <<
      tsi.bary_coord_1[0] << delim <<
      tsi.bary_coord_1[1] << delim <<
      tsi.bary_coord_1[2] << delim <<
      tsi.face_id_2 << delim <<
      tsi.mesh_id_2 << delim <<
      tsi.bary_coord_2[0] << delim <<
      tsi.bary_coord_2[1] << delim <<
      tsi.bary_coord_2[2] <<
      std::endl;
  }

  out.close();

  if(verbose)
    std::cout << "[c3/s_ops/surface_intersect:save_file] "
      << "End of save." << std::endl;

  return 0;
}

int Surface_Intersect::load_file(std::string filename){
  return 0;
}

}//end of namespace c3

