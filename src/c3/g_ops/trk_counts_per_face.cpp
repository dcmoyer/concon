

#include "trk_counts_per_face.hpp"

namespace c3 {

int trk_endpoints_per_face(
  std::unordered_map< std::string, std::string> params
){

  bool verbose = false;
  std::string lh_mesh_string, rh_mesh_string;
  std::string lh_output_string, rh_output_string;
  std::string intersects_file;

  for( std::pair<std::string,std::string> p : params ){

    if(p.first == "verbose"){
      verbose = true;
    }

    if(p.first == "lh_mesh"){
      lh_mesh_string = p.second;
    }

    if(p.first == "rh_mesh"){
      rh_mesh_string = p.second;
    }

    if(p.first == "lh_output"){
      lh_output_string = p.second;
    }

    if(p.first == "rh_output"){
      rh_output_string = p.second;
    }

    if(p.first == "intersects_file"){
      intersects_file = p.second;
    }

  }

  //
  //
  //

  if(verbose){
    std::cout << "[c3/g_ops/trk_counts_per_face] "
      << "Loading meshes." << std::endl;
  }

  MeshLib::Solid lh_mesh,rh_mesh;
  lh_mesh.read( lh_mesh_string.c_str() );
  rh_mesh.read( rh_mesh_string.c_str() );

  //
  //
  //

  if(verbose){
    std::cout << "[c3/g_ops/trk_counts_per_face] "
      << "Loading intersections from " << intersects_file << std::endl;
  }

  //TODO: this needs to be a class structure.

  FILE * sphere_coord_file = fopen( intersects_file.c_str() ,"r");
  char buffer[500];

  //read file header
  int num_tracks;
  fscanf(sphere_coord_file,"#%i\n", &num_tracks);
  int num_nonnull_tracks = num_tracks;

  int ** trk_endpoints_idx = new int*[num_tracks];
  float ** trk_endpoints_coords = new float*[num_tracks];
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
      std::cerr << "[c3/g_ops/trk_counts_per_face] " <<
        " Endpoint read error, aborting." << std::endl;
      exit(1);
    }

    if(trk_endpoints_idx[i][0] == -1 ||
      trk_endpoints_idx[i][2] == -1){
      num_nonnull_tracks--;
    }
  }
  fclose(sphere_coord_file);

  if(verbose){
    std::cout << "[c3/g_ops/trk_counts_per_face] "
      << "Counting." << std::endl;
  }

  int rh_size = rh_mesh.numVertices();
  std::vector<int> mesh_counts(
    rh_mesh.numVertices() + lh_mesh.numVertices(),0);

  int triangle_idx, mesh_idx;
  for(int i = 0; i < num_tracks; ++i){

    if(trk_endpoints_idx[i][0] == -1 ||
      trk_endpoints_idx[i][2] == -1){
      continue;
    }

    triangle_idx = trk_endpoints_idx[i][0];
    mesh_idx = trk_endpoints_idx[i][1];

    //TODO: change to vertex iter
    if(mesh_idx >= 0){
      MeshLib::Face* f;
      if(mesh_idx == 0){
        f = rh_mesh.idFace(triangle_idx);
      } else {
        f = lh_mesh.idFace(triangle_idx);
      }
      for(MeshLib::FaceVertexIterator fv_it(f); !fv_it.end(); ++fv_it){
        mesh_counts[(*fv_it)->id() + mesh_idx * rh_size - 1] += 1;
      }
    }

    triangle_idx = trk_endpoints_idx[i][2];
    mesh_idx = trk_endpoints_idx[i][3];

    if(mesh_idx >= 0){
      MeshLib::Face* f;
      if(mesh_idx == 0){
        f = rh_mesh.idFace(triangle_idx);
      } else {
        f = lh_mesh.idFace(triangle_idx);
      }
      for(MeshLib::FaceVertexIterator fv_it(f); !fv_it.end(); ++fv_it){
        mesh_counts[(*fv_it)->id() + mesh_idx * rh_size - 1] += 1;
      }
    }
  }

  //
  //
  //

  int max_count = 0;
  //int total_count = 0;
  int s = mesh_counts.size();
  for(int i = 0; i < s; ++i){
    if(mesh_counts[i] > max_count){
      max_count = mesh_counts[i];
    }

    //total_count += mesh_counts[i];
  }

  if(verbose){
    std::cout << "[c3/g_ops/trk_counts_per_face] "
      << "Coloring." << std::endl;
  }

  float cold[3] = {50,50,50};
  float hot[3] = {250,50,50};
  float color[3],r;

  //RH
  int i = 0;
  for(MeshLib::SolidVertexIterator viter(&rh_mesh);
      !viter.end(); ++viter){

    r = (float) mesh_counts[i] / (float) max_count;
    for(int j = 0; j < 3; ++j){
      color[j] = r * hot[j] + (1 - r) * cold[j];
    }
    (*viter)->color_vertex(
      color[0]/255.0,
      color[1]/255.0,
      color[2]/255.0
    );

    i += 1;
  }

  //LH
  i = 0;
  for(MeshLib::SolidVertexIterator viter(&lh_mesh);
      !viter.end(); ++viter){

    r = (float) mesh_counts[i + rh_size] / (float) max_count;
    for(int j = 0; j < 3; ++j)
      color[j] = r * hot[j] + (1 - r) * cold[j];
    (*viter)->color_vertex(
      color[0]/255.0,
      color[1]/255.0,
      color[2]/255.0
    );

    i += 1;
  }

  if(verbose){
    std::cout << "[c3/g_ops/trk_counts_per_face] "
      << "Writing." << std::endl;
  }

  lh_mesh.write(lh_output_string.c_str());
  rh_mesh.write(rh_output_string.c_str());

  return 0;
}

}//end of namespace c3










