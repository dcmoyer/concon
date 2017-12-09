

#include "c3/g_ops/bary_to_sphere.hpp"


namespace c3 {

int bary_to_sphere(
  std::unordered_map<std::string,std::string> params,
  std::string subj //defaulted empty
){

  //
  //  Local Param Vars
  //
  int verbose = 0;
  std::string sphere_prefix, lh_sphere_postfix, rh_sphere_postfix;
  std::string bary_coords_prefix, bary_coords_postfix;
  std::string sphere_coords_prefix, sphere_coords_postfix;

  //
  //  Handler for params
  //

  for( std::pair<std::string,std::string> p : params ){

    if(p.first == "verbose"){
      verbose = c3::string_to_int(p.second);
    }

    if(p.first == "sphere_prefix"){
      sphere_prefix = p.second;
    }

    if(p.first == "rh_sphere_postfix"){
      rh_sphere_postfix = p.second;
    }

    if(p.first == "lh_sphere_postfix"){
      lh_sphere_postfix = p.second;
    }

    if(p.first == "bary_coords_prefix"){
      bary_coords_prefix = p.second;
    }

    if(p.first == "bary_coords_postfix"){
      bary_coords_postfix = p.second;
    }

    if(p.first == "sphere_coords_prefix"){
      sphere_coords_prefix = p.second;
    }

    if(p.first == "sphere_coords_postfix"){
      sphere_coords_postfix = p.second;
    }

  }

  //
  //  Load in/Local Vars
  //
  MeshLib::Solid sphere_mesh_R;
  MeshLib::Solid sphere_mesh_L;

  std::string filename;

  //TODO: file exists/accessible checks

  filename = sphere_prefix + subj + rh_sphere_postfix;
  sphere_mesh_R.read(filename.c_str());

  filename = sphere_prefix + subj + lh_sphere_postfix;
  sphere_mesh_L.read(filename.c_str());

  //bary_coords
  filename = bary_coords_prefix + subj + bary_coords_postfix;
  //TODO: replace bary_coord_file reader with standard intersect interface
  FILE *bary_coord_file = fopen( filename.c_str() ,"r");
  char buffer[500];
  int N;

  fscanf(bary_coord_file,"#%i\n", &N);
  int ** idx_mat = new int*[N];
  float ** coord_mat = new float*[N];
  for(int i = 0; i < N; ++i){
    idx_mat[i] = new int[4];
    coord_mat[i] = new float[6];
  }

  //first line burn
  fgets(buffer, 500, bary_coord_file);
  for(int i = 0; i < N; ++i){
    
    fscanf(bary_coord_file,"%i\t%i\t%f\t%f\t%f\t%i\t%i\t%f\t%f\t%f\n",
      &(idx_mat[i][0]),
      &(idx_mat[i][1]),
      &(coord_mat[i][0]),
      &(coord_mat[i][1]),
      &(coord_mat[i][2]),
      &(idx_mat[i][2]),
      &(idx_mat[i][3]),
      &(coord_mat[i][3]),
      &(coord_mat[i][4]),
      &(coord_mat[i][5]));
    //if( i % 1000 == 0){
    //  std::cout << "Loaded " << i << " lines" << std::endl;
    //}
  }
  fclose(bary_coord_file);

  //output
  filename = sphere_coords_prefix + subj + sphere_coords_postfix;

  FILE * world_coord_file = fopen( filename.c_str(), "w");
  fprintf(world_coord_file,"#%i\n",N);
  MeshLib::Solid * mesh_1, * mesh_2;
  MeshLib::Face * current_face;
  MeshLib::Point coord, new_coords;

  for(int i = 0; i < N; ++i){

    if( idx_mat[i][0] < 0 || //No Intersection
        idx_mat[i][1] < 0 ||
        idx_mat[i][2] < 0 ||
        idx_mat[i][3] < 0){
      fprintf(world_coord_file,"%i\t%i\t%f\t%f\t%f\t%i\t%i\t%f\t%f\t%f\n",
        -1,-1,-1.0,-1.0,-1.0,-1,-1,-1.0,-1.0,-1.0);
    } else {

      //FIRST END POINT
      //get appropriate mesh
      if(idx_mat[i][1] == 0){ //right
        mesh_1 = &sphere_mesh_R;
      } else {
        mesh_1 = &sphere_mesh_L;
      }

      //SECOND END POINT
      //get appropriate mesh
      if(idx_mat[i][3] == 0){ //right
        mesh_2 = &sphere_mesh_R;
      } else {
        mesh_2 = &sphere_mesh_L;
      }

      if( idx_mat[i][0] > mesh_1->numFaces() ||
          idx_mat[i][2] > mesh_2->numFaces() ){
        fprintf(world_coord_file,"%i\t%i\t%f\t%f\t%f\t%i\t%i\t%f\t%f\t%f\n",
          -1,-1,-1.0,-1.0,-1.0,-1,-1,-1.0,-1.0,-1.0);
        std::cerr << "[c3/g_ops/bary_to_sphere]: " <<
          "Bad index call, index exceeds numFaces(), aborting." << std::endl;
        exit(1);
      }
 
      //get correct face
      current_face = mesh_1->idFace(idx_mat[i][0]);

      coord = MeshLib::Point(coord_mat[i][0],coord_mat[i][1],coord_mat[i][2]);
      new_coords = current_face->get_world_coords_from_barry( coord );
      new_coords = new_coords / new_coords.norm();
      fprintf(world_coord_file,"%i\t%i\t%f\t%f\t%f\t",
        idx_mat[i][0],
        idx_mat[i][1],
        new_coords[0],new_coords[1],new_coords[2]);

      //get correct face
      current_face = mesh_2->idFace(idx_mat[i][2]);

      coord = MeshLib::Point(coord_mat[i][3],coord_mat[i][4],coord_mat[i][5]);
      new_coords = current_face->get_world_coords_from_barry( coord );
      new_coords = new_coords / new_coords.norm();
      fprintf(world_coord_file,"%i\t%i\t%f\t%f\t%f\n",
        idx_mat[i][2],
        idx_mat[i][3],
        new_coords[0],new_coords[1],new_coords[2]);

    }

  }


  fclose(world_coord_file);
  for(int i = 0; i < N; ++i){
    delete[] idx_mat[i];
    delete[] coord_mat[i];
  }

  delete[] idx_mat;
  delete[] coord_mat;

  return 0;
}

}//end of namespace c3
