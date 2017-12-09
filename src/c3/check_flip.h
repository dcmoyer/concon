
#ifndef CHECK_FLIP_H
#define CHECK_FLIP_H

namespace check_flip {

//
//  This function checks whether or not to flip the order of the coordinates
//  for any set of points on a union of surfaces (meshes).
//
//  INPUT:
//    mesh_ids (integers)
//    coordinates (double arrays with length dim)
//    dim (default=3)
//  OUTPUT:
//    0 if no possible alignment
//    1 if aligned
//    2 if flip of y required
//
//  Here we output a flip of y (though I think it's reversable).
//
template <typename T>
int check_flip(
  int mesh_idx_x1, T* x1, int mesh_idx_x2, T* x2,
  int mesh_idx_y1, T* y1, int mesh_idx_y2, T* y2,
  int dim=3){

  int lock_x1 = -1;
  int lock_x2 = -1;

  //check if we're even on the correct meshes
  //if we're not, no possible alignment
  if( mesh_idx_x1 != mesh_idx_y1 &&
      mesh_idx_x1 != mesh_idx_y2 ){
    return 0;
  } else if( mesh_idx_x1 != mesh_idx_y1 ){
    lock_x1 = 2; //must_flip
  } else if(mesh_idx_x1 != mesh_idx_y2 ){
    lock_x1 = 1; //no flip
  }

  //same as above, for x2
  if( mesh_idx_x2 != mesh_idx_y1 &&
      mesh_idx_x2 != mesh_idx_y2 ){
    return 0;
  } else if( mesh_idx_x2 != mesh_idx_y1 ){
    lock_x2 = 1; //no flip
  } else if( mesh_idx_x2 != mesh_idx_y2 ){
    lock_x2 = 2; //must flip
  }

  //if both must flip
  if( lock_x1 == lock_x2 ){
    if (lock_x1 == 1){//both no flip
      return 1;
    } else if (lock_x1 == 2) {//both flip
      return 2;
    }//else do the rest of this function
  } else if( lock_x1 != -1 &&
      lock_x2 != -1 ){
    return 0; //one must flip, other cannot flip
  } else if( lock_x1 != -1 ||
      lock_x2 != -1 ){ //I think this case is not possible
    if(lock_x1 == -1){ //one must flip, other doesn't care
      return lock_x2;
    } else { 
      return lock_x1;
    }
  }

  //otherwise we're already on the right meshes
  //sanity check though
  if(lock_x1 != -1 && lock_x2){
    std::cout << "ERROR check_flip.h: bad conditions" << std::endl;
    exit(1);
  }

  //remember we're using cosine distance
  double product_x1_y1 = 0;
  double product_x1_y2 = 0;
  double product_x2_y1 = 0;
  double product_x2_y2 = 0;

  for(int i = 0; i < dim; ++i){

    product_x1_y1 += x1[i] * y1[i];
    product_x1_y2 += x1[i] * y2[i];
    product_x2_y1 += x2[i] * y1[i];
    product_x2_y2 += x2[i] * y1[i];

  }

  if(product_x1_y1 + product_x2_y2 > product_x1_y2 + product_x2_y1){
    return 1;
  } else {
    return 2;
  }

}

}

#endif



