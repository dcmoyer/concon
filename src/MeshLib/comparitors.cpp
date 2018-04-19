//Copyright 2017 Yalin Wang, Boris Gutman, and Daniel Moyer
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

#include "MeshLib/comparitors.hpp"

namespace MeshLib {

bool is_equal_in_fv(
  MeshLib::Solid& left, MeshLib::Solid& right,
  bool verbose){

  //Summary Stats
  int num_vertices = left.numVertices();
  if( right.numVertices() != num_vertices){
    if(verbose)
      std::cout << "[MeshLib:comparitors] Not-equal, different numvertices.";
    return false;
  }

  int num_faces = left.numFaces();
  if( right.numFaces() != num_faces){
    if(verbose)
      std::cout << "[MeshLib:comparitors] Not-equal, different num faces.";
    return false;
  }

  if( left.numEdges() != right.numEdges() ){
    if(verbose)
      std::cout << "[MeshLib:comparitors] Not-equal, different num edges.";
    return false;
  }

  //Vertices
  {
    MeshLib::SolidVertexIterator viter_left(&left);
    MeshLib::SolidVertexIterator viter_right(&right);

    MeshLib::Point p_left, p_right;

    for(int i = 0; i < num_vertices; ++i){
      p_left = (*viter_left)->point();
      p_right = (*viter_right)->point();

      for(int j = 0; j < 3; ++j){
        if(p_left[j] != p_right[j]){
          if(verbose)
            std::cout << "[MeshLib:comparitors] Vertex coord mismatch.";
          return false;
        }
      }

      ++viter_left;
      ++viter_right;
    }
  }

  //Faces
  {
    MeshLib::SolidFaceIterator fiter_left(&left);
    MeshLib::SolidFaceIterator fiter_right(&right);

    //TODO: this might be slow? time this and/or find alternative
    std::vector<int> p_left,p_right;
    int s;

    for(int i = 0; i < num_faces; ++i){
      for( MeshLib::FaceVertexIterator fv_iter(*fiter_left);
        !fv_iter.end(); ++fv_iter){
        p_left.push_back(
          left.vertexId( *fv_iter )
        );
      }

      for( MeshLib::FaceVertexIterator fv_iter(*fiter_right);
        !fv_iter.end(); ++fv_iter){
        p_right.push_back(
          right.vertexId( *fv_iter )
        );
      }

      if(p_left.size() != p_right.size()){
        if(verbose)
          std::cout << "[MeshLib:comparitors] Face size mismatch.";
        return false;
      }

      s = p_left.size();
      //vertex equivalence guaranteed by this point (including index eq)
      for(int j = 0; j < s; ++j){
        if(p_left[j] != p_right[j]){
          if(verbose)
            std::cout << "[MeshLib:comparitors] Face mismatch.";
          return false;
        }
      }

      p_left.clear();
      p_right.clear();

      ++fiter_left;
      ++fiter_right;
    }
  }

  return true;
}

}//end of namespace MeshLib

